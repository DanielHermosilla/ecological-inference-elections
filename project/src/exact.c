#include "exact.h"
#include "globals.h"
#include "matrixUtils.h"
#include "memoizationUtil.h"
#include <cblas.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])
/*
typedef struct
{
    size_t n;
    size_t k;
    size_t *data;
} gsl_combination;
*/
typedef struct
{
    size_t **data; // Array of valid combinations
    size_t size;   // Number of valid combinations
} CombinationSet;

typedef struct
{
    CombinationSet **H; // Precomputed H_bf sets
    CombinationSet **K; // Precomputed K_bf sets
} PrecomputedSets;

PrecomputedSets PRECOMPUTED_SETS;
SizeTMatrix *CANDIDATEARRAYS = NULL;

/*
gsl_combination *gsl_combination_calloc(size_t n, size_t k);
size_t *gsl_combination_data(const gsl_combination *c);
*/
// Funciones utiles;
// gsl_sl_fact(unsigned int n) Factorial
//
// Combinatoria: Está en el siguiente structure
// typedef struct
// {
// size_t n;
// size_t k;
// size_t *data;
//} gsl_combination;

/*
gsl_combination *gsl_combination_calloc(size_t n, size_t k)

    This function allocates memory for a new combination with parameters n, k and initializes it to the
lexicographically first combination. A null pointer is returned if insufficient memory is available to create the
combination.
*/

// void gsl_combination_free(gsl_combination *c);

void precomputeH()
{
    PRECOMPUTED_SETS.H = malloc(TOTAL_BALLOTS * sizeof(CombinationSet *));
#pragma omp parallel for schedule(dynamic)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        PRECOMPUTED_SETS.H[b] = malloc(TOTAL_GROUPS * sizeof(CombinationSet));
        for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
        {                                     // For all groups, given a arbitrary ballot box
            int wbf = MATRIX_AT_PTR(W, b, f); // Total votes for group f

            // Initialize the GSL combination generator
            gsl_combination *Hcomb = gsl_combination_calloc(
                wbf + TOTAL_CANDIDATES - 1,
                TOTAL_CANDIDATES - 1); // Generate the combinations, note that the last candidate vote is redundant.
                                       // That would be used for efficiency purpose.

            size_t **combinations = NULL; // Array of pointers to combinations
            size_t count = 0;             // Counter for the number of valid combinations

            do
            {
                combinations = realloc(combinations, (count + 1) * sizeof(size_t *));
                combinations[count] = malloc(TOTAL_CANDIDATES * sizeof(size_t));

                const size_t *elements = gsl_combination_data(Hcomb);

                // Fill the current combination
                for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
                {
                    combinations[count][i] =
                        (i < TOTAL_CANDIDATES - 1)
                            ? elements[i]
                            : (wbf -
                               elements[TOTAL_CANDIDATES - 2]); // The redudancy of using TOTAL_CANDIDATES-1 votes.
                }

                count++;
            } while (gsl_combination_next(Hcomb) == GSL_SUCCESS);

            gsl_combination_free(Hcomb);

            // Store in the precomputed set
            PRECOMPUTED_SETS.H[b][f].data = combinations;
            PRECOMPUTED_SETS.H[b][f].size = count; // Total valid combinations
        }
    }
    printf("Precomputation of H complete.\n");
}

void precomputeK()
{
    PRECOMPUTED_SETS.K = malloc(TOTAL_BALLOTS * sizeof(CombinationSet *));
#pragma omp parallel for schedule(dynamic)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        PRECOMPUTED_SETS.K[b] = malloc(TOTAL_GROUPS * sizeof(CombinationSet));
        for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
        {
            int cumulativeVotes = 0;
            for (uint16_t f_prime = 0; f_prime <= f; f_prime++)
            {
                cumulativeVotes += MATRIX_AT_PTR(W, b, f_prime);
            }

            gsl_combination *Kcomb =
                gsl_combination_calloc(cumulativeVotes + TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);

            size_t **combinations = NULL;
            size_t count = 0;

            do
            {
                const size_t *elements = gsl_combination_data(Kcomb);
                int valid = 1;

                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    size_t kc = (c < TOTAL_CANDIDATES - 1)
                                    ? elements[c]
                                    : (cumulativeVotes -
                                       elements[TOTAL_CANDIDATES -
                                                2]); // If the votes for a candidate is bigger than its actual votes
                    if (kc > MATRIX_AT_PTR(X, c, b))
                    {
                        valid = 0;
                        break;
                    }
                }

                if (valid)
                {
                    combinations = realloc(combinations, (count + 1) * sizeof(size_t *));
                    combinations[count] = malloc(TOTAL_CANDIDATES * sizeof(size_t));

                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    {
                        combinations[count][c] = (c < TOTAL_CANDIDATES - 1)
                                                     ? elements[c]
                                                     : (cumulativeVotes - elements[TOTAL_CANDIDATES - 2]);
                    }
                    count++;
                }
            } while (gsl_combination_next(Kcomb) == GSL_SUCCESS);

            gsl_combination_free(Kcomb);

            PRECOMPUTED_SETS.K[b][f].data = combinations;
            PRECOMPUTED_SETS.K[b][f].size = count;
        }
    }
    printf("Precomputation of K complete.\n");
}

void freePrecomputedSets()
{
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
        {
            for (size_t i = 0; i < PRECOMPUTED_SETS.H[b][f].size; i++)
                free(PRECOMPUTED_SETS.H[b][f].data[i]);
            free(PRECOMPUTED_SETS.H[b][f].data);

            for (size_t i = 0; i < PRECOMPUTED_SETS.K[b][f].size; i++)
                free(PRECOMPUTED_SETS.K[b][f].data[i]);
            free(PRECOMPUTED_SETS.K[b][f].data);
        }
        free(PRECOMPUTED_SETS.H[b]);
        free(PRECOMPUTED_SETS.K[b]);
    }
    free(PRECOMPUTED_SETS.H);
    free(PRECOMPUTED_SETS.K);
}

//    This function frees all the memory used by the combination c.

/**
 * @brief Converts a type `double` array to a `size_t`
 *
 * Given an array that is made with `double` elements, it fills another one with `size_t` elements. Note that the
 * decimal values will be truncated.
 *
 * @param[in] *doubleVector A pointer towards the vector that has the values in `double`.
 * @param[in] *sizeTVector A pointer toward the vector with values of type `size_t`.
 * @param[in] size The size of the array
 *
 * @return void
 *
 */
void convertDoubleToSizeT(const double *doubleVector, size_t *sizeTVector, int size)
{
    for (int i = 0; i < size; i++)
    {
        sizeTVector[i] = (size_t)doubleVector[i];
    }
}

double factorial(double n)
{
    double result = 1.0;

#pragma omp parallel
    {
        double local_result = 1.0;

#pragma omp for schedule(static)
        for (int i = 1; i <= (int)n; i++)
        {
            local_result *= i;
        }

// Combine the results safely
#pragma omp atomic
        result *= local_result;
    }

    return result;
}
/**
 * @brief Calculate the chained product between probabilities as defined in `a`:
 *
 * Given an `H` element and the `f` index it computes the chained product.
 *
 * @param[in] *hElement A pointer towards the element of `H`. This element is an array with the possible combinations
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] f The index of `f`
 *
 * @return double: The result of the product
 *
 */

double prod(const size_t *hElement, const Matrix *probabilities, const int f)
{
    double log_result = 0;
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        double prob = MATRIX_AT_PTR(probabilities, f, c);
        if (prob == 0.0 && hElement[c] > 0)
        {
            return 0.0; // Early exit if probability is zero for a non-zero count
        }
        if (prob > 0)
        {
            log_result += hElement[c] * log(prob);
        }
    }
    return exp(log_result); // Exponentiate the final result
}
/**
 * @brief Calculate `a` given the multinomial coefficient, subject to the `H` set:
 *
 * Given an `H` element, it computes its multinomial coefficient subject to the total amount of votes a group received.
 * It yields the following calculation:
 * $$\binom{w_{bf}}{h_1,\dots,h_C}=\frac{w_{bf}!}{h_1!\cdots h_{C}!}$$
 *
 * @param [in] b The index of the ballot.
 * @param [in] f The index of the group
 * @param[in] *hElement A pointer towards the element of `H`. This element is an array with the possible combinations
 *
 * @return double: The result of the multinomial coefficient
 *
 */
double multinomialCoeff(const int b, const int f, const size_t *hElement)
{
    double result = gsl_sf_lngamma((int)MATRIX_AT_PTR(W, b, f) + 1); // ln(w_bf!)
    printf("\nThe result of ln factorial is\t%1.f\n", gsl_sf_lngamma(MATRIX_AT_PTR(W, b, f)));
    // double result = factorial(MATRIX_AT_PTR(W, b, f)); // Calculate w_bf! GSL_SF_FACT HAVE A MAXIMUM VALUE OF
    // 170.
    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    {
        printf("\nThe hElement is\t%1.zu\n", hElement[i]);
        result -= gsl_sf_lngamma(hElement[i] + 1); // Divide by each h_i!
    }

    return exp(result);
}

/**
 * @brief Checks if all elements of `K` are bigger than the elements from `H` to ensure non negativity.
 *
 * Given vectors `H` and `K` it iterates between all its values to see the sign of the difference
 *
 * @param[in] *hElement A pointer to the array `H`
 * @param[in] *kElement A pointer to the array `K`
 *
 * @return Boolean value that tells if the condition is met.
 *
 */

bool ifAllElements(const size_t *hElement, const size_t *kElement)
{
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        if (hElement[c] > kElement[c])
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Computes the `A` term from the pseudocode
 *
 * Given the vector `H` and the index from the ballot box and cummulative group, it computes the value from `A`, defined
 * as a product of inner products and binomial coefficients.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the f group (starting at 0).
 * @param[in] *hElement Pointer to the `H` vector
 * @param[in] *probabilities Pointer to the probabilities matrix.
 *
 * @return The value of `A`
 *
 */

double computeA(const int b, const int f, const size_t *hElement, const Matrix *probabilities)
{
    printf("\nWhen passing to computeA, the multinomial coefficient is %1.f and the product is %1.f",
           multinomialCoeff(b, f, hElement), prod(hElement, probabilities, f));
    return multinomialCoeff(b, f, hElement) * prod(hElement, probabilities, f);
}

/**
 * @brief Calculate a given `H` set.
 *
 * Given the index of the set, it returns the `H` set as a gsl_combination struct. Note that it could be used with the
 * gsl_combination functions since it doesn't have restrictions.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the f group (starting at 0).
 *
 * @return gsl_combination *: A pointer towards the struct of a gsl_combination.
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */
gsl_combination *getH(int b, int f)
{
    if (MATRIX_AT_PTR(W, b, f) < TOTAL_CANDIDATES)
    {
        gsl_combination *result1 = gsl_combination_calloc((int)TOTAL_CANDIDATES, (int)TOTAL_CANDIDATES);
        return result1;
    }
    gsl_combination *result = gsl_combination_calloc((int)MATRIX_AT_PTR(W, b, f), (int)TOTAL_CANDIDATES); // Order: n, k
    return result;
}

/**
 * @brief Tells if a combination is valid according the `K` definition
 *
 * Given the ballot box identificator, it checks if the condition for joining `K` is met.
 *
 * @param[in] *hVector
 * @param[in] b Index that represents the ballot box.
 *
 * @return A boolean that tells if the condition is met
 *
 *
 */

bool filterCombinations(const size_t *hVector, const int b)
{
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        if (hVector[c] > (size_t)MATRIX_AT_PTR(X, b, c))
        {
            return false;
        }
    }
    return true;
}

/**
 * @brief Calculate the difference between two vectors.
 *
 * Utility function to get the difference between two 1-dimensional arrays.
 *
 * @param[in] *K A pointer to the first array
 * @param[in] *H A pointer to the second array
 * @param[in, out] *arr The result of the operation
 *
 * @return void
 *
 */

void vectorDiff(const size_t *K, const size_t *H, size_t *arr)
{
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        arr[c] = K[c] - H[c];
    }
}

/**
 * @brief Calculate the value of `u_{bfgc}(k)`.
 *
 * Given the index from `u` it returns its value. If it hasn't been calculated yet, it starts iterating the recursion
 * for getting a final answer.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the first `f` groups (starting at 0).
 * @param[in] g Index that represents the `g` group (starting at 0).
 * @param[in] c Index that represents the `c` candidate (starting at 0).
 * @param[in] *vector A pointer to the vector that functions as a key.
 * @param[in, out] *memo A pointer to the hash table.
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 *
 * @return double: The value from `u_{bfgc}(k)`
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */

void printSizeTVector(size_t *vector, int size)
{
    if (!vector)
    {
        printf("Vector is NULL\n");
        return;
    }

    printf("Vector (size = %d): [", size);
    for (int i = 0; i < size; i++)
    {
        printf("%zu", vector[i]);
        if (i < size - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
}

bool checkNull(const size_t *vector, size_t size)
{
    for (size_t i = 0; i < size; i++)
    {
        if (vector[i] != 0)
        {
            return false;
        }
    }
    return true;
}

void recursion2(MemoizationTable *memo, const Matrix *probabilities)
{
#pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t f = 1; f < TOTAL_GROUPS; f++)
        {
            for (size_t k = 0; k < PRECOMPUTED_SETS.K[b][f].size; k++)
            {
                size_t *currentK = PRECOMPUTED_SETS.K[b][f].data[k];
                printf("\nThe current combination vector is:\n");
                printSizeTVector(currentK, TOTAL_CANDIDATES);
                printf("\n");
                for (size_t k2 = 0; k2 < PRECOMPUTED_SETS.K[b][f].size; k2++)
                {
                    size_t *currentK2 = PRECOMPUTED_SETS.K[b][f].data[k2];
                    printf("\nThe current combination2 vector is:\n");
                    printSizeTVector(currentK2, TOTAL_CANDIDATES);
                    printf("\n");
                    if (!ifAllElements(currentK, currentK2)) // Maybe could be optimized
                        continue;
                    // ---- Recursion: compute the value of `a` ----
                    double a = computeA(b, f, currentK2, probabilities);
                    printf("\nThe value a is\t%1.f", a);
                    double valueBefore[TOTAL_CANDIDATES][TOTAL_GROUPS];
                    if (f == 1 && checkNull(currentK2, (size_t)TOTAL_CANDIDATES))
                    {
                        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                        {
                            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                            {
                                valueBefore[c][g] = 1.0;
                            }
                        }
                    }
                    else if (f == 1 && !checkNull(currentK2, (size_t)TOTAL_CANDIDATES))
                    {
                        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                        {
                            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                            {
                                valueBefore[c][g] = 0;
                            }
                        }
                    }
                    else
                    {
                        size_t *substractionVector = malloc((size_t)TOTAL_CANDIDATES * sizeof(size_t));
                        vectorDiff(currentK, currentK2, substractionVector);
                        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                        {
                            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                            {
                                valueBefore[c][g] =
                                    getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                            }
                        }
                        free(substractionVector);
                    }
                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    {
                        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                        {
                            double valueNow = getMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES);
                            printf("\nThe value retrieved when searching is\t%1.f\n", valueNow);
                            if (valueNow == INVALID)
                            {
                                valueNow = 0.0;
                            }
                            if (f == g)
                            {
                                valueNow += valueBefore[c][g] * a * currentK2[c] /
                                            (MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f));
                            }
                            else
                            {
                                valueNow += valueBefore[c][g] * a;
                            }
#pragma omp critical
                            {
                                setMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES, valueNow);
                                printf("\nThe value to add is:\t%1.f\n", valueNow);
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Calculate the value of `q_{bgc}`.
 *
 * It calculates all of the values for q_{bgc} by the definition on the paper. It returns the array of type
 * `double`.
 *
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 *
 * @return *double: A pointer toward the array.
 *
 * @note: A single pointer is used to store the array continously. This is for using cBLAS operations later.
 *
 */
double *computeQExact(const Matrix *probabilities)
{

    // ---- Initialize CANDIDATEARRAYS, which corresponds as a copy of `X` but in size_t.
    if (CANDIDATEARRAYS == NULL)
    {
        // ---- Define the memory for the matrix and its parameters
        CANDIDATEARRAYS = malloc(sizeof(SizeTMatrix));
        CANDIDATEARRAYS->cols = X->cols;
        CANDIDATEARRAYS->rows = X->rows;
        CANDIDATEARRAYS->data = malloc(TOTAL_CANDIDATES * TOTAL_BALLOTS * sizeof(size_t));
// ---- Parallelization will consider the two inner loops as only one
#pragma omp parallel for collapse(2)
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // ---- For each ballot box
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                MATRIX_AT_PTR(CANDIDATEARRAYS, c, b) = (size_t)MATRIX_AT_PTR(X, c, b);
            }
        }
    }
    // ---- Initialize the `H` set, for recicling its values
    if (PRECOMPUTED_SETS.H == NULL && PRECOMPUTED_SETS.K == NULL)
    {
        printf("Starting the heavy precomputation\n");
        precomputeH();
        precomputeK();
    }

    printf("Starting the recursion\n");
    // ---- Initialize the hash table
    MemoizationTable *table = initMemo();
    // ---- Initialize the array to return
    double *array2 = (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double));
    recursion2(table, probabilities);

// ---- The parallelization would be made on the outer loops. It might cause a problem if the recursion tries to
// write at the same location, to be checked.
#pragma omp parallel for schedule(dynamic) collapse(2)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot box
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // ---- Initialize the denominator variable to avoid multiple computations
            double den = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
              // ---- Add the values of the denominator
                den += getMemoValue(table, b, TOTAL_GROUPS, g, c, &MATRIX_AT_PTR(CANDIDATEARRAYS, c, b),
                                    TOTAL_CANDIDATES) *
                       MATRIX_AT_PTR(probabilities, g, c);
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                // ---- Compute the numerator; the recursion shouldn't be ran since the denominator initialized the
                // hashed values.
                double num = getMemoValue(table, b, TOTAL_GROUPS, g, c, &MATRIX_AT_PTR(CANDIDATEARRAYS, c, b),
                                          TOTAL_CANDIDATES) *
                             MATRIX_AT_PTR(probabilities, g, c);
                // ---- Store the resulting values as q_{bgc}
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = num / den;
                printf("The result stored is %.1f\n", num / den);
            }
        }
    }
    freeMemo(table);
    return array2;
}

__attribute__((destructor)) void cleanUp()
{
    // ---- Destroy the candidate array of size_t
    if (CANDIDATEARRAYS != NULL)
    {
        free(CANDIDATEARRAYS->data);
        free(CANDIDATEARRAYS);
        CANDIDATEARRAYS = NULL;
    }
    // ---- Destroy the array of precomputed H
    freePrecomputedSets();
}
