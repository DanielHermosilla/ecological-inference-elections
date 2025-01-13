#include "exact.h"
#include "globals.h"
#include "matrixUtils.h"
#include "memoizationUtil.h"
#include <cblas.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>

#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])

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

typedef struct
{
    uint32_t b;
    uint16_t g;
    size_t **data;
    size_t size;
} Set;

Set **HSETS = NULL; // Global pointer to store all Hsets
Set **KSETS = NULL; // Global pointer to store all Ksets

PrecomputedSets PRECOMPUTED_SETS;
SizeTMatrix *CANDIDATEARRAYS = NULL;

// Recursive function to distribute votes
void generateConfigurations(int b, size_t *votes, int position, int remainingVotes, int numCandidates,
                            size_t ***results, size_t *count)
{
    if (position == numCandidates - 1)
    {
        // Assign remaining votes to the last candidate
        votes[position] = remainingVotes;

        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // Case where the last candidate exceeds the maximum possible amount.
            return;
        }

        // Store the result
        (*results) = realloc(*results, (*count + 1) * sizeof(size_t *));
        (*results)[*count] = malloc(numCandidates * sizeof(size_t));
        memcpy((*results)[*count], votes, numCandidates * sizeof(size_t));
        (*count)++;
        return;
    }

    for (int i = 0; i <= remainingVotes; i++)
    {
        votes[position] = i;

        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // This would imply that all of the subsequent elements aren't useful.
            return;
        }
        generateConfigurations(b, votes, position + 1, remainingVotes - i, numCandidates, results, count);
    }
}

// Count is a pointer that indicates the total amount of elements
size_t **generateAllConfigurations(int b, int totalVotes, int numCandidates, size_t *count)
{
    size_t **results = NULL;
    *count = 0;

    size_t *votes = malloc(numCandidates * sizeof(size_t));
    // Votes is the
    generateConfigurations(b, votes, 0, totalVotes, numCandidates, &results, count);

    free(votes);
    return results;
}

void generateSets()
{
    HSETS = malloc(TOTAL_BALLOTS * sizeof(Set *));
    KSETS = malloc(TOTAL_BALLOTS * sizeof(Set *));

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        HSETS[b] = malloc(TOTAL_GROUPS * sizeof(Set));
        KSETS[b] = malloc(TOTAL_GROUPS * sizeof(Set));
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            // Initialize H set
            HSETS[b][g].b = b;
            HSETS[b][g].g = g;
            // Initialize K set
            KSETS[b][g].b = b;
            KSETS[b][g].g = g;
            // Parameters for the function
            size_t total = (size_t)MATRIX_AT_PTR(W, b, g);
            // Generate configurations
            size_t count = 0;
            size_t **configurations = generateAllConfigurations(b, total, TOTAL_CANDIDATES, &count);

            // Store configurations and size
            HSETS[b][g].data = configurations;
            HSETS[b][g].size = count;

            if (g == 0)
            {
                KSETS[b][g].data = configurations;
                KSETS[b][g].size = count;
            }
            else
            {
                for (uint16_t f = 0; f < g; f++)
                {
                    total += (size_t)(MATRIX_AT_PTR(W, b, f));
                }
                size_t Kcount = 0;
                size_t **Kconfigurations = generateAllConfigurations(b, total, TOTAL_CANDIDATES, Kcount);
                KSETS[b][g].data = Kconfigurations;
                KSETS[b][g].size = Kcount;
            }
        }
    }
}

void precomputeH()
{
    PRECOMPUTED_SETS.H = malloc(TOTAL_BALLOTS * sizeof(CombinationSet *));

#pragma omp parallel for schedule(dynamic)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For all ballot boxes
        PRECOMPUTED_SETS.H[b] = malloc(TOTAL_GROUPS * sizeof(CombinationSet));

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For all groups, given a arbitrary ballot box

            int wbg = (int)MATRIX_AT_PTR(W, b, g); // Total votes for group f
            size_t **combinations = NULL;          // Array of pointers to combinations
            size_t count = 0;                      // Counter for the number of valid combinations

            // We know the combinations aren't defined if wbg+1<TOTAL_CANDIDATES.
            if (wbg + 1 < TOTAL_CANDIDATES)
            {
                PRECOMPUTED_SETS.H[b][g].data = combinations; // Add the arrays of the while-loop
                PRECOMPUTED_SETS.H[b][g].size = count;        // Total valid combinations
                continue;                                     // Go to the next group in case of an invalid combination.
            }

            // Initialize the `gsl_combination` struct with all possible combinations -> it's still needed to make a
            // verification given the total amount of votes the candidates had.
            gsl_permutation *Hcomb = gsl_permutation_alloc(TOTAL_CANDIDATES); //

            do
            { // ---- Loop over all possible combinations.
                bool valid = true;
                size_t *elements = gsl_combination_data(Hcomb); // An array with combinations.
                size_t sumOfCombinations = 0;

                // Check if the combination is possible
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    sumOfCombinations += elements[c];
                    if (elements[c] >
                        MATRIX_AT_PTR(X, c, b)) // If a candidate has more votes than its real amount of votes.
                        valid = false;
                }

                if (sumOfCombinations != (size_t)wbg)
                { // If the sum of votes it's not equivalent to the total amount of votes from the group.
                    valid = false;
                }

                // On case the combination is possible; add memory for the new combination.
                if (valid)
                {
                    // Memory to add the new array
                    combinations = realloc(combinations, (count + 1) * sizeof(size_t *));
                    // Memory to add each element of the new array
                    combinations[count] = malloc(TOTAL_CANDIDATES * sizeof(size_t));
                    memcpy(combinations[count], elements, TOTAL_CANDIDATES * sizeof(size_t));
                    count++;
                }
            } while (gsl_combination_next(Hcomb) == GSL_SUCCESS);

            // Free the combination for H_bf, it's not needed anymore
            gsl_combination_free(Hcomb);

            // Store in the precomputed set
            PRECOMPUTED_SETS.H[b][g].data = combinations; // Add the arrays of the while-loop
            PRECOMPUTED_SETS.H[b][g].size = count;        // Total valid combinations
        }
    }
    printf("Precomputation of H complete.\n");
}

void precomputeK()
{
    PRECOMPUTED_SETS.K = malloc(TOTAL_BALLOTS * sizeof(CombinationSet *));
    // #pragma omp parallel for schedule(dynamic)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        PRECOMPUTED_SETS.K[b] = malloc(TOTAL_GROUPS * sizeof(CombinationSet));
        // The combinations for the 0th group is equivalent to the H_b0, it MUST be created before.
        PRECOMPUTED_SETS.K[b][0].data = PRECOMPUTED_SETS.H[b][0].data;
        PRECOMPUTED_SETS.K[b][0].size = PRECOMPUTED_SETS.H[b][0].size;

        for (uint16_t g = 1; g < TOTAL_GROUPS; g++)
        {
            size_t **combinations = NULL; // Array of pointers to combinations
            size_t count = 0;             // Counter for the number of valid combinations
            int wbf = 0;
            for (uint16_t i = 0; i < g; i++)
            {
                wbf += (int)MATRIX_AT_PTR(W, b, g); // Total votes for the first `f` groups.
            }
            // We know the combinations aren't defined if wbg+1>TOTAL_CANDIDATES.
            printf("the wbf value is %d\n", wbf);
            if (wbf + 1 < TOTAL_CANDIDATES)
            {
                PRECOMPUTED_SETS.K[b][g].data = combinations; // Add the arrays of the while-loop
                PRECOMPUTED_SETS.K[b][g].size = count;        // Total valid combinations
                continue;                                     // Go to the next group in case of an invalid combination.
            }

            gsl_combination *Kcomb = gsl_combination_calloc(wbf + 1, TOTAL_CANDIDATES); //
            do
            {
                bool valid = true;
                size_t *elements = gsl_combination_data(Kcomb); // An array with combinations.
                size_t sumOfCombinations = 0;

                // Check if the combination is possible
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    sumOfCombinations += elements[c];
                    if (elements[c] >
                        MATRIX_AT_PTR(X, c, b)) // If a candidate has more votes than its real amount of votes.
                        valid = false;
                }

                if (sumOfCombinations != (size_t)wbf)
                { // If the sum of votes it's not equivalent to the total amount of votes from the group.
                    valid = false;
                }

                // On case the combination is possible; add memory for the new combination.
                if (valid)
                {
                    printf("Adding the combination\n{");
                    gsl_combination_fprintf(stdout, Kcomb, " %u");
                    printf(" }\n");
                    // Memory to add the new array
                    combinations = realloc(combinations, (count + 1) * sizeof(size_t *));
                    // Memory to add each element of the new array
                    combinations[count] = malloc(TOTAL_CANDIDATES * sizeof(size_t));
                    memcpy(combinations[count], elements, TOTAL_CANDIDATES * sizeof(size_t));
                    count++;
                }
            } while (gsl_combination_next(Kcomb) == GSL_SUCCESS);

            gsl_combination_free(Kcomb);

            PRECOMPUTED_SETS.K[b][g].data = combinations;
            PRECOMPUTED_SETS.K[b][g].size = count;
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
    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    {
        result -= gsl_sf_lngamma(hElement[i] + 1); // Divide by each h_i!
    }

    return exp(result);
}

/**
 * @brief Checks if all elements of `H` are bigger than the elements from `K` to ensure non negativity.
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
    printf("\nWhen passing to computeA, the multinomial coefficient is %1.f and the product is %.7f",
           multinomialCoeff(b, f, hElement), prod(hElement, probabilities, f));
    return multinomialCoeff(b, f, hElement) * prod(hElement, probabilities, f);
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
    printf("entering...");
    // #pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t f = 1; f < TOTAL_GROUPS; f++)
        {
            printf("THE FIRST F VALUE IS %d\n", f);
            printf("PRECOMPUTED_SETS.K[b][%d].size: %zu\n", f, PRECOMPUTED_SETS.K[b][f].size);
            sleep(5);

            for (size_t k = 0; k < PRECOMPUTED_SETS.K[b][f].size; k++)
            {
                printf("PRECOMPUTED_SETS.K[b][f].data: %p, PRECOMPUTED_SETS.K[b][f].data[k]: %p, f = %d\n",
                       PRECOMPUTED_SETS.K[b][f].data,
                       k < PRECOMPUTED_SETS.K[b][f].size ? PRECOMPUTED_SETS.K[b][f].data[k] : NULL, f);
                if (!PRECOMPUTED_SETS.K[b][f].data || !PRECOMPUTED_SETS.K[b][f].data[k])
                {
                    printf("here for %d and %d on the first if\n", b, f);
                    continue;
                }
                size_t *currentK = PRECOMPUTED_SETS.K[b][f].data[k];
                for (size_t h = 0; h < PRECOMPUTED_SETS.H[b][f].size; h++)
                {
                    size_t *currentH = PRECOMPUTED_SETS.H[b][f].data[h];
                    if (!ifAllElements(currentH, currentK))
                    { // Maybe could be optimized
                        printf("here for %d and %d on the second if\n", b, f);
                        continue;
                    }
                    // ---- Recursion: compute the value of `a` ----
                    double a = computeA(b, f, currentH, probabilities);
                    double valueBefore[TOTAL_CANDIDATES][TOTAL_GROUPS];

                    // ---- Get the value from u_b,g-1,g,c(k-h) ---- //
                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    {
                        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                        {
                            printf("THE F VALUE IS %d and the boolean returns %d\n", f,
                                   checkNull(currentK, (size_t)TOTAL_CANDIDATES));
                            sleep(1);
                            if (f == 1 && checkNull(currentK, (size_t)TOTAL_CANDIDATES))
                            {
                                valueBefore[c][g] = 1.0;
                                setMemoValue(memo, b, 0, g, c, currentK, TOTAL_CANDIDATES, 1.0);
                                printf("IM HERE");
                                sleep(4);
                            }
                            else if (f == 1 && !checkNull(currentK, (size_t)TOTAL_CANDIDATES))
                            {
                                valueBefore[c][g] = 0;
                                setMemoValue(memo, b, 0, g, c, currentK, TOTAL_CANDIDATES, 1.0);
                                printf("IM HERE, SECOND LOOP");
                                sleep(4);
                            }
                            else
                            {
                                size_t *substractionVector = malloc((size_t)TOTAL_CANDIDATES * sizeof(size_t));
                                vectorDiff(currentK, currentH, substractionVector);
                                printf("The vector diff is:\n");
                                printSizeTVector(substractionVector, TOTAL_CANDIDATES);
                                valueBefore[c][g] =
                                    getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                                free(substractionVector);
                            }
                        }
                    }
                    // ---...--- //

                    // ---- Set the new value ---- //
                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    {
                        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                        {
                            double valueNow = getMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES);
                            printf("When searching for u on index %d, %d, %d, %d the value now is %.4f\n", b, f, g, c,
                                   valueNow);
                            if (valueNow == INVALID)
                            { // If there wasn't a value before
                                valueNow = 0.0;
                            }
                            printf("The current value is:\t%.8f\nThe before value is:\t%.8f\n", valueNow,
                                   valueBefore[c][g]);
                            if (f == g)
                            {
                                valueNow += valueBefore[c][g] * a * currentH[c] /
                                            (MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f));
                            }
                            else
                            {
                                valueNow += valueBefore[c][g] * a;
                            }
#pragma omp critical
                            {
                                setMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES, valueNow);
                            }
                        }
                    }
                    // ---...--- //
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
    // #pragma omp parallel for schedule(dynamic) collapse(2)
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
                // printf("\nThe result stored for ballot %d, group %d and candidate %d is:\t%.8f\n", b, g, c, num /
                // den);
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
