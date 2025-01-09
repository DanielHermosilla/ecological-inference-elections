#include "exact.h"
#include "globals.h"
#include "matrixUtils.h"
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
SizeTMatrix *CANDIDATEARRAYS = NULL;
gsl_combination ***PRECOMPUTED_H = NULL;

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

/*
void gsl_combination_free(gsl_combination *c)

    This function frees all the memory used by the combination c.
*/

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
    double result = 1;
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        result *= pow(MATRIX_AT_PTR(probabilities, f, c), hElement[c]);
    }
    return result;
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
    double result = factorial(MATRIX_AT_PTR(W, b, f)); // Calculate w_bf! GSL_SF_FACT HAVE A MAXIMUM VALUE OF 170.
    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    {
        result /= factorial(hElement[i]); // Divide by each h_i!
    }

    return result;
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
 * @brief Calculate a given `K` set.
 *
 * Given the index of the set, it returns the `K` set and its size. Note that it will be an array of arrays. Both arrays
 * are size_t.
 *
 * @param[in] b Index that represents the corresponding ballot.
 * @param[in] f Index that represents the first `f` groups (starting at 0).
 * @param[in, out] size The size of the new array.
 *
 * @return size_t **: A pointer towards the elements of the array of arrays.
 *
 * @note TODO: The while loop could be optimized.
 *
 * @warning
 * - The `f` index cannot be bigger than the total amount of groups.
 * - The `b` index cannot be bigger than the total ballots.
 *
 */

size_t **getK(const int b, const int f, size_t *size) // combination[i][j]; i the combination, j the element
{
    if (f >= (int)TOTAL_GROUPS || b >= (int)TOTAL_BALLOTS)
    {
        printf("An incorrect index was passed to the `K` set");
        exit(EXIT_FAILURE);
    }

    size_t **restriction = NULL;
    *size = 0;

    int totalVotes = 0; // The total votes accumulated for the first `f` groups
    for (int i = 0; i < f; i++)
    {
        totalVotes += MATRIX_AT_PTR(W, b, i);
    }
    gsl_combination *noRestriction = getH(totalVotes, TOTAL_CANDIDATES);

    do // It could be optimized by stepping on a big amount of elements.
    {
        const size_t *elementArray = gsl_combination_data(noRestriction);

        if (filterCombinations(elementArray, b))
        {
            restriction = realloc(restriction, (*size + 1) * sizeof(size_t *)); // The size of the whole array
            restriction[*size] = malloc(TOTAL_CANDIDATES * sizeof(size_t)); // The size of a single element (constant)

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                restriction[*size][c] = elementArray[c]; // Copy the element into the array
            }
        }

    } while (gsl_combination_next(noRestriction) == GSL_SUCCESS);

    gsl_combination_free(noRestriction);

    return restriction;
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
double recursion2(const int b, const int f, const int g, const int c, size_t *vector, MemoizationTable *memo,
                  const Matrix *probabilities)
{

    double result = 0;
    for (int16_t f = 0; f < TOTAL_GROUPS; f++)
    {

        gsl_combination *Hlocal = gsl_combination_alloc(PRECOMPUTED_H[b][f]->n, PRECOMPUTED_H[b][f]->k);
        gsl_combination_memcpy(Hlocal, PRECOMPUTED_H[b][f]);
        do // If it's too heavy, parallelize
        {
            // ---- Recursion: get a pointer of an element from H ----
            size_t *hElement = gsl_combination_data(Hlocal);
            // ---- Recursion: handle case when substraction could be negative (*) ----
            if (!ifAllElements(hElement, vector))
                continue;

            double a = computeA(b, f, hElement, probabilities);

            if (f == 1 && checkNull(hElement, (size_t)TOTAL_CANDIDATES))
            {
                double valueBefore = 1.0;
            }
            else if (f == 1 && !checkNull(hElement, (size_t)TOTAL_CANDIDATES))
            {
                double valueBefore = 0.0;
            }
            else
            {
                size_t *substractionVector = malloc((size_t)TOTAL_CANDIDATES * sizeof(size_t));
                vectorDiff(vector, hElement, substractionVector);
                double valueBefore = getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                free(substractionVector);
            }
            for (int16_t g = 0; g < TOTAL_GROUPS; g++)
            {
                double result = valueBefore +
            }

        } while (gsl_combination_next(Hlocal) == GSL_SUCCESS); // Recursion: condition to loop over each H element
    }
    return result;
    // ---- End recursion ----
}

double recursion(const int b, const int f, const int g, const int c, size_t *vector, MemoizationTable *memo,
                 const Matrix *probabilities)
{
    // If it's already defined, return the defined value instead; at each iteration this condition will be met more
    // times.
    printSizeTVector(vector, TOTAL_CANDIDATES);
    printf("\nUnder the indices of (%d, %d, %d, %d)", b, f, g, c);
    double currentVal = getMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES);
    if (currentVal != INVALID)
    {
        return currentVal;
    }

    // ---- Base case ----
    if (f == 0)
    {
        printf("\nIt's the base case");
        for (uint16_t k = 0; k < TOTAL_CANDIDATES; k++)
        { // ---- For each candidate
            if (vector[k] != 0)
            { // ---- If the vector contains an non zero element; return 0.0
#pragma omp critical
                { // ---- Race condition isn't a problem but it might cause memory corruption for the hash table
                    setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, 0.0);
                }
                return 0.0;
            }
        }
// ---- Up to this point, the vector should be full of zeros.
#pragma omp critical
        {
            setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, 1.0);
        }
        return 1.0;
    }
    // ---- End base case ----

    // ---- Recursion ----
    // ---- Recursion: definition of H ----
    // ---- There will be a copy for handling multithreaded loops, hence, the parallelization would be made starting
    // this point.
    double result = 0;
    gsl_combination *Hlocal = gsl_combination_alloc(PRECOMPUTED_H[b][f]->n, PRECOMPUTED_H[b][f]->k);
    gsl_combination_memcpy(Hlocal, PRECOMPUTED_H[b][f]);
    // ---- Recursion: definition of the resulting variable ----
    // ---- Recursion: loop over H, it must be a summatory ----
    // ---- The while loop can be parralelized if each thread get its own local copy of Hlocal
    do // If it's too heavy, parallelize
    {
        // ---- Recursion: get a pointer of an element from H ----
        size_t *hElement = gsl_combination_data(Hlocal);
        if (!hElement)
        {
            fprintf(stderr, "Error: hElement is NULL in b=%d, f=%d\n", b, f);
            exit(EXIT_FAILURE);
        }
        // ---- Recursion: handle case when substraction could be negative (*) ----
        if (!ifAllElements(hElement, vector))
            continue;
        // ---- Recursion: compute the value of `a` ----
        double a = computeA(b, f, hElement, probabilities);
        // ---- Recursion: define the vector k-h ----
        size_t *substractionVector = malloc((size_t)TOTAL_CANDIDATES * sizeof(size_t));
        vectorDiff(vector, hElement, substractionVector);
        printf("\nThe substracting vector is:\n");
        printSizeTVector(substractionVector, TOTAL_CANDIDATES);

        // ---- Recursion: handle case of the additional multiplication ----
        if (f == g)
        {
            double division = MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f);
            result += (recursion(b, f - 1, g, c, substractionVector, memo, probabilities) * a) / division;
        }
        else
        {
            result += recursion(b, f - 1, g, c, substractionVector, memo, probabilities) * a;
        }
        // ---- Recursion: Free the substracting vector
        free(substractionVector);
        // ---- Note: the element is defined as a summatory from all elements from H_bf that accomplish (*) ----

    } while (gsl_combination_next(Hlocal) == GSL_SUCCESS); // Recursion: condition to loop over each H element

    // ---- Recursion: free the copy of `H`
    gsl_combination_free(Hlocal);
// ---- Recursion: final step, set the hash table to the final value and return it
#pragma omp critical
    {
        setMemoValue(memo, b, f, g, c, vector, TOTAL_CANDIDATES, result);
    }
    // ---- Note: the vector that was used as a key will be free'd with the table, since it's a pointer.
    return result;
    // ---- End recursion ----
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
    if (PRECOMPUTED_H == NULL)
    {
        printf("Starting the heavy precomputation\n");
        // ---- Define memory for the b index
        PRECOMPUTED_H = malloc(TOTAL_BALLOTS * sizeof(gsl_combination **));
        // ---- Parallelization would be dynamic according the different workload
        // #pragma omp parallel for collapse(2) schedule(dynamic)
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // ---- For each ballot box
            // ---- Define memory for the f index
            PRECOMPUTED_H[b] = malloc(TOTAL_GROUPS * sizeof(gsl_combination *));
            for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
            { // ---- For each group
                // ---- Handle error case where n < k, on that case, just give a "1" vector. Maybe it could change
                if (MATRIX_AT_PTR(W, b, f) < TOTAL_CANDIDATES)
                {
                    PRECOMPUTED_H[b][f] = gsl_combination_calloc((int)TOTAL_CANDIDATES, (int)TOTAL_CANDIDATES);
                }
                else
                {
                    PRECOMPUTED_H[b][f] = gsl_combination_calloc(
                        MATRIX_AT_PTR(W, b, f),
                        (int)TOTAL_CANDIDATES); // Maybe, instead of casting as int use the candidate arrays
                }
            }
        }
    }

    printf("Starting the recursion\n");
    // ---- Initialize the hash table
    MemoizationTable *table = initMemo();
    // ---- Initialize the array to return
    double *array2 = (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double));
    // ---- Initialize the loop for computing each element of q_{bgc}
    // ---- The parallelization would be made on the outer loops. It might cause a problem if the recursion tries to
    // write at the same location, to be checked. #pragma omp parallel for schedule(dynamic) collapse(2)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot box
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // ---- Initialize the denominator variable to avoid multiple computations
            double den = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                // ---- Add the values of the denominator
                den +=
                    recursion(b, TOTAL_GROUPS - 1, g, c, &MATRIX_AT_PTR(CANDIDATEARRAYS, c, b), table, probabilities) *
                    MATRIX_AT_PTR(probabilities, g, c);
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
                // ---- Compute the numerator; the recursion shouldn't be ran since the denominator initialized the
                // hashed values.
                double num =
                    recursion(b, TOTAL_GROUPS - 1, g, c, &MATRIX_AT_PTR(CANDIDATEARRAYS, c, b), table, probabilities) *
                    MATRIX_AT_PTR(probabilities, g, c);
                // ---- Store the resulting values as q_{bgc}
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = num / den;
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
    if (PRECOMPUTED_H != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
            {
                gsl_combination_free(PRECOMPUTED_H[b][f]);
            }
            free(PRECOMPUTED_H[b]);
        }
        free(PRECOMPUTED_H);
        PRECOMPUTED_H = NULL;
    }
}