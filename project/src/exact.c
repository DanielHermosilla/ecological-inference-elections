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
    uint32_t b;
    uint16_t g;
    size_t **data;
    size_t size;
} Set;

Set **HSETS = NULL;       // Global pointer to store all Hsets
Set **KSETS = NULL;       // Global pointer to store all Ksets
size_t **CANDIDATEARRAYS; // 2D array indexed by [c][b]

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

void generateHSets()
{
    HSETS = malloc(TOTAL_BALLOTS * sizeof(Set *));

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        HSETS[b] = malloc(TOTAL_GROUPS * sizeof(Set));
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            // Initialize H set
            HSETS[b][g].b = b;
            HSETS[b][g].g = g;

            // Parameters for the function
            size_t total = (size_t)MATRIX_AT_PTR(W, b, g);
            // Generate configurations
            size_t count = 0;
            size_t **configurations = generateAllConfigurations(b, total, TOTAL_CANDIDATES, &count);

            // Store configurations and size
            HSETS[b][g].data = configurations;
            HSETS[b][g].size = count;
        }
    }
}

void generateKSets()
{
    KSETS = malloc(TOTAL_BALLOTS * sizeof(Set *));

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        KSETS[b] = malloc(TOTAL_GROUPS * sizeof(Set));
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            // Initialize K set
            KSETS[b][g].b = b;
            KSETS[b][g].g = g;

            size_t total = 0;
            for (uint16_t f = 0; f <= g; f++)
            {
                total += (size_t)MATRIX_AT_PTR(W, b, f);
            }

            // Generate configurations
            size_t count = 0;
            size_t **configurations = generateAllConfigurations(b, total, TOTAL_CANDIDATES, &count);

            // Store configurations and size
            KSETS[b][g].data = configurations;
            KSETS[b][g].size = count;
        }
    }
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
            size_t **configurations = generateAllConfigurations(b, total, TOTAL_CANDIDATES, &count, true);

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
                int ktotal = 0;
                for (uint16_t f = 0; f <= g; f++)
                {
                    ktotal += (size_t)(MATRIX_AT_PTR(W, b, f));
                }
                size_t Kcount = 0;
                size_t **Kconfigurations = generateAllConfigurations(b, ktotal, TOTAL_CANDIDATES, &Kcount, true);
                KSETS[b][g].data = Kconfigurations;
                KSETS[b][g].size = Kcount;
            }
            printf("\nFor the last array on ballot %d and group %d we've got:\nHSET:\n", b, g);
            printSizeTVector(HSETS[b][g].data[HSETS[b][g].size], TOTAL_CANDIDATES);
            printf("\nKSET\n");
            printSizeTVector(KSETS[b][g].data[KSETS[b][g].size], TOTAL_CANDIDATES);
            sleep(1);
        }
    }
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
    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    {
        result -= gsl_sf_lngamma(hElement[i] + 1); // Divide by each h_i!
    }

    return exp(result);
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

void recursion2(MemoizationTable *memo, const Matrix *probabilities)
{
#pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
        {
            for (size_t k = 0; k < KSETS[b][f].size; k++)
            {
                if (!KSETS[b][f].data || !KSETS[b][f].data[k])
                {
                    continue;
                }

                size_t *currentK = KSETS[b][f].data[k];
                for (size_t h = 0; h < HSETS[b][f].size; h++)
                {
                    size_t *currentH = HSETS[b][f].data[h];
                    if (!ifAllElements(currentH, currentK))
                    { // Maybe could be optimized
                        continue;
                    }

                    // ---- Compute the values that are independent from c and g ---- //
                    double a = computeA(b, f, currentH, probabilities);
                    double valueBefore;
                    size_t *substractionVector = malloc((size_t)TOTAL_CANDIDATES * sizeof(size_t));
                    vectorDiff(currentK, currentH, substractionVector);
                    // --- ... --- //

                    // ---- Get the value from u_b,g-1,g,c(k-h) ---- //
                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    {
                        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                        {

                            // ---- Get the value from the last iteration ---- //
                            if (f == 0 && checkNull(substractionVector, TOTAL_CANDIDATES))
                            {
                                valueBefore = 1.0;
                            }
                            else if (f == 0)
                            {
                                valueBefore = 0.0;
                            }
                            else
                            {
                                valueBefore = getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                            }
                            // --- ... --- //

                            // ---- Get the current value ---- //
                            double valueNow = getMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES);
                            if (valueNow == INVALID) // ---- If there's not a value created, set it as zero. ----
                                valueNow = 0.0;
                            if (f == g)
                            {
                                if (MATRIX_AT_PTR(W, b, f) == 0)
                                    continue;
                                valueNow += valueBefore * a * currentH[c] /
                                            (MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f));
                            }
                            else
                            {
                                valueNow += valueBefore * a;
                            }
                            // --- ... --- //

                            // ---- Store the value ---- //
#pragma omp critical
                            {
                                setMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES, valueNow);
                            }
                            // --- ... --- //
                        }
                    }
                    free(substractionVector);

                    // ---...--- //

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
        CANDIDATEARRAYS = malloc(TOTAL_BALLOTS * sizeof(size_t *));
        // ---- Parallelization will consider the two inner loops as only one
        for (uint16_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            CANDIDATEARRAYS[b] = malloc(TOTAL_CANDIDATES * sizeof(size_t));
            for (uint32_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                CANDIDATEARRAYS[b][c] = (size_t)MATRIX_AT_PTR(X, c, b);
            }
        }
    }
    // ---- Initialize the `H` set, for recicling its values
    if (HSETS == NULL && KSETS == NULL)
    {
        printf("Starting the heavy precomputation\n");
        generateHSets();
        printf("\nHSETS is generated.\n");
        generateKSets();
        printf("\nKSETS is generated.\n");
        // printSets();
        printf("The heavy precomputation is finished\n");
    }

    printf("Starting the recursion\n");
    // ---- Initialize the hash table
    MemoizationTable *table = initMemo();
    // ---- Initialize the array to return
    double *array2 = (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double));
    recursion2(table, probabilities);
    printf("Finished the recursion\n");
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
#pragma omp critical
                {
                    den += getMemoValue(table, b, TOTAL_GROUPS - 1, g, c, CANDIDATEARRAYS[b], TOTAL_CANDIDATES) *
                           MATRIX_AT_PTR(probabilities, g, c);
                }
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate
              // ---- Compute the numerator; the recursion shouldn't be ran since the denominator initialized the
              // hashed values.
#pragma omp critical
                {
                    double num = getMemoValue(table, b, TOTAL_GROUPS - 1, g, c, CANDIDATEARRAYS[b], TOTAL_CANDIDATES) *
                                 MATRIX_AT_PTR(probabilities, g, c);
                    // ---- Store the resulting values as q_{bgc}
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = num / den;
                }
            }
        }
    }
    freeMemo(table);
    return array2;
}

void freeHSet()
{
    // Free HSETS
    if (HSETS != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            if (HSETS[b] != NULL)
            {
                for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                {
                    if (HSETS[b][g].data != NULL)
                    {
                        // Free each configuration in the Hset
                        for (size_t i = 0; i < HSETS[b][g].size; i++)
                        {
                            free(HSETS[b][g].data[i]);
                        }
                        free(HSETS[b][g].data); // Free the array of configurations
                    }
                }
                free(HSETS[b]); // Free the array of Hsets for this ballot
            }
        }
        free(HSETS); // Free the global array
        HSETS = NULL;
    }
}

void freeKSet()
{
    // Free KSETS
    if (KSETS != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            if (KSETS[b] != NULL)
            {
                for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                {
                    if (KSETS[b][g].data != NULL)
                    {
                        // Free each configuration in the Bset
                        for (size_t i = 0; i < KSETS[b][g].size; i++)
                        {
                            free(KSETS[b][g].data[i]);
                        }
                        free(KSETS[b][g].data); // Free the array of configurations
                    }
                }
                free(KSETS[b]); // Free the array of Bsets for this ballot
            }
        }
        free(KSETS); // Free the global array
        KSETS = NULL;
    }
}

__attribute__((destructor)) void cleanUp()
{
    // ---- Destroy the candidate array of size_t
    if (CANDIDATEARRAYS != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            free(CANDIDATEARRAYS[b]);
        }
        free(CANDIDATEARRAYS);
    }
    // ---- Destroy the array of precomputed H
    freeHSet();
    freeKSet();
}
