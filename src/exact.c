/*
Copyright (c) 2025 fastei team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "exact.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h> // for R_CheckUserInterrupt()
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#ifndef Realloc
#define Realloc(p, n, t) ((t *)R_chk_realloc((void *)(p), (size_t)((n) * sizeof(t))))
#endif

size_t **CANDIDATEARRAYS;

// Helper to free the [B][C] candidate arrays
static inline void freeCandidateArrays(uint32_t B)
{
    if (CANDIDATEARRAYS != NULL)
    {
        for (uint32_t b = 0; b < B; ++b)
        {
            if (CANDIDATEARRAYS[b] != NULL)
            {
                Free(CANDIDATEARRAYS[b]);
                CANDIDATEARRAYS[b] = NULL;
            }
        }
        Free(CANDIDATEARRAYS);
        CANDIDATEARRAYS = NULL;
    }
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
 * @brief Checks if the vector is null.
 *
 * @param[in] *vector A pointer towards the vector to check the condition.
 * @param[in] size The size of the vector.
 *
 * @return Boolean value that tells if the vector is null (true) or not (false).
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

/**
 * @brief Calculate the possible configurations via a recursion.
 *
 * Given the initial parameters, it will calculate an array with all of the possible voting outcomes, having the
 * constraint on the votes that a candidate had scored. The constraint will put an upper limit towards the dimension of
 * the given candidate.
 *
 * @param[in] b. The index of the ballot box that is going to be calculated.
 * @param[in] *votes. The array that will be storing every vote.
 * @param[in] position. The current position that indicates the index of the candidate that's being calculated.
 * @param[in] remainingVotes. The remaining votes to distribute over each dimension.
 * @param[in] numCandidates. The amount of candidates. It will also determine the dimensions of the array.
 * @param[in, out] ***results. The main array for storing each possible combination.
 * @param[in, out] *count. A counter that will have the amount of combinations
 *
 * @return void.
 *
 */

void generateConfigurations(EMContext *ctx, int b, size_t *votes, int position, int remainingVotes, int numCandidates,
                            size_t ***results, size_t *count)
{
    Matrix *X = &ctx->X; // Get the X matrix
    // ---- Base case: we're on the last candidate ---- //
    if (position == numCandidates - 1)
    {
        // ---- Assign remaining votes to the last candidate ----
        votes[position] = remainingVotes;

        // ---- If the last candidate actually had less votes, ditch that combination ----
        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // ---- Exit the recursion and don't save anything
            return;
        }

        // ---- Store the result ---- //
        // ---- Up to this point, the combination is valid, hence, the results will be stored.
        (*results) = Realloc(*results, (*count + 1), size_t *);
        (*results)[*count] = Calloc(numCandidates, size_t);
        memcpy((*results)[*count], votes, numCandidates * sizeof(size_t));
        (*count)++;
        return;
        // ---...--- //
    }
    // ---...--- //

    // ---- Loop over all the remaining votes ---- //
    for (int i = 0; i <= remainingVotes; i++)
    { // ---- For each remaining vote
        // ---- Assing that amount of votes to the candidate in the given position ----
        votes[position] = i;

        // ---- If the candidate actually had less votes, ditch that combination ----
        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // ---- Exit the recursion and dont save anything
            return;
        }
        // ---- Call the recursion ----
        generateConfigurations(ctx, b, votes, position + 1, remainingVotes - i, numCandidates, results, count);
    }
    // ---...--- //
}

/**
 * @brief Main function for generating all combinations possible.
 *
 * Given that the combinations are created towards a recursion, the main function is to add as a wrapper towards the
 * recursion function. It will return a pointer with all of the possible configurations. The combinations have an upper
 * constraint given the results that a candidate have gotten.
 *
 * @param[in] b. The index of the ballot box that is going to be calculated.
 * @param[in] totalVotes. The total amount of votes to handle. For example, if a group had `10` votes, the sum of each
 * element will be 10.
 * @param[in] numCandidates. The amount of candidates. It will also determine the dimensions of the array.
 * @param[in, out] *count. Pointer that will store the total amount of combinations. Useful for iterating over a set.
 *
 * @return size_t **: A pointer that will store arrays of arrays, having all of the possible combinations.
 *
 */
size_t **generateAllConfigurations(EMContext *ctx, int b, int totalVotes, int numCandidates, size_t *count)
{
    // ---- Initialize parameters ---- //
    size_t **results = NULL;
    *count = 0;
    size_t *votes = Calloc(numCandidates, size_t);
    // --- ... --- //

    // ---- Call the recursion ---- //
    generateConfigurations(ctx, b, votes, 0, totalVotes, numCandidates, &results, count);
    // --- ... --- //
    Free(votes);
    return results;
}

/**
 * @brief Precalculates the `H` set for every index.
 *
 * Given that the set wasn't calculated before, it calculates the H set defined as every possible combination for a
 * given `g` group
 *
 * @return void: Results written at the global variable HSETS.
 *
 */
void generateHSets(EMContext *ctx)
{
    Matrix *W = &ctx->W; // Get the W matrix
                         // ---- Allocate memory for the `b` index ----
    ctx->hset = Calloc(TOTAL_BALLOTS * TOTAL_GROUPS, Set);
    Set *HSETS = ctx->hset;

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box

        // ---- Allocate memory for the `g` index ----

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group given a ballot box

            // ---- Initialize H set and set initial parameters ---- //
            HSET(ctx, b, g)->b = b;
            HSET(ctx, b, g)->g = g;
            // ---- Parameters for the function ----
            size_t total = (size_t)MATRIX_AT_PTR(W, b, g);
            // --- ... --- //

            // ---- Compute the set combinations ---- //
            size_t count = 0;
            size_t **configurations = generateAllConfigurations(ctx, b, total, TOTAL_CANDIDATES, &count);

            // ---- Store configurations and size ----
            HSET(ctx, b, g)->data = configurations;
            HSET(ctx, b, g)->size = count;
            // --- ... --- //
        }
    }
}

/**
 * @brief Precalculates the `K` set for every index.
 *
 * Given that the set wasn't calculated before, it calculates the K set defined as every possible combination as a
 * cummulative set given the first `f` groups.
 *
 * @return void: Results written at the global variable KSETS.
 *
 */

void generateKSets(EMContext *ctx)
{
    Matrix *W = &ctx->W; // Get the W matrix
    // ---- Allocate memory for the `b` index ----
    ctx->kset = Calloc(TOTAL_BALLOTS * TOTAL_GROUPS, Set);
    Set *KSETS = ctx->kset;

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box

        // ---- Allocate memory for the `f` index ----
        // KSETS[b] = Calloc(TOTAL_GROUPS, Set);

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group given a ballot box

            // ---- Initialize K set and set initial parameters ---- //
            KSET(ctx, b, g)->b = b;
            KSET(ctx, b, g)->g = g;
            // ---- Parameters for the function ----
            size_t total = 0;
            for (uint16_t f = 0; f <= g; f++)
            { // --- For each cummulative group
                // ---- Sum the amount of votes from the first `f` groups ----
                total += (size_t)MATRIX_AT_PTR(W, b, f);
            }
            // --- ... --- //

            // ---- Compute the set combinations ---- //
            size_t count = 0;
            size_t **configurations = generateAllConfigurations(ctx, b, total, TOTAL_CANDIDATES, &count);

            // ---- Store configurations and size ----
            KSET(ctx, b, g)->data = configurations;
            KSET(ctx, b, g)->size = count;
            // --- ... --- //
        }
    }
}

/**
 * @brief Calculate the chained product between probabilities as defined in `a`:
 *
 * Given an `H` element and the `f` index it computes the chained product. It will use logarithmics for reducing
 * complexity.
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
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate

        // ---- Get the matrix value at `f` and `c` ----
        double prob = MATRIX_AT_PTR(probabilities, f, c);

        // ---- If the multiplication gets to be zero, it would be undefined in the logarithm (and zero for all of the
        // chain). Do an early stop if that happens ----
        if (prob == 0.0 && hElement[c] > 0)
        {
            return 0.0; // ---- Early stop ----
        }

        // ---- Ensures the probability is greater than zero for not getting a undefined logarithm. Anyways, that
        // shouldn't happen. ----
        if (prob > 0)
        {
            // ---- Add the result to the logarithm. Remember that by logarithm properties, the multiplication of the
            // arguments goes as a summatory ----
            log_result += hElement[c] * log(prob);
        }
    }
    // --- ... --- //

    // ---- Exponetiate the final result ----
    return exp(log_result);
}

/**
 * @brief Calculate the multinomial coefficient, subject to the `H` set:
 *
 * Given an `H` element, it computes its multinomial coefficient subject to the total amount of votes a group received.
 * It yields the following calculation:
 * $$\binom{w_{bf}}{h_1,\dots,h_C}=\frac{w_{bf}!}{h_1!\cdots h_{C}!}$$
 * It will use logarithms for the factorials, and later be converted to their exponential equivalent.
 *
 * @param [in] b The index of the ballot.
 * @param [in] f The index of the group
 * @param[in] *hElement A pointer towards the element of `H`. This element is an array with the possible combinations
 *
 * @return double: The result of the multinomial coefficient
 *
 */
double multinomialCoeff(EMContext *ctx, const int b, const int f, const size_t *hElement)
{
    Matrix *W = &ctx->W; // Get the W matrix
    // --- Compute ln(w_bf!). When adding by one, it considers the last element too ---
    double result = lgamma1p((int)MATRIX_AT_PTR(W, b, f));

    for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
    { // ---- For each candidate
        // ---- Divide by each h_i! ----
        result -= lgamma1p(hElement[i]);
    }
    // ---- Return the original result by exponentiating ----
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
double computeA(EMContext *ctx, const int b, const int f, const size_t *hElement)
{
    Matrix *probabilities = &ctx->probabilities; // Get the probabilities matrix
    return multinomialCoeff(ctx, b, f, hElement) * prod(hElement, probabilities, f);
}

/**
 * @brief Calculates the main loop of the pseudocode and store the values on the hash table.
 *
 * It will loop over every value and add the results in the memoization table. It will be parallelized just on the outer
 * loop, benchmarking has shown that parallelizing over `k` and `h` doesn't win anything. Anyways, it could be reverted.
 *
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 * @param[in, out] *memo A pointer towards the hash table.
 *
 * @return void. Results to be written on the hash table
 *
 */
void recursion(EMContext *ctx, MemoizationTable *memo)
{
    Matrix *probabilities = &ctx->probabilities; // Get the probabilities matrix
    Matrix *W = &ctx->W;                         // Get the W matrix
    Set *KSETS = ctx->kset;
    Set *HSETS = ctx->hset;

    // #pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot box
        // if (b % 5 == 0) // Checks condition every 5 iterations
        // R_CheckUserInterrupt();
        for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
        { // ---- For each group, given a ballot box
            // #pragma omp parallel for collapse(2)
            // It actually worsened the performance
            for (size_t k = 0; k < KSET(ctx, b, f)->size; k++)
            { // ---- For each element from the K_bf set
                // ---- If there's no existing combination, skip the loop ----
                if (!KSET(ctx, b, f)->data || !(KSET(ctx, b, f)->data[k]) || (KSET(ctx, b, f)->data[k]) == NULL)
                {
                    continue;
                }

                // ---- Define the current element from the K set ----
                // size_t *currentK = KSETS[b][f].data[k];
                size_t *currentK = KSET(ctx, b, f)->data[k];

                for (size_t h = 0; h < HSET(ctx, b, f)->size; h++)
                { // ---- For each element from the H_bf set
                    // if (HSET(ctx, b, f)->size > 5000 && h % 250 == 0)
                    // R_CheckUserInterrupt();
                    size_t *currentH = HSET(ctx, b, f)->data[h];
                    // ---- If the element from h isn't smaller than the one from k ----
                    // ---- Note that, when generating the H set, the restriction from candidate votes was also imposed,
                    // so it excluded "trivial" cases ----
                    if (!ifAllElements(currentH, currentK))
                    { // ---- Maybe could be optimized
                        continue;
                    }

                    // ---- Compute the values that are independent from c and g ---- //
                    // ---- The value `a` from the pseudocode. Check `computeA` for more information ----
                    double a = computeA(ctx, b, f, currentH);
                    // ---- Initialize the variable that will store the past iteration ----
                    double valueBefore;
                    // ---- Substract the Kth element with the Hth element (k-h) ----
                    size_t *substractionVector = Calloc((size_t)TOTAL_CANDIDATES, size_t); // ---- Allocates memory
                    vectorDiff(currentK, currentH, substractionVector);
                    // --- ... --- //

                    // ---- Retrieve the initial and current values  ---- //
                    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                    { // ----  For each candidate
                        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                        { // ---- For each group given a candidate

                            // ---- Get the value from the last iteration ---- //

                            // ---- Border case, if f = 0 we assume that it's initialized as 1 if the substracting
                            // vector is one ----
                            if (f == 0 && checkNull(substractionVector, TOTAL_CANDIDATES))
                            {
                                valueBefore = 1.0;
                            }
                            // ---- Another border case, if the substracting vector isn't 0 and f = 0, then the
                            // initialized value is zero ----
                            else if (f == 0)
                            {
                                valueBefore = 0.0;
                            }
                            // ---- No more border cases, at this point, every past value should had been defined ----
                            else
                            {
                                // #ifdef _OPENMP
                                // #pragma omp critical
                                // #endif
                                {
                                    valueBefore =
                                        getMemoValue(memo, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                                }
                            }
                            // --- ... --- //
                            double valueNow;
                            // ---- Get the current value ---- //
                            // #ifdef _OPENMP
                            // #pragma omp critical
                            // #endif
                            {
                                valueNow = getMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES);
                            }
                            // ---- If there's not a value created, set it as zero. ----
                            if (valueNow == INVALID)
                                valueNow = 0.0;
                            // --- ... --- //

                            // ---- Calculate the new value ---- //

                            if (f == g)
                            {
                                // ---- Precompute the denominator to avoid divitions by zero. ----
                                double den = MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f);
                                if (den == 0)
                                    continue;
                                // ---- Add the new value ----
                                valueNow += valueBefore * a * currentH[c] / den;
                            }
                            else
                            {
                                valueNow += valueBefore * a;
                            }
                            // --- ... --- //

                            // ---- Store the value ---- //
                            // ---- We set a critical point in case there's an error when inserting two values at the
                            // same time. Remember that this loop is parallelized ----
#ifdef _OPENMP
#pragma omp critical
#endif
                            {
                                setMemoValue(memo, b, f, g, c, currentK, TOTAL_CANDIDATES, valueNow);
                            }
                            // --- ... --- //
                        } // ---- End g loop
                    } // ---- End c loop
                    Free(substractionVector);
                    // ---...--- //
                } // ---- End H set loop
            } // ---- End K set loop
        } // ---- End `f` loop
    } // ---- End ballot boxes loop
}

/*
 * Gets the log-likelihood of the exact method, assuming the computation of the exact method is ongoing.
 *
 * TODO: Make an alternative version where we do not assume the exact method is being computed.
 */
static inline double exactLL_one_b(MemoizationTable *memo_b, uint32_t b)
{
    double px =
        getMemoValue(memo_b, b, TOTAL_GROUPS - 1, 0, TOTAL_CANDIDATES - 1, CANDIDATEARRAYS[b], TOTAL_CANDIDATES);
    return (fabs(px) > 0 ? log(fabs(px)) : 0.0);
}

/**
 * @brief Calculates the main loop of the pseudocode for a single ballot box.
 *
 * It will loop over every value and add the results in the memoization table, but only for a single `b`.
 * This version is meant to be executed in parallel over `b`, each one with its own `MemoizationTable *memo_b`,
 * avoiding any race conditions between threads.
 *
 * @param[in]  *ctx   A pointer to the EMContext with all matrices and parameters.
 * @param[in]   b     Index of the current ballot box to process.
 * @param[in,out] *memo_b A pointer toward the hash table for this specific ballot box.
 *
 * @return void. Results to be written on the hash table.
 *
 */
void recursion_one_b(EMContext *ctx, uint32_t b, MemoizationTable *memo_b)
{
    Matrix *probabilities = &ctx->probabilities; // Get the probabilities matrix
    Matrix *W = &ctx->W;                         // Get the W matrix
    Set *KSETS = ctx->kset;
    Set *HSETS = ctx->hset;

    for (uint16_t f = 0; f < TOTAL_GROUPS; f++)
    { // ---- For each group, given a ballot box
        for (size_t k = 0; k < KSET(ctx, b, f)->size; k++)
        { // ---- For each element from the K_bf set
            // ---- If there's no existing combination, skip the loop ----
            if (!KSET(ctx, b, f)->data || !(KSET(ctx, b, f)->data[k]) || (KSET(ctx, b, f)->data[k]) == NULL)
            {
                continue;
            }

            // ---- Define the current element from the K set ----
            size_t *currentK = KSET(ctx, b, f)->data[k];

            for (size_t h = 0; h < HSET(ctx, b, f)->size; h++)
            { // ---- For each element from the H_bf set
                size_t *currentH = HSET(ctx, b, f)->data[h];
                // ---- If the element from h isn't smaller than the one from k ----
                // ---- Note that, when generating the H set, the restriction from candidate votes was also imposed,
                // so it excluded "trivial" cases ----
                if (!ifAllElements(currentH, currentK))
                { // ---- Maybe could be optimized
                    continue;
                }

                // ---- Compute the values that are independent from c and g ---- //
                // ---- The value `a` from the pseudocode. Check `computeA` for more information ----
                double a = computeA(ctx, b, f, currentH);
                // ---- Initialize the variable that will store the past iteration ----
                double valueBefore;
                // ---- Substract the Kth element with the Hth element (k-h) ----
                size_t *substractionVector = Calloc((size_t)TOTAL_CANDIDATES, size_t); // ---- Allocates memory
                vectorDiff(currentK, currentH, substractionVector);
                // --- ... --- //

                // ---- Retrieve the initial and current values  ---- //
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // ----  For each candidate
                    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                    { // ---- For each group given a candidate

                        // ---- Get the value from the last iteration ---- //

                        // ---- Border case, if f = 0 we assume that it's initialized as 1 if the substracting
                        // vector is one ----
                        if (f == 0 && checkNull(substractionVector, TOTAL_CANDIDATES))
                        {
                            valueBefore = 1.0;
                        }
                        // ---- Another border case, if the substracting vector isn't 0 and f = 0, then the
                        // initialized value is zero ----
                        else if (f == 0)
                        {
                            valueBefore = 0.0;
                        }
                        // ---- No more border cases, at this point, every past value should had been defined ----
                        else
                        {
                            valueBefore = getMemoValue(memo_b, b, f - 1, g, c, substractionVector, TOTAL_CANDIDATES);
                        }
                        // --- ... --- //
                        double valueNow;
                        // ---- Get the current value ---- //
                        valueNow = getMemoValue(memo_b, b, f, g, c, currentK, TOTAL_CANDIDATES);
                        // ---- If there's not a value created, set it as zero. ----
                        if (valueNow == INVALID)
                            valueNow = 0.0;
                        // --- ... --- //

                        // ---- Calculate the new value ---- //
                        if (f == g)
                        {
                            // ---- Precompute the denominator to avoid divisions by zero. ----
                            double den = MATRIX_AT_PTR(probabilities, f, c) * MATRIX_AT_PTR(W, b, f);
                            if (den == 0)
                                continue;
                            // ---- Add the new value ----
                            valueNow += valueBefore * a * currentH[c] / den;
                        }
                        else
                        {
                            valueNow += valueBefore * a;
                        }
                        // --- ... --- //

                        // ---- Store the value ---- //
                        // ---- We do NOT set a critical section here since each b runs on its own memo table ----
                        setMemoValue(memo_b, b, f, g, c, currentK, TOTAL_CANDIDATES, valueNow);
                        // --- ... --- //
                    } // ---- End g loop
                } // ---- End c loop
                Free(substractionVector);
                // ---...--- //
            } // ---- End H set loop
        } // ---- End K set loop
    } // ---- End `f` loop
}

void computeQExact(EMContext *ctx, QMethodInput params, double *ll)
{
    Matrix *probabilities = &ctx->probabilities;
    Matrix *X = &ctx->X;
    double *q = ctx->q;

    // ---------- Build CANDIDATEARRAYS once ----------
    if (CANDIDATEARRAYS == NULL)
    {
        CANDIDATEARRAYS = Calloc(TOTAL_BALLOTS, size_t *);
        for (uint16_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            CANDIDATEARRAYS[b] = Calloc(TOTAL_CANDIDATES, size_t);
            for (uint32_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                CANDIDATEARRAYS[b][c] = (size_t)MATRIX_AT_PTR(X, c, b);
            }
        }
    }

    double ll_sum = 0.0;

    // ---------- Parallelize over b with owner-computes ----------
    // Each thread creates, uses, and frees its own memo table.
    // No critical sections needed; reduction for LL.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+ : ll_sum)
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        MemoizationTable *memo_b = initMemo();

        // Fill memo for this b
        recursion_one_b(ctx, b, memo_b);

        // Compute q_{b,g,c} using this b's memo
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            double den = 0.0;
            double num[TOTAL_CANDIDATES];

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double m = getMemoValue(memo_b, b, TOTAL_GROUPS - 1, g, c, CANDIDATEARRAYS[b], TOTAL_CANDIDATES);

                num[c] = m * MATRIX_AT_PTR(probabilities, g, c);
                den += num[c];
            }

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double result = (den != 0.0) ? (num[c] / den) : 0.0;
                if (isnan(result) || isinf(result))
                    result = 0.0;

                Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = result;
            }
        }

        // Optional: add to LL
        if (params.computeLL)
        {
            ll_sum += exactLL_one_b(memo_b, b);
        }

        freeMemo(memo_b);
    } // end parallel for b

    if (params.computeLL)
        *ll = ll_sum;

    freeCandidateArrays(TOTAL_BALLOTS);
}
