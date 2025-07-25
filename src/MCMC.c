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

#include "MCMC.h"
#include "globals.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h> // for R_CheckUserInterrupt()
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

// ctx->omegaset **OMEGASET = NULL; // Global pointer to store all H sets
// double **multinomialVals = NULL;
// double *logGammaArr = NULL;
// double **Qconstant = NULL;

int lessThanColRow(IntMatrix mat, int b, int g, int c, int candidateVotes, int groupVotes)
{
    int groupSum = 0;
    int canSum = 0;
    for (uint16_t i = 0; i < TOTAL_GROUPS; i++)
    {
        canSum += MATRIX_AT(mat, i, c);
    }
    for (uint16_t j = 0; j < TOTAL_CANDIDATES; j++)
    {
        groupSum += MATRIX_AT(mat, g, j);
    }
    int slackC = candidateVotes - canSum;
    int slackG = groupVotes - groupSum;

    return MIN(slackC, slackG);
}

IntMatrix startingPoint3(EMContext *ctx, int b)
{
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    // ---- Retrieve the initial variables ---- //
    IntMatrix toReturn = createMatrixInt(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *groupVotes = getRow(W, b);
    double *candidateVotes = getColumn(X, b);

    // ---- Calculate the expected value ---- //
    double totalC = 0;
    double totalG = 0;
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double mult = groupVotes[g] * candidateVotes[c];

            if (g == 0)
                totalC += candidateVotes[c];

            // In case of mismatch, we divide for the maximum

            MATRIX_AT(toReturn, g, c) = mult;
        }
        totalG += groupVotes[g];
    }
    // ---...--- //

    // ---- Division for mismatchs ---- //
    double divide = MAX(ctx->ballots_votes[b], totalC);
    divide = MAX(divide, totalG);

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double newValue = MATRIX_AT(toReturn, g, c) / divide;
            double floored = floor(newValue);
            MATRIX_AT(toReturn, g, c) = floored;
        }
    }

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            int groupRestriction = groupVotes[g];
            int candidateRestriction = candidateVotes[c];

            int m = lessThanColRow(toReturn, b, g, c, candidateRestriction, groupRestriction);
            if (m > 0)
            {
                MATRIX_AT(toReturn, g, c) += m;
            }
        }
    }

    // ---...--- //
    Free(groupVotes);
    Free(candidateVotes);
    return toReturn;
}

// 0 y C*(C-1), G(G-1)-1, o sea, sin los divididos por 2, la razón, es porque creo que es más sencillo encodear y
// decodear de la segunda manera
void allocateRandoms(int M, int S, uint8_t **c1, uint8_t **c2, uint8_t **g1, uint8_t **g2)
{
    uint32_t size = M * S;
    // Allocate memory correctly
    *c1 = (uint8_t *)Calloc(size, uint8_t);
    *c2 = (uint8_t *)Calloc(size, uint8_t);
    *g1 = (uint8_t *)Calloc(size, uint8_t);
    *g2 = (uint8_t *)Calloc(size, uint8_t);

    GetRNGstate(); // Ensure R's RNG is properly initialized
                   // Fill arrays with random indices
    int allow_repeat = (TOTAL_CANDIDATES <= 1 || TOTAL_GROUPS <= 1);

    for (int i = 0; i < size; i++)
    {
        if (i % 400 == 0) // Checks condition every 400 iterations
            R_CheckUserInterrupt();
        (*c1)[i] = (uint8_t)(unif_rand() * TOTAL_CANDIDATES);
        (*g1)[i] = (uint8_t)(unif_rand() * TOTAL_GROUPS);
        do
        {
            (*c2)[i] = (uint8_t)(unif_rand() * TOTAL_CANDIDATES);
            (*g2)[i] = (uint8_t)(unif_rand() * TOTAL_GROUPS);
        } while (!allow_repeat && ((*c2)[i] == (*c1)[i] || (*g2)[i] == (*g1)[i]));
        // 	} while(!allow_repeat);
    }
    PutRNGstate(); // Finalize RNG state to prevent repeatability
}

/*
 * @brief Precomputes the sets used for the simulation.
 *
 * Precomputes the sets that are independent from each EM iteration. It is made with parallelism (NOT SUPPORTED)
 * towards the ballot boxes and with a static assignment for ensuring reproducibility.
 *
 * @param[in] M. The step size between consecutive samples. Note that the direction is assigned randomly.
 * @param[in] S. The amount of samples for each ballot box.
 *
 * @return void. Written on the global variable.
 */
void generateOmegaSet(EMContext *ctx, int M, int S)
{
    // ---- Allocate memory for the `b` index ----
    ctx->omegaset = Calloc(TOTAL_BALLOTS, OmegaSet *);
    uint8_t *c1 = NULL;
    uint8_t *c2 = NULL;
    uint8_t *g1 = NULL;
    uint8_t *g2 = NULL;

    uint32_t arraySize = M * S;

    allocateRandoms(M, S, &c1, &c2, &g1, &g2);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

// ---- Perform the main iterations ---- //
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // if (b % 5 == 0)             // Checks condition every 5 iterations
        // R_CheckUserInterrupt(); // This might be fatal, since it doesn't free global memory
        // ---- Define a seed, that will be unique per thread ----
        //    unsigned int seed = rand_r(&seedNum) + omp_get_thread_number();
        // ---- Allocate memory for the ctx->omegaset ---- //
        ctx->omegaset[b] = Calloc(1, OmegaSet);
        ctx->omegaset[b]->b = b;
        ctx->omegaset[b]->size = S;
        ctx->omegaset[b]->data = Calloc(S, IntMatrix);
        // ctx->omegaset[b]->data = Calloc(S, IntMatrix *);
        // ---...--- //
        // ---- The `base` element used as a starting point ----
        IntMatrix startingZ = startingPoint3(ctx, b);
        int ballotShift = floor(((double)b / TOTAL_BALLOTS) * (M * S));

        // Impose the first step
        // IntMatrix *append = Calloc(1, IntMatrix);
        // *append = copMatrixI(&startingZ);
        // ctx->omegaset[b]->data[0] = append;A
        ctx->omegaset[b]->data[0] = copMatrixI(&startingZ);
        freeMatrixInt(&startingZ);

        for (int s = 1; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            // IntMatrix pastMatrix = ctx->omegaset[b]->data[s - 1];
            // IntMatrix steppingZ = copMatrixI(&pastMatrix);
            IntMatrix steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            for (int m = 0; m < M; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (M * S);
                uint8_t randomCDraw = c1[shiftIndex];
                uint8_t randomCDraw2 = c2[shiftIndex];
                uint8_t randomGDraw = g1[shiftIndex];
                uint8_t randomGDraw2 = g2[shiftIndex];

                // decode(randomCDraw, TOTAL_CANDIDATES, &c1, &c2);
                //  decode(randomGDraw, TOTAL_GROUPS, &g1, &g2);

                // ---- Check non negativity condition ---- //
                double firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                double secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                    continue;
                // ---...--- //

                // ---- Asign changes on the new matrix ---- //
                MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= 1;
                MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= 1;
                MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += 1;
                MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += 1;
                //  ---...--- //
            } // --- End the step size loop
            // ---- Add the combination to the ctx->omegaset ---- //
            // IntMatrix *append = Calloc(1, IntMatrix);
            // *append = copMatrixI(&steppingZ);
            ctx->omegaset[b]->data[s] = steppingZ;
            // freeMatrixInt(&steppingZ);
            // ---...--- //
        } // --- End the sample loop
        // freeMatrix(&startingZ);
    } // --- End the ballot box loop
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void encode(EMContext *ctx)
{
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        OmegaSet *set = ctx->omegaset[b];
        int S = set->size;

        MatrixHash *hashTable = NULL;

        // First pass: identify and count unique matrices
        for (int s = 0; s < S; s++)
        {
            IntMatrix current = set->data[s];
            unsigned int key = computeMatrixKey(&current);

            MatrixHash *entry = NULL;
            HASH_FIND(hh, hashTable, &key, sizeof(unsigned int), entry);

            int found = 0;
            for (MatrixHash *iter = entry; iter != NULL; iter = iter->hh.next)
            {
                if (iter->key == key && matricesAreEqualI(iter->matrix, &current))
                {
                    iter->count++;
                    found = 1;
                    break;
                }
            }

            if (!found)
            {
                MatrixHash *newEntry = Calloc(1, MatrixHash);
                newEntry->key = key;
                newEntry->matrix = &current;
                newEntry->count = 1;
                HASH_ADD(hh, hashTable, key, sizeof(unsigned int), newEntry);
            }
        }

        // Allocate counts and assign total count per sample
        if (set->counts == NULL)
        {
            set->counts = Calloc(S, int);
        }

        for (int s = 0; s < S; s++)
        {
            IntMatrix current = set->data[s];
            unsigned int key = computeMatrixKey(&current);

            MatrixHash *entry = NULL;
            HASH_FIND(hh, hashTable, &key, sizeof(unsigned int), entry);

            for (MatrixHash *iter = entry; iter != NULL; iter = iter->hh.next)
            {
                if (iter->key == key && matricesAreEqualI(iter->matrix, &current))
                {
                    set->counts[s] = iter->count; // All duplicates now get the same count
                    break;
                }
            }
            if (set->counts[s] == 0) // Hash collision case
                set->counts[s] = 1;
        }

        // Clean up
        MatrixHash *entry, *tmp;
        HASH_ITER(hh, hashTable, entry, tmp)
        {
            HASH_DEL(hashTable, entry);
            Free(entry);
        }
    }
}
/**
 * @brief Computes the pre-computable values of the expression that doesn't depend on EM iterations
 *
 * Given a ballot box index and a matrix represent an element from the Hit and Run ctx->omegaset, it computes the
 * following:
 *
 * $$\Prod_{g'\in G}\binom{w_{bg'}}{z_{bg'1}\cdots z_{bg'C}}$$
 *
 * @param[in] b The ballot box index
 * @param[in] *currentMatrix A pointer towards the current matricial element, of size GxC.
 *
 * double The result of the calculation.
 */
double preMultinomialCoeff(EMContext *ctx, const int b, IntMatrix *currentMatrix)
{
    Matrix *W = &ctx->W;          // Get the W matrix
    IntMatrix *intW = &ctx->intW; // Get the W matrix in integer form

    double result = 0;
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        // --- Compute ln(w_bf!). When adding by one, it considers the last element too ---
        // result += gsl_sf_lngamma((int)MATRIX_AT_PTR(W, b, g) + 1);
        result += lgamma1p(MATRIX_AT_PTR(intW, b, g)); // TODO: Could be saved

        for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate
            // ---- Divide by each h_i! ----
            // result -= gsl_sf_lngamma(MATRIX_AT_PTR(currentMatrix, g, i) + 1);
            result -= lgamma1p(MATRIX_AT_PTR(currentMatrix, g, i));
        }
    }
    // ---- Return the original result by exponentiating ----
    return result; // TODO: This can be avoided!
}

/**
 * @brief Calculates the last term of the multiplication ctx->omegaset
 *
 * Given a probability matrix, a ballot index and a ctx->omegaset index, it calculates:
 *
 * $$\Prod_{g\in G}\Prod_{c\in C}p_{gc}^{z_{bgc}}$$
 *
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] b The index of the ballot box
 * @param[in] setIndex The index of the ctx->omegaset
 *
 * @return double: The result of the product
 *
 */
double logarithmicProduct(EMContext *ctx, const int b, const int setIndex)
{
    Matrix *probabilities = &ctx->probabilities; // Get the probabilities matrix
    // ---- Define initial parameters ---- //
    double log_result = 0;
    IntMatrix currentMatrix = ctx->omegaset[b]->data[setIndex]; // Z
    // ---...--- //
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // Cambiar a log(1) if probabilidad 0
            // Producto punto
            double prob = MATRIX_AT_PTR(probabilities, g, c);
            log_result += (prob > 0.0) ? MATRIX_AT(currentMatrix, g, c) * log(prob) : 0.0;
        }
    }
    // --- ... --- //
    return log_result;
}

/**
 * @brief Precomputes the multinomial multiplication that is independent for each EM iteration.
 *
 * Calls the main function for computing all of the calculations related with the final result that are independent
 * from each EM call. Specifically, for each ballot box and its simulations, the following is calculated:
 *
 * $$\Prod_{g\in G}\binom{w_{bg}}{z_{bg1},\cdots, z_{bgC}}$$
 *
 * It transform the main product and the factorial to logarithmic scale for making efficient calculations.
 *
 * @return. Results written at the global variable
 */
void preComputeMultinomial(EMContext *ctx)
{
    // ---- Initialize space for storing all of the simulations ---- //
    ctx->multinomial = Calloc(TOTAL_BALLOTS, double *);
    // ---...--- //
    // ---- Compute the simulated combinations for each ctx->omegaset ---- //
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        // ---- Define the current ctx->omegaset and allocate memory for saving its size ----
        OmegaSet *currentSet = ctx->omegaset[b];
        ctx->multinomial[b] = Calloc(currentSet->size, double);
        for (size_t s = 0; s < currentSet->size; s++)
        { // --- For each simulation given a balot box
            ctx->multinomial[b][s] = preMultinomialCoeff(ctx, b, &currentSet->data[s]);
        }
    }
}

void precomputeQConstant(EMContext *ctx, int size)
{
    IntMatrix *intW = &ctx->intW; // Get the W matrix in integer form
    ctx->Qconstant = Calloc(TOTAL_BALLOTS, double *);

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        // ---- Define the current ctx->omegaset and allocate memory for saving its size ----
        OmegaSet *currentSet = ctx->omegaset[b];
        ctx->Qconstant[b] = Calloc(currentSet->size, double);
        for (size_t s = 0; s < currentSet->size; s++)
        { // --- For each simulation given a balot box
            IntMatrix sampleValue = currentSet->data[s];
            ctx->Qconstant[b][s] = 0;
            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            {
                int W_bg = MATRIX_AT_PTR(intW, b, g);
                ctx->Qconstant[b][s] += ctx->logGamma[W_bg];
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    ctx->Qconstant[b][s] -= ctx->logGamma[(int)MATRIX_AT(sampleValue, g, c)];
                }
            }
        }
    }
}

/*
 * Computes the big 'Q' for the log-likelihood.
 * This value needs to be aggregated to the log-likelihood
 */
double computeQ(EMContext *ctx)
{
    Matrix *probabilities = &ctx->probabilities; // Get the probabilities matrix
    double *q = ctx->q;                          // Get the q array
    Matrix *W = &ctx->W;                         // Get the W matrix
    IntMatrix *intW = &ctx->intW;                // Get the W matrix in integer form
    // double thirdTerm = underflowSum(q);
    // double total = -thirdTerm;
    double total = 0;
    double borrar = 0;
    double borrar2 = 0;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            int w_bg = MATRIX_AT_PTR(intW, b, g);
            borrar += w_bg == 0 ? 0 : ctx->logGamma[w_bg];
            // total += w_bg == 0 ? 0 : logGammaArr[w_bg]; // Second term
            double qsum = 0;
            double firstTerm = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double q_bgc = Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES);
                double p_gc = MATRIX_AT_PTR(probabilities, g, c);
                firstTerm += (p_gc == 0.0 || q_bgc == 0.0) ? 0 : q_bgc * log(MATRIX_AT_PTR(probabilities, g, c));
            }
            // First term
            total += firstTerm * w_bg;
            borrar2 += firstTerm * w_bg;
        }
    }
    // Rprintf("----------\nQ: %f + %f + %f\n", borrar2, borrar, -thirdTerm);

    return total;
}

void computeQHitAndRun(EMContext *ctx, QMethodInput params, double *ll)
{
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    IntMatrix *intX = &ctx->intX;
    IntMatrix *intW = &ctx->intW;
    double *q = ctx->q;
    Matrix *probabilities = &ctx->probabilities;
    // ---...--- //

    // ---- Use a static assignment since the workload is even between threads ----

    *ll = 0;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        OmegaSet *currentSet = ctx->omegaset[b];
        double *multiplicationValues = (currentSet->size <= 10000) ? (double[10000]){1}                // Stack
                                                                   : Calloc(currentSet->size, double); // Heap
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            double W_bg = MATRIX_AT_PTR(W, b, g);
            if (W_bg == 0)
            {
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                }
                continue;
            }

            // --- Precompute multiplicationValues for this (b, g) combination ---
            double max = -DBL_MAX;
            for (size_t s = 0; s < currentSet->size; s++)
            {                                                  // --- For each sample given a group and a ballot box
                IntMatrix currentMatrix = currentSet->data[s]; // Z
                double a_i = logarithmicProduct(ctx, b, s) +
                             ctx->multinomial[b][s]; // $\sum_{g\in G}\sum_{c\in C}z_{bgc}*\log(p_{gc})+$
                multiplicationValues[s] = a_i;

                max = multiplicationValues[s] > max ? multiplicationValues[s] : max;
            }
            // ---- Shift the values by the maximum
            for (size_t s = 0; s < currentSet->size; s++)
            {
                multiplicationValues[s] -= max;
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box
                // ---- Obtain the summatory over all of the values ---- //
                double sum_exp_den = 0.0;
                double sum_exp_num = 0.0;
                double v = 0;
                for (int i = 0; i < currentSet->size; i++)
                {                                     // --- For each sample
                    v = exp(multiplicationValues[i]); // exp(ls - m)
                    sum_exp_num += v;
                    IntMatrix currentMatrix = currentSet->data[i];
                    sum_exp_den += (MATRIX_AT(currentMatrix, g, c) / W_bg) * v;
                }
                double log_sum_exp_num = log(sum_exp_num);
                double result = exp(log(sum_exp_den) - log_sum_exp_num);
                Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) =
                    !isnan(result) && !isinf(result) ? result : 0;
                // ---...--- //
            } // --- End candidate loop
        } // --- End group loop

        // ---- Calculate the log-likelihood ---- //
        double sum_exp_num = 0.0;
        double v = 0;
        for (int i = 0; i < currentSet->size; i++)
        { // --- For each sample
            v = exp(multiplicationValues[i]);
            sum_exp_num += v;
        }
        for (int i = 0; i < currentSet->size; i++)
        { // --- For each sample
            double val = exp(multiplicationValues[i]) / sum_exp_num;
            *ll -= val * log(currentSet->counts[i] * val);
            *ll += ctx->Qconstant[b][i] * val; // New term
        }
        // ---...--- //

        // ---- Free allocated memory ---- //
        if (currentSet->size > 10000)
            Free(multiplicationValues);
        // ---...--- //
    } // --- End ballot box loop

    // Calculo Q
    double toprint = computeQ(ctx);

    // *ll += computeQ(array2, probabilities);
    *ll += toprint;
}
