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

#include "metropolis.h"
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

static const double EPS = 1e-300;

static int lessThanColRow(IntMatrix mat, int b, int g, int c, int candidateVotes, int groupVotes)
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

static IntMatrix startingPoint3(EMContext *ctx, int b)
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

static void allocateRandoms(int M, int S, uint8_t **c1, uint8_t **c2, uint8_t **g1, uint8_t **g2, double **MS)
{
    uint32_t size = M * S;
    // Allocate memory correctly
    *c1 = (uint8_t *)Calloc(size, uint8_t);
    *c2 = (uint8_t *)Calloc(size, uint8_t);
    *g1 = (uint8_t *)Calloc(size, uint8_t);
    *g2 = (uint8_t *)Calloc(size, uint8_t);
    *MS = (double *)Calloc(size, double);

    GetRNGstate(); // Ensure R's RNG is properly initialized
                   // Fill arrays with random indices
    int allow_repeat = (TOTAL_CANDIDATES <= 1 || TOTAL_GROUPS <= 1);

    for (int i = 0; i < size; i++)
    {
        if (i % 400 == 0) // Checks condition every 400 iterations
            R_CheckUserInterrupt();
        (*MS)[i] = (double)(unif_rand()); // The uniform draw
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
 * Obtain the `q` values that are constant within iterations.
 */
void qMid(EMContext *ctx)
{
    Matrix *logProbabilities = &ctx->metropolisProbability;
    // bMatrix *qMetropolis = Calloc(TOTAL_BALLOTS * ctx->omegaset[0]->size, Matrix);
    Matrix qMetropolis = createMatrix(TOTAL_BALLOTS, ctx->omegaset[0]->size);
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        int size = ctx->omegaset[b]->size;
        for (int s = 0; s < size; s++)
        {
            IntMatrix sample = ctx->omegaset[b]->data[s]; // Z_b
            double a = 0;
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                for (int c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    int elementZ = MATRIX_AT(sample, g, c); // z_bgc
                    a += elementZ * MATRIX_AT_PTR(logProbabilities, g, c);
                }
            }
            MATRIX_AT(qMetropolis, b, s) = a; // q_bgs
        }
    }
    ctx->qMetropolis = copMatrix(&qMetropolis);
}

void calculateLogP(EMContext *ctx)
{
    Matrix *oldProbabilities = &ctx->metropolisProbability;
    for (int g = 0; g < TOTAL_GROUPS; g++)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double oldValue = MATRIX_AT_PTR(oldProbabilities, g, c);
            MATRIX_AT_PTR(oldProbabilities, g, c) = oldValue != 0 ? log(MATRIX_AT_PTR(oldProbabilities, g, c)) : -1e10;
        }
    }
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
void generateOmegaSetMetropolis(EMContext *ctx, int M, int S, int burnInSteps)
{
    // ---- Allocate memory for the `b` index ----
    if (ctx->omegaset != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            if (ctx->omegaset[b] != NULL)
            {
                for (int s = 1; s < ctx->omegaset[b]->size; s++)
                {
                    freeMatrixInt(&ctx->omegaset[b]->data[s]);
                }
                // Free(ctx->omegaset[b]->data);
                // Free(ctx->omegaset[b]);
            }
        }
        Free(ctx->omegaset);
    }
    int iteration = ctx->iteration;
    Matrix *probabilities = &ctx->probabilities;           // Get the probabilities matrix
    ctx->metropolisProbability = copMatrix(probabilities); // Copy the probabilities matrix for the metropolis method
    // calculateLogP(ctx); // Calculate the logarithm of the probabilities
    ctx->omegaset = Calloc(TOTAL_BALLOTS, OmegaSet *);
    uint8_t *c1 = NULL;
    uint8_t *c2 = NULL;
    uint8_t *g1 = NULL;
    uint8_t *g2 = NULL;
    double *MS = NULL;

    uint32_t arraySize = M * S;

    allocateRandoms(M, S, &c1, &c2, &g1, &g2, &MS);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases
    int a = 0;

// ---- Perform the main iterations ---- //
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {   // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
        int k = 0;
        if (ctx->omegaset[b] == NULL)
            ctx->omegaset[b] = Calloc(1, OmegaSet);
        if (ctx->omegaset[b]->data == NULL)
            Free(ctx->omegaset[b]->data);
        // ---- Set the ballot box index and size ---- //
        ctx->omegaset[b]->data = Calloc(S, IntMatrix);
        ctx->omegaset[b]->b = b;
        ctx->omegaset[b]->size = S;
        // ---...--- //
        int ballotShift = floor(((double)b / TOTAL_BALLOTS) * (arraySize));

        // M_save =
        int Mactual = iteration == 0 ? burnInSteps : M;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                    steppingZ = startingPoint3(ctx, b);
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                // ----...---- //
            }
            else
            {
                steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            }

            // for (int m = 0; m < Mactual; m++)
            for (int m = 0; m < Mactual; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (arraySize);
                uint8_t randomCDraw = c1[shiftIndex];
                uint8_t randomCDraw2 = c2[shiftIndex];
                uint8_t randomGDraw = g1[shiftIndex];
                uint8_t randomGDraw2 = g2[shiftIndex];

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                    continue;
                // ---...--- //
                double transitionProbNum = (MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) + 1) *
                                           (MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) + 1) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw2);

                double transitionProbDen = (MATRIX_AT(steppingZ, randomGDraw, randomCDraw)) *
                                           (MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2)) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw2) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw);

                double prob = transitionProbDen / transitionProbNum;

                if (MS[shiftIndex] < prob)
                {
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += 1;
                }
                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            a += k;
            // ---...--- //
        } // --- End the sample loop
    } // --- End the ballot box loop
    Rprintf("Did %d movements in total, meaning %.4f movements per ballot box.\n", a, (double)a / TOTAL_BALLOTS);
    calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

static double computeQ(EMContext *ctx)
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

    return total;
}

void computeQhastingIteration(EMContext *ctx, double *ll)
{
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    IntMatrix *intX = &ctx->intX;
    IntMatrix *intW = &ctx->intW;
    double *q = ctx->q;
    Matrix *probabilities = &ctx->probabilities;

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        OmegaSet *currentSet = ctx->omegaset[b];
        double *multiplicationValues = (currentSet->size <= 10000) ? (double[10000]){1}                // Stack
                                                                   : Calloc(currentSet->size, double); // Heap
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            int W_bg = MATRIX_AT_PTR(intW, b, g);
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box
                // ---- Obtain the summatory over all of the values ---- //
                int num = 0;
                for (int s = 0; s < currentSet->size; s++)
                { // --- For each sample
                    IntMatrix currentMatrix = currentSet->data[s];
                    num += (MATRIX_AT(currentMatrix, g, c));
                }
                int den = W_bg * currentSet->size;
                double result = den != 0 ? (double)num / (double)den : 0;
                Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = result;
                // ---...--- //
            } // --- End candidate loop
        } // --- End group loop
        if (currentSet->size > 10000)
            Free(multiplicationValues);
    }
}

// Between iterations
void computeQhastingMidIteration(EMContext *ctx, double *ll)
{
    Matrix *W = &ctx->W;
    IntMatrix *intW = &ctx->intW;
    Matrix *Pold = &ctx->metropolisProbability;
    Matrix *Pnew = &ctx->probabilities;
    double *q = ctx->q;
    Matrix const *a = &ctx->qMetropolis;

    *ll = 0.0;

    Matrix logPnew = createMatrix(Pnew->rows, Pnew->cols);
    for (int g = 0; g < TOTAL_GROUPS; g++)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double currentP = MATRIX_AT_PTR(Pnew, g, c);
            MATRIX_AT(logPnew, g, c) = currentP != 0 ? log(currentP) : -1e10; // Avoid log(0)
        }
    }
    Matrix bMatrix = createMatrix(TOTAL_BALLOTS, ctx->omegaset[0]->size);
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        OmegaSet *S = ctx->omegaset[b];
        int Ssz = S->size;

        // build log‐arrays for old and new
        double max_logw = -INFINITY;
        for (int s = 0; s < Ssz; s++)
        {
            IntMatrix sample = ctx->omegaset[b]->data[s]; // Z_b
            double a2 = 0;
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                for (int c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    int elementZ = MATRIX_AT(sample, g, c); // z_bgc
                    a2 += elementZ * MATRIX_AT(logPnew, g, c);
                }
            }
            MATRIX_AT(bMatrix, b, s) = a2 - MATRIX_AT_PTR(a, b, s); // q_bgs
            if (MATRIX_AT(bMatrix, b, s) > max_logw)
                max_logw = MATRIX_AT(bMatrix, b, s);
        }

        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            int w_bg = MATRIX_AT_PTR(intW, b, g);
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {
                if (w_bg == 0)
                {
                    Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = 0;
                    continue;
                }

                double num = 0.0;
                double denom = 0.0;

                for (int s = 0; s < Ssz; s++)
                {
                    double logw = MATRIX_AT(bMatrix, b, s);
                    double w = exp(logw - max_logw); // log-sum-exp trick
                    int z = MATRIX_AT(S->data[s], g, c);
                    num += w * z;
                    denom += w;
                }

                Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = denom != 0 ? num / (denom * w_bg) : 0;
            }
        }
    }
    freeMatrix(&bMatrix);
    freeMatrix(&logPnew);
}

void computeQMetropolis(EMContext *ctx, QMethodInput params, double *ll)
{
    int *SamplingIter = &params.iters;
    int *currentIter = &ctx->iteration;
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    IntMatrix *intX = &ctx->intX;
    IntMatrix *intW = &ctx->intW;
    double *q = ctx->q;
    Matrix *probabilities = &ctx->probabilities;
    // ---...--- //

    *ll = 0;
    if (ctx->iteration % params.iters == 0)
    {
        generateOmegaSetMetropolis(ctx, params.M, params.S, 10000);
        // encode(ctx);
        // preComputeMultinomial(ctx);
        computeQhastingIteration(ctx, ll);
        return;
    }
    else
    {
        computeQhastingMidIteration(ctx, ll);
        return;
    }

    // ---- Calculate the log-likelihood ---- //
    /*
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
     */
    // ---...--- //

    // ---- Free allocated memory ---- //
    // ---...--- //
} // --- End ballot box loop

// Calculo Q
/*
double toprint = computeQ(ctx);
*/
// *ll += computeQ(array2, probabilities);
