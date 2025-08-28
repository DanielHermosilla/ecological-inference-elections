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

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

static const double EPS = 1e-300;
double *oldQ;
double *logQ;

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

static IntMatrix startingPoint4(EMContext *ctx, int b)
{
    Matrix *W = &ctx->W;
    Matrix *X = &ctx->X;
    double *Q = ctx->q;
    double *votes = ctx->predicted_votes;
    IntMatrix toReturn = createMatrixInt(TOTAL_GROUPS, TOTAL_CANDIDATES);

    int count = 0;
    for (int c = 0; c < TOTAL_CANDIDATES; c++)
    {
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            int vote = (int)Q_3D(votes, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
            MATRIX_AT(toReturn, g, c) = vote;
            count += vote;
        }
    }
    int repartir = ctx->ballots_votes[b] - count;
    if (repartir == 0)
        return toReturn;
    else
    {
        double *groupVotes = getRow(W, b);
        double *candidateVotes = getColumn(X, b);
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            for (int g = 0; g < TOTAL_GROUPS; g++)
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
        return toReturn;
    }
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

static void allocateRandoms_weighted(int M, int S, uint8_t **c1, uint8_t **c2, uint8_t **g1, uint8_t **g2, double **MS,
                                     const uint32_t *group_votes,    // tamaño TOTAL_GROUPS
                                     const uint32_t *candidate_votes // tamaño TOTAL_CANDIDATES
)
{
    uint32_t size = M * S;
    *c1 = Calloc(size, uint8_t);
    *c2 = Calloc(size, uint8_t);
    *g1 = Calloc(size, uint8_t);
    *g2 = Calloc(size, uint8_t);
    *MS = Calloc(size, double);

    // 1) Construir pesos acumulados de grupos
    double cum_g[TOTAL_GROUPS];
    cum_g[0] = group_votes[0];
    for (int g = 1; g < TOTAL_GROUPS; g++)
        cum_g[g] = cum_g[g - 1] + group_votes[g];
    double total_g = cum_g[TOTAL_GROUPS - 1];

    // 2) Pesos acumulados de candidatos
    double cum_c[TOTAL_CANDIDATES];
    cum_c[0] = candidate_votes[0];
    for (int c = 1; c < TOTAL_CANDIDATES; c++)
        cum_c[c] = cum_c[c - 1] + candidate_votes[c];
    double total_c = cum_c[TOTAL_CANDIDATES - 1];

    GetRNGstate();
    for (int i = 0; i < size; i++)
    {
        if (i % 400 == 0)
            R_CheckUserInterrupt();

        // muestreo del uniform draw
        (*MS)[i] = unif_rand();

        // muestreo grupo g1
        double u = unif_rand() * total_g;
        int gsel = 0;
        while (u > cum_g[gsel])
            gsel++;
        (*g1)[i] = (uint8_t)gsel;

        // muestreo candidato c1
        double v = unif_rand() * total_c;
        int csel = 0;
        while (v > cum_c[csel])
            csel++;
        (*c1)[i] = (uint8_t)csel;

        // igual para g2,c2 con allow_repeat
        // int allow_repeat = (TOTAL_CANDIDATES <= 1 || TOTAL_GROUPS <= 1);
        int allow_repeat = 1;
        int g2sel, c2sel;
        do
        {
            double u2 = unif_rand() * total_g;
            g2sel = 0;
            while (u2 > cum_g[g2sel])
                g2sel++;

            double v2 = unif_rand() * total_c;
            c2sel = 0;
            while (v2 > cum_c[c2sel])
                c2sel++;
        } while (!allow_repeat && (g2sel == gsel || c2sel == csel));

        (*g2)[i] = (uint8_t)g2sel;
        (*c2)[i] = (uint8_t)c2sel;
    }
    PutRNGstate();
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
                   // int allow_repeat = (TOTAL_CANDIDATES <= 1 || TOTAL_GROUPS <= 1);
    int allow_repeat = 1;

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

void calculateLogQ(EMContext *ctx)
{
    Matrix *oldProbabilities = &ctx->metropolisProbability;
    double *q = ctx->q; // Get the q array

    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {

                double oldValue = MATRIX_AT_PTR(oldProbabilities, g, c);
                Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) =
                    oldValue != 0 ? log(Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES)) : -1e10;
            }
        }
    }
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

static inline int flatIdx(const IntMatrix *m, uint8_t row, uint8_t col)
{
    return col * m->rows + row;
}

// Move data[k] += delta for all k >= start
static void updateFlattenCumsum(IntMatrix *cum, int start, int delta)
{
    int N = cum->rows * cum->cols;
    for (int k = start; k < N; ++k)
    {
        cum->data[k] += delta;
    }
}

void sampleFromCumulative(const IntMatrix *cum, int *out_g, int *out_c, int draw)
{
    int R = cum->rows;
    int C = cum->cols;
    int N = R * C;
    // 1. total = last element
    int total = cum->data[N - 1];
    // 2. Uniform draw in [1..total]
    // int draw = (int)(unif_rand() * total) + 1;

    // 3. Binary search for the first index >= draw
    int lo = 0, hi = N - 1;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if (cum->data[mid] < draw)
            lo = mid + 1;
        else
            hi = mid;
    }
    int idx = lo;

    // 4. Map back to (row, col) in col-major storage
    *out_g = idx % R;
    *out_c = idx / R;
}

IntMatrix getCumulativeVotes(const IntMatrix *votes)
{
    int R = votes->rows, C = votes->cols, N = R * C;
    IntMatrix cum = createMatrixInt(R, C);
    // arr apunta a los datos de votes, también col-major
    int *arrV = votes->data;
    int *arrC = cum.data;
    if (N > 0)
    {
        arrC[0] = arrV[0];
        for (int k = 1; k < N; ++k)
        {
            arrC[k] = arrC[k - 1] + arrV[k];
        }
    }
    return cum;
}

void sampleFromCumulative2(IntMatrix sampleMat, int *out_g, int *out_c, double draw, int Coriginal, int Goriginal)
{
    int G = sampleMat.rows;
    int C = sampleMat.cols;
    int N = C * G;

    IntMatrix Z = copMatrixI(&sampleMat);

    for (int c = 0; c < C; c++)
    {
        for (int g = 0; g < G; g++)
        {
            if (c == Coriginal || Goriginal == g)
                MATRIX_AT(Z, g, c) = 0;
        }
    }

    IntMatrix accumulatedVotes = getCumulativeVotes(&Z);
    int Ncells = accumulatedVotes.rows * accumulatedVotes.cols;
    int total1 = accumulatedVotes.data[Ncells - 1];
    int drawFinal = (int)(draw * total1) + 1;
    sampleFromCumulative(&accumulatedVotes, out_g, out_c, drawFinal);
    freeMatrixInt(&Z);
    freeMatrixInt(&accumulatedVotes);
}

int fzh(uint8_t g1, uint8_t g2, uint8_t c1, uint8_t c2, IntMatrix z, int h, double perc)
{

    int primerodos = MIN(MATRIX_AT(z, g1, c1) - h, MATRIX_AT(z, g2, c2) - h);

    int ultimosdos = MIN(MATRIX_AT(z, g2, c1) + h, MATRIX_AT(z, g1, c2) + h);

    int minimo = (int)floor(MIN(primerodos, ultimosdos) * perc);

    return MAX(1, minimo);
}

void generateOmegaSetMetropolis6(EMContext *ctx, int M, int S, int burnInSteps, char *initial_value)
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
    double *Q = ctx->q;                                    // Get the probabilities matrix
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

    calculateLogP(ctx); // Calculate the logarithm of the probabilities, guarda en metropolisProbability
    allocateRandoms(M, S, &c1, &c2, &g1, &g2, &MS);

    // allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    uint64_t atest = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
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
        uint64_t ktest = 0;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                {
                    if (strcmp(initial_value, "proportional") == 0)
                        steppingZ = startingPoint3(ctx, b);
                    else
                        steppingZ = startingPoint4(ctx, b);
                }
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                // ----...---- //
            }
            else
            {
                steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            }

            IntMatrix accumulatedVotes = getCumulativeVotes(&steppingZ);

            // for (int m = 0; m < Mactual; m++)
            for (int m = 0; m < Mactual; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (arraySize);
                int shiftIndex1 = (s * M + ballotShift + m + 42) % (arraySize);

                int randomGDraw, randomCDraw;
                int randomCDraw2;
                int randomGDraw2;

                int Ncells = accumulatedVotes.rows * accumulatedVotes.cols;
                int total1 = accumulatedVotes.data[Ncells - 1];

                int draw1 = (int)(MS[shiftIndex] * total1) + 1;
                double draw2 = MS[shiftIndex1];

                sampleFromCumulative(&accumulatedVotes, &randomGDraw, &randomCDraw, draw1);
                sampleFromCumulative2(steppingZ, &randomGDraw2, &randomCDraw2, draw2, randomCDraw, randomGDraw);

                if (randomCDraw == randomCDraw2 || randomGDraw == randomGDraw2)
                {
                    Rprintf("Se sampleo c1 = %d, c2 = %d, g1 = %d, g2 = %d\n", randomCDraw, randomCDraw2, randomGDraw,
                            randomGDraw2);
                    continue;
                }

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                {
                    continue;
                }

                int f = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, 0, 0.1);
                int h = (int)(f * MS[shiftIndex1]) + 1;
                int fh = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, h, 0.1);

                if (h > fh)
                    continue;

                Matrix *logProb = &ctx->metropolisProbability;
                double firstLine = h * (-MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw) -
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw2) +
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw) +
                                        MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw2));

                double secondLine = lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw) - 1) -
                                    lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw) - h);

                double thirdLine = lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) - 1) -
                                   lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) - h);

                double fourthLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) + h - 1) +
                                    lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw));

                double fifthLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) + h - 1) +
                                   lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2));

                double sixthLine = log((double)f / (double)fh);

                Matrix W = ctx->W;
                Matrix X = ctx->X;

                int commonTerm = 2 * ctx->total_votes - MATRIX_AT(W, b, randomGDraw) - MATRIX_AT(W, b, randomGDraw2) -
                                 MATRIX_AT(X, randomCDraw, b) - MATRIX_AT(X, randomCDraw2, b);

                double seventhLine = log(2 * h + commonTerm + MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) +
                                         MATRIX_AT(steppingZ, randomGDraw2, randomCDraw));
                double eigthLine =
                    -log(ctx->total_votes - MATRIX_AT(W, b, randomGDraw) - MATRIX_AT(X, randomCDraw2, b) +
                         MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) + h);

                double nineLine = -log(ctx->total_votes - MATRIX_AT(W, b, randomGDraw2) - MATRIX_AT(X, randomCDraw, b) +
                                       MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) + h);

                double tenthLine = -log(commonTerm + MATRIX_AT(steppingZ, randomGDraw, randomCDraw) +
                                        MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2));

                double eleventhLine =
                    log(ctx->total_votes - MATRIX_AT(W, b, randomGDraw) - MATRIX_AT(X, randomCDraw, b) +
                        MATRIX_AT(steppingZ, randomGDraw, randomCDraw));

                double twelveLine =
                    log(ctx->total_votes - MATRIX_AT(W, b, randomGDraw2) - MATRIX_AT(X, randomCDraw2, b) +
                        MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2));

                double total = firstLine + secondLine + thirdLine + fourthLine + fifthLine + sixthLine + seventhLine +
                               eigthLine + nineLine + tenthLine + eleventhLine + twelveLine;

                // Rprintf("%.4f\n", total);
                // if (randomCDraw == 8 && randomGDraw == 4)
                // Rprintf("%.4f\n", exp(total));
                if (log(MS[shiftIndex]) < total)
                {
                    ktest += 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= h;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += h;

                    int s1 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw);
                    int s2 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw2);
                    int s3 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw2);
                    int s4 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw);

                    updateFlattenCumsum(&accumulatedVotes, s1, -h);
                    updateFlattenCumsum(&accumulatedVotes, s2, -h);
                    updateFlattenCumsum(&accumulatedVotes, s3, h);
                    updateFlattenCumsum(&accumulatedVotes, s4, h);
                }

                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
            atest / (double)(TOTAL_BALLOTS * S * M));
    // calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void generateOmegaSetMetropolis5(EMContext *ctx, int M, int S, int burnInSteps, char *initial_value)
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
    double *Q = ctx->q;                                    // Get the probabilities matrix
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

    calculateLogP(ctx); // Calculate the logarithm of the probabilities, guarda en metropolisProbability
    allocateRandoms(M, S, &c1, &c2, &g1, &g2, &MS);

    // allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    uint64_t atest = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
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
        uint64_t ktest = 0;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                {
                    if (strcmp(initial_value, "proportional") == 0)
                        steppingZ = startingPoint3(ctx, b);
                    else
                        steppingZ = startingPoint4(ctx, b);
                }
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                // ----...---- //
            }
            else
            {
                steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            }

            IntMatrix accumulatedVotes = getCumulativeVotes(&steppingZ);

            // for (int m = 0; m < Mactual; m++)
            for (int m = 0; m < Mactual; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (arraySize);
                int shiftIndex1 = (s * M + ballotShift + m + 42) % (arraySize);

                int randomGDraw, randomCDraw;

                int Ncells = accumulatedVotes.rows * accumulatedVotes.cols;
                int total1 = accumulatedVotes.data[Ncells - 1];
                if (total1 <= 0)
                {
                    // sin masa para samplear: sáltate este m o elige un fallback
                    continue;
                }

                int draw1 = (int)(MS[shiftIndex] * total1) + 1;
                int draw2 = (int)(MS[shiftIndex1] * total1) + 1;
                // clamp por seguridad
                if (draw1 > total1)
                    draw1 = total1;
                if (draw2 > total1)
                    draw2 = total1;

                sampleFromCumulative(&accumulatedVotes, &randomGDraw, &randomCDraw, draw1);
                // uint8_t randomCDraw = c1[shiftIndex];
                //  uint8_t randomCDraw2 = c2[shiftIndex];
                //  uint8_t randomGDraw = g1[shiftIndex];
                //  uint8_t randomGDraw2 = g2[shiftIndex]
                int randomCDraw2;
                int randomGDraw2;
                sampleFromCumulative(&accumulatedVotes, &randomGDraw2, &randomCDraw2, draw2);

                // Rprintf("Se sampleo c1 = %d, c2 = %d, g1 = %d, g2 = %d\n", randomCDraw, randomCDraw2, randomGDraw,
                //        randomGDraw2);
                if (randomCDraw == randomCDraw2 || randomGDraw == randomGDraw2)
                    continue;

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                {
                    continue;
                }

                int f = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, 0, 0.1);
                int h = (int)(f * MS[shiftIndex1]) + 1;
                int fh = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, h, 0.1);

                if (h > fh)
                    continue;

                Matrix *logProb = &ctx->metropolisProbability;
                double firstLine = h * (-MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw) -
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw2) +
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw) +
                                        MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw2));

                double secondLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw) - h);

                double thirdLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) - h);

                double fourthLine = -lgamma1p(h) + lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw));

                double fifthLine = -lgamma1p(h) + lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2));

                double sixthLine = log((double)f / (double)fh);

                double total = firstLine + secondLine + thirdLine + fourthLine + fifthLine + sixthLine;

                // Rprintf("%.4f\n", total);
                // Rprintf("%.4f\n", exp(total));
                if (log(MS[shiftIndex]) < total)
                {
                    ktest += 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= h;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += h;

                    int s1 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw);
                    int s2 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw2);
                    int s3 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw2);
                    int s4 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw);

                    updateFlattenCumsum(&accumulatedVotes, s1, -h);
                    updateFlattenCumsum(&accumulatedVotes, s2, -h);
                    updateFlattenCumsum(&accumulatedVotes, s3, h);
                    updateFlattenCumsum(&accumulatedVotes, s4, h);
                }

                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
            atest / (double)(TOTAL_BALLOTS * S * M));
    // calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void generateOmegaSetMetropolis4(EMContext *ctx, int M, int S, int burnInSteps, char *initial_value,
                                 char *sampling_method)
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
    double *Q = ctx->q;                                    // Get the probabilities matrix
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

    if (strcmp(sampling_method, "uniform") == 0)
        allocateRandoms(M, S, &c1, &c2, &g1, &g2, &MS);
    else
        allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);

    calculateLogP(ctx); // Calculate the logarithm of the probabilities, guarda en metropolisProbability

    // allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    uint64_t atest = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
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
        uint64_t ktest = 0;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                    if (strcmp(initial_value, "proportional") == 0)
                        steppingZ = startingPoint3(ctx, b);
                    else
                        steppingZ = startingPoint4(ctx, b);
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                // ----...---- //
            }
            else
            {
                steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            }

            // IntMatrix accumulatedVotes = getCumulativeVotes(&steppingZ);

            // for (int m = 0; m < Mactual; m++)
            for (int m = 0; m < Mactual; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (arraySize);
                int shiftIndex1 = (s * M + ballotShift + m + 42) % (arraySize);
                uint8_t randomCDraw = c1[shiftIndex];
                uint8_t randomCDraw2 = c2[shiftIndex];
                uint8_t randomGDraw = g1[shiftIndex];
                uint8_t randomGDraw2 = g2[shiftIndex];

                // int randomCDraw2;
                // int randomGDraw2;
                // sampleFromCumulative(&accumulatedVotes, &randomGDraw2, &randomCDraw2);

                if (c1[shiftIndex] == c2[shiftIndex] || g1[shiftIndex] == g2[shiftIndex])
                    continue;

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                {
                    continue;
                }

                int f = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, 0, 0.1);
                int h = (int)(f * MS[shiftIndex1]) + 1;
                int fh = fzh(randomGDraw, randomGDraw2, randomCDraw, randomCDraw2, steppingZ, h, 0.1);

                if (h > fh)
                    continue;

                Matrix *logProb = &ctx->metropolisProbability;
                double firstLine = h * (-MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw) -
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw2) +
                                        MATRIX_AT_PTR(logProb, randomGDraw2, randomCDraw) +
                                        MATRIX_AT_PTR(logProb, randomGDraw, randomCDraw2));

                double secondLine = lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw)) -
                                    lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw) - h);

                double thirdLine = lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2)) -
                                   lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) - h);

                double fourthLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) + h) +
                                    lgamma1p(MATRIX_AT(steppingZ, randomGDraw2, randomCDraw));

                double fifthLine = -lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) + h) +
                                   lgamma1p(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2));

                double sixthLine = log((double)f / (double)fh);

                double total = firstLine + secondLine + thirdLine + fourthLine + fifthLine + sixthLine;

                {
                    int minz1 = MIN(MATRIX_AT(steppingZ, randomGDraw, randomCDraw),
                                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2));
                    int minz2 = MIN(MATRIX_AT(steppingZ, randomGDraw, randomCDraw2),
                                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw));
                    int minimo = MIN(minz1, minz2);
                    double pct = (double)f / (double)minimo;
                    /*
                    Rprintf("minimo = %d, Z1: %d, Z2: %d:, Z3: %d, Z4: %d, C1: %d, C2: %d, G1: %d, G2: %d, h: %d, f: "
                            "%d, fh: %d, f_pct = %.4f, ratio = "
                            "%.4f\n",
                            minimo, MATRIX_AT(steppingZ, randomGDraw, randomCDraw),
                            MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2),
                            MATRIX_AT(steppingZ, randomGDraw, randomCDraw2),
                            MATRIX_AT(steppingZ, randomGDraw2, randomCDraw), randomCDraw, randomCDraw2, randomGDraw,
                            randomGDraw2, h, f, fh, pct, exp(total));
                    */
                }

                // Rprintf("%.4f\n", total);
                // Rprintf("%.4f\n", exp(total));
                if (log(MS[shiftIndex]) < total)
                {
                    // Rprintf("El ratio uniforme es %.4f\n", MS[shiftIndex]);
                    ktest += 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= h;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += h;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += h;
                    /*
                    if (b == 25 &&
                            ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)) ||
                        (b == 25 &&
                         ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2))))

                    {
                        printMatrixInt(&steppingZ);
                        Rprintf("-------\n\n");
                    }
                    */
                }

                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
            atest / (double)(TOTAL_BALLOTS * S * M));
    // calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void generateOmegaSetMetropolis3(EMContext *ctx, int M, int S, int burnInSteps)
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
    // double *Q = ctx->q;                                    // Get the probabilities matrix
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

    allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);

    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    IntMatrix **votingProbabilities = Calloc(TOTAL_BALLOTS, IntMatrix *);
    uint64_t atest = 0;
    // #ifdef _OPENMP
    // #pragma omp parallel for
    // #endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
        if (ctx->omegaset[b] == NULL)
            ctx->omegaset[b] = Calloc(1, OmegaSet);
        if (ctx->omegaset[b]->data == NULL)
            Free(ctx->omegaset[b]->data);
        // ---- Set the ballot box index and size ---- //
        ctx->omegaset[b]->data = Calloc(S, IntMatrix);
        ctx->omegaset[b]->b = b;
        ctx->omegaset[b]->size = S;
        int total_votes = ctx->ballots_votes[b];
        // ---...--- //
        int ballotShift = floor(((double)b / TOTAL_BALLOTS) * (arraySize));

        // M_save =
        int Mactual = iteration == 0 ? burnInSteps : M;
        uint64_t ktest = 0;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                    steppingZ = startingPoint4(ctx, b);
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                if (b == 25)
                {
                    // Rprintf("Partiendo con\n");
                    //  printMatrixInt(&steppingZ);
                }
                // ----...---- //
            }
            else
            {
                steppingZ = copMatrixI(&ctx->omegaset[b]->data[s - 1]);
            }

            IntMatrix accumulatedVotes = getCumulativeVotes(&steppingZ);

            // for (int m = 0; m < Mactual; m++)
            for (int m = 0; m < Mactual; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (arraySize);
                int randomGDraw, randomCDraw;
                // sampleFromCumulative(&accumulatedVotes, &randomGDraw, &randomCDraw);
                // uint8_t randomCDraw = c1[shiftIndex];
                //  uint8_t randomCDraw2 = c2[shiftIndex];
                //  uint8_t randomGDraw = g1[shiftIndex];
                //  uint8_t randomGDraw2 = g2[shiftIndex]
                int randomCDraw2;
                int randomGDraw2;
                // sampleFromCumulative(&accumulatedVotes, &randomGDraw2, &randomCDraw2);

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (b == 25 && ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)))
                {
                    Rprintf("Se tiene que c1 = %d, c2 = %d, g1 = %d y g2 = %d. ACá SE DEBERíA QUITAR VOTOS\n",
                            randomCDraw, randomCDraw2, randomGDraw, randomGDraw2);
                }
                if (b == 25 && ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2)))
                {
                    Rprintf("Se tiene que c1 = %d, c2 = %d, g1 = %d y g2 = %d, ACá SE DEBERIA AÑADIR VOTOS\n",
                            randomCDraw, randomCDraw2, randomGDraw, randomGDraw2);
                }

                if (firstSubstraction <= 0 || secondSubstraction <= 0 || randomCDraw == randomCDraw2 ||
                    randomGDraw == randomGDraw2)
                {

                    if (b == 25 && ((randomCDraw == randomCDraw2) || (randomGDraw == randomGDraw2)))
                    {
                        // Rprintf("Se escoge la misma columna y/o fila. Se tiene g1 = %d, g2 = %d, c1 = %d y c2 =
                        // %d\n",
                        //        randomGDraw, randomGDraw2, randomCDraw, randomCDraw2);
                        // printMatrixInt(&steppingZ);
                    }

                    if (b == 25 &&
                        ((randomCDraw == 2 && randomGDraw == 0) || (randomCDraw2 == 2 && randomGDraw2 == 0)) &&
                        (randomCDraw == randomCDraw2 || randomGDraw == randomGDraw2))
                    {
                        Rprintf("Se elige la misma columna y/o fila dado que estamos en (3,1)\n");
                        continue;
                    }
                    if (b == 25 && ((randomCDraw == 2 && randomGDraw == 0) || (randomCDraw2 == 2 && randomGDraw2 == 0)))
                    {
                        Rprintf("Resta negativa, seguir\n\n");
                        // printMatrixInt(&steppingZ);
                    }

                    continue;
                }

                // ---...--- //
                double transitionProbNum = MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw2);

                double transitionProbDen = MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw2) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw);

                double prob = transitionProbDen / transitionProbNum;

                if (b == 25 && ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)))
                    Rprintf("El ratio de transición para QUITAR votos es %.4f\n", prob);

                if (b == 25 && ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2)))
                    Rprintf("El ratio de transición para  aÑADIR votos es %.4f\n", prob);

                if (MS[shiftIndex] < prob)
                {
                    ktest += 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += 1;
                    if (b == 25 &&
                            ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)) ||
                        (b == 25 &&
                         ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2))))

                    {
                        printMatrixInt(&steppingZ);
                        Rprintf("-------\n\n");
                    }
                    int s1 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw);
                    int s2 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw2);
                    int s3 = flatIdx(&accumulatedVotes, randomGDraw, randomCDraw2);
                    int s4 = flatIdx(&accumulatedVotes, randomGDraw2, randomCDraw);

                    // 2) Push deltas into the cumsum array
                    updateFlattenCumsum(&accumulatedVotes, s1, -1);
                    updateFlattenCumsum(&accumulatedVotes, s2, -1);
                    updateFlattenCumsum(&accumulatedVotes, s3, +1);
                    updateFlattenCumsum(&accumulatedVotes, s4, +1);
                }

                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
            atest / (double)(TOTAL_BALLOTS * S * M));
    calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void generateOmegaSetMetropolis2(EMContext *ctx, int M, int S, int burnInSteps)
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
    double *Q = ctx->q;                                    // Get the probabilities matrix
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
    // allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    uint64_t atest = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
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
        uint64_t ktest = 0;

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            IntMatrix steppingZ;
            if (s == 0)
            {
                // ---- The `base` element used as a starting point ----
                if (ctx->omegaset[b]->data[0].data == NULL)
                    steppingZ = startingPoint4(ctx, b);
                else
                    steppingZ = ctx->omegaset[b]->data[0];
                if (b == 25)
                {
                    // Rprintf("Partiendo con\n");
                    //  printMatrixInt(&steppingZ);
                }
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

                if (c1[shiftIndex] == c2[shiftIndex] || g1[shiftIndex] == g2[shiftIndex])
                    continue;

                // ---- Check non negativity condition ---- //
                int firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                int secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                /*
                if (b == 25 && ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)))
                {
                    Rprintf("Se tiene que c1 = %d, c2 = %d, g1 = %d y g2 = %d. ACá SE DEBERíA QUITAR VOTOS\n",
                            randomCDraw, randomCDraw2, randomGDraw, randomGDraw2);
                }
                if (b == 25 && ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2)))
                {
                    Rprintf("Se tiene que c1 = %d, c2 = %d, g1 = %d y g2 = %d, ACá SE DEBERIA AÑADIR VOTOS\n",
                            randomCDraw, randomCDraw2, randomGDraw, randomGDraw2);
                }
*/
                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                {
                    /*
                    if (b == 25 && ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)))
                    {
                        Rprintf("La resta da negativa\n");
                        printMatrixInt(&steppingZ);
                    }
                    */
                    if (b == 25 && ((randomCDraw == 2 && randomGDraw == 0) || (randomCDraw2 == 2 && randomGDraw2 == 0)))
                    {
                        // Rprintf("Resta negativa, seguir\n\n");
                        // printMatrixInt(&steppingZ);
                    }

                    continue;
                }
                // ---...--- //
                double transitionProbNum = (MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) + 1) *
                                           (MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) + 1) *
                                           // Q_3D(Q, b, randomGDraw, randomCDraw, TOTAL_GROUPS, TOTAL_CANDIDATES) *
                                           // Q_3D(Q, b, randomGDraw2, randomCDraw2, TOTAL_GROUPS, TOTAL_CANDIDATES);
                                           MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw2);

                double transitionProbDen = (MATRIX_AT(steppingZ, randomGDraw, randomCDraw)) *
                                           (MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2)) *
                                           // Q_3D(Q, b, randomGDraw, randomCDraw2, TOTAL_GROUPS, TOTAL_CANDIDATES) *
                                           // Q_3D(Q, b, randomGDraw2, randomCDraw, TOTAL_GROUPS, TOTAL_CANDIDATES);
                                           MATRIX_AT_PTR(probabilities, randomGDraw, randomCDraw2) *
                                           MATRIX_AT_PTR(probabilities, randomGDraw2, randomCDraw);

                double prob = transitionProbDen / transitionProbNum;

                // if (b == 25 && ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)))
                //  Rprintf("El ratio de transición para QUITAR votos es %.4f\n", prob);

                // if (b == 25 && ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2)))
                // Rprintf("El ratio de transición para  aÑADIR votos es %.4f\n", prob);

                if (MS[shiftIndex] < prob)
                {
                    // Rprintf("El ratio uniforme es %.4f\n", MS[shiftIndex]);
                    ktest += 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= 1;
                    MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += 1;
                    MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += 1;
                    /*
                    if (b == 25 &&
                            ((randomCDraw == 0 && randomGDraw == 2) || (randomCDraw2 == 0 && randomGDraw2 == 2)) ||
                        (b == 25 &&
                         ((randomCDraw == 0 && randomGDraw2 == 2) || (randomCDraw2 == 0 && randomGDraw == 2))))

                    {
                        printMatrixInt(&steppingZ);
                        Rprintf("-------\n\n");
                    }
                    */
                }

                //  ---...--- //
            } // --- End the step size loop
            Mactual = M;
            // ---- Add the combination to the ctx->omegaset ---- //
            ctx->omegaset[b]->data[s] = steppingZ;
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
            atest / (double)(TOTAL_BALLOTS * S * M));
    calculateLogP(ctx); // Calculate the logarithm of the probabilities
    qMid(ctx);
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
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

    allocateRandoms_weighted(M, S, &c1, &c2, &g1, &g2, &MS, ctx->group_votes, ctx->candidates_votes);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

    // ---- Perform the main iterations ---- //
    uint64_t atest = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        // ---- Allocate memory for the ctx->omegaset ---- //
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
        uint64_t ktest = 0;

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
                    ktest += 1;
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
            // ---...--- //
        } // --- End the sample loop
        atest += ktest;
    } // --- End the ballot box loop
    // Rprintf("Did %llu metropolis steps, equivalent to a %.4f percent \n", atest,
    //      atest / (double)(TOTAL_BALLOTS * S * M));
    // calculateLogP(ctx); // Calculate the logarithm of the probabilities
    // qMid(ctx);
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
    // Rprintf("El primer elemento de q es: %.4f\n", Q_3D(q, 0, 0, 0, TOTAL_GROUPS, TOTAL_CANDIDATES));
    Matrix *probabilities = &ctx->probabilities;
    // Rprintf("Se tiene que el omegaset es de tamaño %d\n", (int)ctx->omegaset[0]->size);

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
                // if (g == 3 && c == 1 && b % 5 == 0)
                // Rprintf("\niteración con sampleo para b = %d !\n", b);
                for (int s = 0; s < currentSet->size; s++)
                { // --- For each sample
                    IntMatrix currentMatrix = currentSet->data[s];
                    num += (MATRIX_AT(currentMatrix, g, c));
                    // if (g == 3 && c == 1 && b % 5 == 0)
                    //    Rprintf("%d, ", MATRIX_AT(currentMatrix, g, c));
                }
                // if (g == 3 && c == 1 && (b % 5 == 0))
                //  Rprintf("\n");
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

void computeQhastingMidIteration2(EMContext *ctx, double *ll)
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

                // if (g == 3 && c == 1 && (b % 5 == 0))
                // Rprintf("\niteración sin sampleo para b = %d!\n", b);
                for (int s = 0; s < Ssz; s++)
                {
                    double logw = MATRIX_AT(bMatrix, b, s);
                    double w = exp(logw - max_logw); // log-sum-exp trick
                    int z = MATRIX_AT(S->data[s], g, c);
                    num += w * z;
                    denom += w;
                    //  if (g == 3 & c == 1 && (b % 5 == 0))
                    //     Rprintf("%d, ", z);
                }
                // if (g == 3 && c == 1 && (b % 5 == 0))
                //  Rprintf("\n");

                Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = denom != 0 ? num / (denom * w_bg) : 0;
            }
        }
    }
    freeMatrix(&bMatrix);
    freeMatrix(&logPnew);
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

                // if (g == 3 && c == 1 && (b % 5 == 0))
                // Rprintf("\niteración sin sampleo para b = %d!\n", b);
                for (int s = 0; s < Ssz; s++)
                {
                    double logw = MATRIX_AT(bMatrix, b, s);
                    double w = exp(logw - max_logw); // log-sum-exp trick
                    int z = MATRIX_AT(S->data[s], g, c);
                    num += w * z;
                    denom += w;
                    //  if (g == 3 & c == 1 && (b % 5 == 0))
                    //     Rprintf("%d, ", z);
                }
                // if (g == 3 && c == 1 && (b % 5 == 0))
                //  Rprintf("\n");

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

    Rprintf("Running metropolis with %s sampling method and %s starting point\n", params.sampling_method,
            params.initial_value);

    if (ctx->iteration == 0)
    {
        if (strcmp(params.sampling_method, "sample_proportional") == 0)
            generateOmegaSetMetropolis5(ctx, params.M, params.S, params.burnInSteps, params.initial_value);
        else if (strcmp(params.sampling_method, "sample_proportional2") == 0)
            generateOmegaSetMetropolis6(ctx, params.M, params.S, params.burnInSteps, params.initial_value);
        else
            generateOmegaSetMetropolis4(ctx, params.M, params.S, params.burnInSteps, params.initial_value,
                                        params.sampling_method);

        computeQhastingIteration(ctx, ll);
        return;
    }
    if (ctx->iteration % params.iters == 0)
    {
        if (strcmp(params.sampling_method, "sample_proportional") == 0)
            generateOmegaSetMetropolis5(ctx, params.M, params.S, params.burnInSteps, params.initial_value);
        else if (strcmp(params.sampling_method, "sample_proportional2") == 0)
            generateOmegaSetMetropolis6(ctx, params.M, params.S, params.burnInSteps, params.initial_value);
        else
            generateOmegaSetMetropolis4(ctx, params.M, params.S, params.burnInSteps, params.initial_value,
                                        params.sampling_method);
        // encode(ctx);
        // preComputeMultinomial(ctx);
        computeQhastingIteration(ctx, ll);
        // computeQhastingMidIteration(ctx, ll);
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
