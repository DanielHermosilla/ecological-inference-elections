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

#include "multivariate-pdf.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

typedef struct
{
    Matrix muR;      // G x (C-1)
    Matrix **sigma;  // array G of (C-1 x C-1) matrices
    double *feature; // size C (temp: X[,b])
    double *muG;     // size C-1 (temp: muR[g,])
    double *QC;      // size C (per-group numerators)
    double *maha;    // size G*C (mahalanobis distances, row-major)
} Arena;

// Create and initialize an Arena
static Arena Arena_init(int G, int C)
{
    Arena A = {0};
    A.muR = createMatrix(G, C - 1);

    A.sigma = (Matrix **)Calloc(G, Matrix *);
    for (int g = 0; g < G; ++g)
    {
        A.sigma[g] = (Matrix *)Calloc(1, Matrix);
        *(A.sigma[g]) = createMatrix(C - 1, C - 1);
    }

    A.feature = (double *)Calloc(C, double);
    A.muG = (double *)Calloc(C - 1, double);
    A.QC = (double *)Calloc(C, double);
    A.maha = (double *)Calloc((size_t)G * (size_t)C, double);

    return A;
}

// Frees an Arena
static void Arena_free(Arena *A, int G)
{
    if (!A)
        return;
    for (int g = 0; g < G; ++g)
    {
        if (A->sigma && A->sigma[g])
        {
            freeMatrix(A->sigma[g]);
            Free(A->sigma[g]);
        }
    }
    if (A->sigma)
        Free(A->sigma);
    freeMatrix(&A->muR);
    if (A->feature)
        Free(A->feature);
    if (A->muG)
        Free(A->muG);
    if (A->QC)
        Free(A->QC);
    if (A->maha)
        Free(A->maha);
    memset(A, 0, sizeof(*A));
}

static inline void getColumn_into(const Matrix *M, int col, double *out)
{
    for (int r = 0; r < M->rows; r++)
        out[r] = MATRIX_AT_PTR(M, r, col);
}

static inline void getRow_into(const Matrix *M, int row, double *out)
{
    for (int c = 0; c < M->cols; c++)
        out[c] = MATRIX_AT_PTR(M, row, c);
}

/**
 * @brief Computes the `q` values for a given ballot box.
 *
 * Given a ballot box index, probabilities and the reduced version (with C-1 candidates) of the probabilities matrix, it
 * calculates the `q` values in a flattened way
 *
 * @param[in] b. The index of the ballot box
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 * @param[in] *probabilitiesReduced. A pointer towards the reduced probabilities matrix.
 *
 * @return A (g x c) matrix with the values of `q` according the candidate and group index.
 *
 */
void computeQforABallot(EMContext *ctx, int b, const Matrix *probabilities, const Matrix *probabilitiesReduced,
                        double *ll, QMethodInput params, Arena *A)
{
    const int G = (int)ctx->G;
    const int C = (int)ctx->C;
    Matrix *X = &ctx->X;

    // ---- Fill muR and sigma for this ballot ---- //
    getAverageConditional(ctx, b, probabilitiesReduced, &A->muR, A->sigma);

    // Invert sigmas
    for (int g = 0; g < G; ++g)
        inverseSymmetricPositiveMatrix(A->sigma[g]);

    // ---- Compute determinant (for loglikelihood normalization) ---- //
    double normalizeConstant = 1.0;
    if (params.computeLL)
    {
        double det = 1.0;
        for (int c = 0; c < C - 1; ++c)
            det *= MATRIX_AT_PTR(A->sigma[0], c, c);
        det = 1.0 / (det * det);
        normalizeConstant = R_pow(R_pow_di(M_2_PI, C - 1) * det, 0.5);
    }

    // ---- Feature vector (candidate results) ----
    getColumn_into(X, b, A->feature);

    // ---- Mahalanobis per group ----
    for (int g = 0; g < G; ++g)
    { // --- For each group
        getRow_into(&A->muR, g, A->muG);
        double *dst = &A->maha[(size_t)g * (size_t)C];
        // --- Mahalanobis distance for each candidate given a group, accounting the mean
        getMahanalobisDist(A->feature, A->muG, A->sigma[g], dst, C - 1, false);
    }
    // ---- build q’s directly into ctx->q ----
    for (int g = 0; g < G; ++g)
    { // --- For each group
        double *ma = &A->maha[(size_t)g * (size_t)C];

        // ---- Work in log-space to avoid underflow ---- //
        double logw[C]; // log of the unnormalized weights numerator
        double logw_max = -INFINITY;

        for (int c = 0; c < C; ++c)
        {
            double prior = MATRIX_AT_PTR(probabilities, g, c);
            double logP = (prior > 0.0) ? log(prior) : -INFINITY;
            logw[c] = -0.5 * ma[c] + logP; // numerator in log-space
                                           // -0.5 * ma[c] + log(P) <=> exp(-0.5 * ma[c]) * P
            if (isfinite(logw[c]) && logw[c] > logw_max)
                logw_max = logw[c];
        }

        double den = 0.0;
        for (int c = 0; c < C; ++c)
        {
            double val = isfinite(logw[c]) ? exp(logw[c] - logw_max) : 0.0; // Here we go back from log-space
            A->QC[c] = val;
            den += val;
        }

        // ---- log-likelihood contribution ---- //
        if (g == 0 && params.computeLL && den > 0.0 && isfinite(logw_max))
        {
            double logden = logw_max + log(den);
            *ll += logden * log(normalizeConstant);
        }

        // ---- Normalize and store q ---- //
        for (int c = 0; c < C; ++c)
        {
            double qgc = (den > 0.0) ? (A->QC[c] / den) : 0.0;
            Q_3D(ctx->q, b, g, c, G, C) = qgc;
        }
    }
}

/**
 * @brief Computes the `q` values for all the ballot boxes given a probability matrix. Uses the Multivariate PDF method.
 *
 * Given a probability matrix with, it returns a flattened array with estimations of the conditional probability. The
 * array can be accesed with the macro `Q_3D` (it's a flattened tensor).
 *
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 *
 * @return A pointer towards the flattened tensor.
 *
 */
void computeQMultivariatePDF(EMContext *ctx, QMethodInput params, double *ll)
{
    *ll = 0.0;

    Matrix *probabilities = &ctx->probabilities;
    Matrix probabilitiesReduced = removeLastColumn(probabilities);

    // ---- Make only one big allocation ---- //
    Arena A = Arena_init((int)ctx->G, (int)ctx->C);

    for (uint32_t b = 0; b < ctx->B; ++b)
    {
        computeQforABallot(ctx, (int)b, probabilities, &probabilitiesReduced, ll, params, &A);
    }

    // ---- Free the arena ---- //
    Arena_free(&A, (int)ctx->G);
    freeMatrix(&probabilitiesReduced);

    if (isnan(*ll) || isinf(*ll))
        *ll = 0.0;
}
