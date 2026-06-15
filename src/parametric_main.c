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

#include "parametric_main.h"
#include "LP.h"
#include "globals.h"
#include "multivariate-cdf.h"
#include "multivariate-pdf.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <R_ext/Utils.h> // for R_CheckUserInterrupt
#include <Rinternals.h>
#include <Rmath.h>
#include <dirent.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#ifdef beta
#undef beta
#endif
#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#ifndef CLOCK_MONOTONIC_RAW
#define CLOCK_MONOTONIC_RAW 4
#endif
#ifndef BLAS_INT
#define BLAS_INT int
#endif
#undef I

typedef struct
{
    int B, G, C, A, D;
    Matrix *prob;  // length B, each G x C
    Matrix S_bc;   // B x C
    Matrix *q_bgc; // length B, each G x C
    Matrix VxA;    // B x (C-1)
    Matrix alpha;
    Matrix beta;
    Matrix grad_alpha; // (C-1) x A
    Matrix grad_beta;  // G x (C-1)
    Matrix H;          // D x D Hessian
    double *gvec;      // length D
    double *vvec;      // length D
    double *ballot_loglik;
} EMBuffers;

void init_EMBuffers(EMBuffers *buf, int B, int G, int Cminus1, int A)
{
    buf->B = B;
    buf->G = G;
    buf->C = Cminus1 + 1;
    buf->A = A;
    buf->D = Cminus1 * A + G * Cminus1;

    // Preallocate probability tensor
    buf->prob = (Matrix *)Calloc(B, Matrix);
    for (int b = 0; b < B; ++b)
    {
        buf->prob[b] = createMatrix(G, buf->C);
    }

    // Preallocate S_bc
    buf->S_bc = createMatrix(B, buf->C);

    // Preallocate q_bgc
    buf->q_bgc = (Matrix *)Calloc(B, Matrix);
    for (int b = 0; b < B; ++b)
    {
        buf->q_bgc[b] = createMatrix(G, buf->C);
    }

    // Preallocate VxA
    buf->VxA = createMatrix(B, Cminus1);

    // Gradients and Hessian

    buf->grad_alpha = createMatrix(Cminus1, A);
    buf->grad_beta = createMatrix(G, Cminus1);
    buf->H = createMatrix(buf->D, buf->D);

    // Vectors
    buf->gvec = (double *)Calloc(buf->D, double);
    buf->vvec = (double *)Calloc(buf->D, double);
    buf->ballot_loglik = (double *)Calloc(B, double);
}

void free_EMBuffers(EMBuffers *buf)
{
    for (int b = 0; b < buf->B; ++b)
    {
        if (buf->prob != NULL)
            freeMatrix(&buf->prob[b]);
        if (buf->q_bgc != NULL)
            freeMatrix(&buf->q_bgc[b]);
    }
    Free(buf->prob);
    Free(buf->q_bgc);

    freeMatrix(&buf->VxA);
    freeMatrix(&buf->S_bc);

    freeMatrix(&buf->H);
    freeMatrix(&buf->grad_alpha);
    freeMatrix(&buf->grad_beta);
    Free(buf->gvec);
    Free(buf->vvec);
    Free(buf->ballot_loglik);
}

// ---- Copy helpers for returning results ----
static Matrix *alloc_matrix_array(int B, int rows, int cols)
{
    Matrix *arr = (Matrix *)Calloc(B, Matrix);
    for (int b = 0; b < B; ++b)
    {
        arr[b] = createMatrix(rows, cols);
    }
    return arr;
}

static void copy_matrix_array(Matrix *dest, const Matrix *src, int B)
{
    for (int b = 0; b < B; ++b)
    {
        size_t n = (size_t)src[b].rows * (size_t)src[b].cols;
        memcpy(dest[b].data, src[b].data, n * sizeof(double));
    }
}

static void compute_expected_outcome(const Matrix *W, Matrix *q_bgc, Matrix *expected, int B, int G, int C)
{
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            double w = MATRIX_AT_PTR(W, b, g);
            for (int c = 0; c < C; ++c)
            {
                MATRIX_AT(expected[b], g, c) = w * MATRIX_AT(q_bgc[b], g, c);
            }
        }
    }
}

static bool hasMismatch(const Matrix *X, const Matrix *W)
{
    int B = X->rows;
    int C = X->cols;
    int G = W->cols;
    for (int b = 0; b < B; b++)
    {
        double sumW = 0.0;
        double sumX = 0.0;
        for (int g = 0; g < G; g++)
            sumW += MATRIX_AT_PTR(W, b, g);
        for (int c = 0; c < C; c++)
            sumX += MATRIX_AT_PTR(X, b, c);
        if (fabs(sumW - sumX) > 1e-6)
            return true;
    }
    return false;
}

static void precomputeScaleFactors(double *scale_factors, const Matrix *X, const Matrix *W)
{
    int B = X->rows;
    int C = X->cols;
    int G = W->cols;
    for (int b = 0; b < B; b++)
    {
        double sum_x = 0.0;
        double sum_w = 0.0;
        for (int c = 0; c < C; c++)
            sum_x += MATRIX_AT_PTR(X, b, c);
        for (int g = 0; g < G; g++)
            sum_w += MATRIX_AT_PTR(W, b, g);
        scale_factors[b] = sum_x / sum_w;
    }
}

static Matrix precomputeNorm(const double *scale_factors, const Matrix *W)
{
    int B = W->rows;
    int G = W->cols;
    Matrix norm = createMatrix(B, G);
    for (int b = 0; b < B; b++)
    {
        double sum = 0.0;
        for (int g = 0; g < G; g++)
        {
            double w = MATRIX_AT_PTR(W, b, g);
            sum += w * w;
        }
        sum *= scale_factors[b];
        for (int g = 0; g < G; g++)
            MATRIX_AT(norm, b, g) = MATRIX_AT_PTR(W, b, g) / sum;
    }
    return norm;
}

static void projectQ(const Matrix *X, const Matrix *W, EMBuffers *buf, const Matrix *norm, const double *scale_factors)
{
    int B = buf->B;
    int G = buf->G;
    int C = buf->C;

    Matrix temp = createMatrix(B, C);
    for (int b = 0; b < B; b++)
    {
        for (int c = 0; c < C; c++)
        {
            double sum = 0.0;
            for (int g = 0; g < G; g++)
            {
                sum += MATRIX_AT(buf->q_bgc[b], g, c) * MATRIX_AT_PTR(W, b, g) * scale_factors[b];
            }
            MATRIX_AT(temp, b, c) = sum;
        }
    }

    for (int b = 0; b < B; b++)
    {
        for (int g = 0; g < G; g++)
        {
            for (int c = 0; c < C; c++)
            {
                double predictedVote = MATRIX_AT(temp, b, c);
                MATRIX_AT(buf->q_bgc[b], g, c) = MATRIX_AT(buf->q_bgc[b], g, c) -
                                                 (predictedVote - MATRIX_AT_PTR(X, b, c)) * MATRIX_AT_PTR(norm, b, g);
            }
        }
    }

    for (int b = 0; b < B; b++)
    {
        for (int g = 0; g < G; g++)
        {
            for (int c = 0; c < C; c++)
            {
                double q = MATRIX_AT(buf->q_bgc[b], g, c);
                if (q < 0 || q > 1)
                {
                    LPW(X, W, buf->q_bgc, b);
                }
            }
        }
    }

    freeMatrix(&temp);
}

// Calculates a B x G X C tensor with the probabilities of each district
/*
void getProbability(EMBuffers *buf, Matrix *V, Matrix *alpha, Matrix *beta)
{
    int B = buf->B, G = buf->G, Cminus1 = alpha->rows, C = buf->C;

    // ---- Generate needed matrices
    Matrix alphaTransposed = transposeMatrix(alpha);

    // ---- Multiply V and alpha transposed
    Matrix VxA = matrixMultiplication(V, &alphaTransposed);

    // ---- Exponentiate
    for (int b = 0; b < B; b++)
    { // --- For each district
        for (int g = 0; g < G; g++)
        { // --- For each group
            double sum = 0.0;
            for (int c = 0; c < Cminus1; c++)
            { // --- For each candidate
                // Obtain the exponential of the linear combination
                double u = MATRIX_AT_PTR(beta, g, c) + MATRIX_AT(VxA, b, c);
                double ex = exp(u);
                MATRIX_AT(buf->prob[b], g, c) = exp(u);
                sum += ex;
            }

            // Base line candidate
            MATRIX_AT(buf->prob[b], g, Cminus1) = 1;
            sum += 1;

            for (int c = 0; c < C; c++)
            { // --- For each candidate
                // Normalize
                MATRIX_AT(buf->prob[b], g, c) /= sum;
            }
        }
    }
    // Free matrices
    freeMatrix(&alphaTransposed);
    freeMatrix(&VxA);
}
*/
// Compute and normalize buf->prob[b][g][c] = softmax_c( beta[g,c] + (V\times \alpha^T)[b,c] )
// buf->C must be \alpha->rows+1
void getProbability(EMBuffers *buf,
                    Matrix *V,           // B×A
                    const Matrix *alpha, // (C-1)×A
                    const Matrix *beta)  // G×(C-1)
{
    int B = buf->B;
    int G = buf->G;
    int Cminus1 = alpha->rows;
    int C = buf->C; // = Cminus1+1

    // 1) compute V × αᵀ into the preallocated buffer
    char transA = 'N';
    char transB = 'T';
    double one = 1.0, zero = 0.0;
    BLAS_INT m = B;
    BLAS_INT n = Cminus1;
    BLAS_INT k = V->cols;
    BLAS_INT lda = V->rows;
    BLAS_INT ldb = alpha->rows;
    BLAS_INT ldc = buf->VxA.rows;

    F77_CALL(dgemm)
    (&transA, &transB, &m, &n, &k, &one, V->data, &lda, alpha->data, &ldb, &zero, buf->VxA.data, &ldc FCONE FCONE);

    // 2) exponentiate & normalize into buf->prob (stabilized softmax)
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            double max_u = 0.0; // include baseline logit = 0
            for (int c = 0; c < Cminus1; ++c)
            {
                double u = MATRIX_AT_PTR(beta, g, c) + MATRIX_AT(buf->VxA, b, c);
                if (u > max_u)
                    max_u = u;
            }

            double sum = exp(-max_u); // baseline term
            for (int c = 0; c < Cminus1; ++c)
            {
                double u = MATRIX_AT_PTR(beta, g, c) + MATRIX_AT(buf->VxA, b, c);
                double ex = exp(u - max_u);
                MATRIX_AT(buf->prob[b], g, c) = ex;
                sum += ex;
            }
            MATRIX_AT(buf->prob[b], g, Cminus1) = exp(-max_u);

            // normalize all C entries
            for (int c = 0; c < C; ++c)
            {
                MATRIX_AT(buf->prob[b], g, c) /= sum;
            }
        }
    }
}

static void compute_parametric_ballot_scores(Matrix *X, Matrix *W, EMBuffers *buf)
{
    int B = buf->B, G = buf->G, C = buf->C;
    const double eps = 1e-12;

    for (int b = 0; b < B; ++b)
    {
        double wsum = 0.0;
        for (int g = 0; g < G; ++g)
            wsum += MATRIX_AT_PTR(W, b, g);
        if (!isfinite(wsum) || wsum <= eps)
            wsum = 1.0;

        double ll_b = 0.0;
        for (int c = 0; c < C; ++c)
        {
            double s = 0.0;
            for (int g = 0; g < G; ++g)
                s += MATRIX_AT_PTR(W, b, g) * MATRIX_AT(buf->prob[b], g, c);
            MATRIX_AT(buf->S_bc, b, c) = s;

            double marginal = s / wsum;
            if (!isfinite(marginal) || marginal < eps)
                marginal = eps;
            ll_b += MATRIX_AT_PTR(X, b, c) * log(marginal);
        }
        buf->ballot_loglik[b] = ll_b;
    }
}

static void E_step_mult(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf)
{

    int B = buf->B, G = buf->G, C = buf->C;

    // ---- Get the probabilities
    getProbability(buf, V, &buf->alpha, &buf->beta);

    compute_parametric_ballot_scores(X, W, buf);

    // --- Compute q_bgc
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            double denom = 0.0;
            // --- first compute numerators into local array
            for (int c = 0; c < C; ++c)
            {
                double n = MATRIX_AT(buf->prob[b], g, c) * MATRIX_AT_PTR(X, b, c);
                double d = MATRIX_AT(buf->S_bc, b, c) - MATRIX_AT(buf->prob[b], g, c);
                if (!isfinite(d) || fabs(d) < 1e-12)
                    d = (d < 0.0) ? -1e-12 : 1e-12;
                double v = n / d;
                if (!isfinite(v) || v < 0.0)
                    v = 0.0;
                MATRIX_AT(buf->q_bgc[b], g, c) = v;
                denom += v;
            }
            // --- normalize
            bool fallback_q = false;
            if (!isfinite(denom) || denom <= 0.0)
            {
                denom = 1.0;
                fallback_q = true;
            }
            for (int c = 0; c < C; ++c)
            {
                if (fallback_q)
                    MATRIX_AT(buf->q_bgc[b], g, c) = MATRIX_AT(buf->prob[b], g, c);
                else
                    MATRIX_AT(buf->q_bgc[b], g, c) /= denom;
            }
        }
    }
}

static void E_step_mvn(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf, bool use_cdf, QMethodInput *q_params)
{
    int B = buf->B, G = buf->G, C = buf->C;

    getProbability(buf, V, &buf->alpha, &buf->beta);

    double *q_flat = (double *)Calloc((size_t)B * (size_t)G * (size_t)C, double);
    Matrix X_ctx = createMatrix(C, B);
    for (int b = 0; b < B; ++b)
        for (int c = 0; c < C; ++c)
            MATRIX_AT(X_ctx, c, b) = MATRIX_AT_PTR(X, b, c);

    EMContext ctx = {0};
    ctx.X = X_ctx;
    ctx.W = *W;
    ctx.probabilities = buf->prob[0];
    ctx.q = q_flat;
    ctx.ballot_loglik = buf->ballot_loglik;
    ctx.B = (uint32_t)B;
    ctx.G = (uint16_t)G;
    ctx.C = (uint16_t)C;

    uint32_t old_total_ballots = TOTAL_BALLOTS;
    uint16_t old_total_candidates = TOTAL_CANDIDATES;
    uint16_t old_total_groups = TOTAL_GROUPS;
    TOTAL_BALLOTS = (uint32_t)B;
    TOTAL_CANDIDATES = (uint16_t)C;
    TOTAL_GROUPS = (uint16_t)G;

    QMethodInput params = {0};
    if (q_params != NULL)
        params = *q_params;
    params.computeLL = true;
    if (use_cdf)
    {
        if (params.monteCarloIter <= 0)
            params.monteCarloIter = 1000;
        if (!(params.errorThreshold > 0.0) || !isfinite(params.errorThreshold))
            params.errorThreshold = 1e-6;
        if (params.simulationMethod == NULL)
            params.simulationMethod = "genz";
        allocateSeed(&ctx, params);
    }

    double ll = 0.0;
    if (use_cdf)
        computeQMultivariateCDFByBallot(&ctx, buf->prob, params, &ll);
    else
        computeQMultivariatePDFByBallot(&ctx, buf->prob, params, &ll);

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int c = 0; c < C; ++c)
                MATRIX_AT(buf->q_bgc[b], g, c) = Q_3D(q_flat, b, g, c, G, C);
        }
    }

    if (ctx.cdf_seeds != NULL)
        Free(ctx.cdf_seeds);
    TOTAL_BALLOTS = old_total_ballots;
    TOTAL_CANDIDATES = old_total_candidates;
    TOTAL_GROUPS = old_total_groups;
    freeMatrix(&X_ctx);
    Free(q_flat);
}

static void E_step_method(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf, const char *q_method, QMethodInput *q_params)
{
    if (q_method != NULL && strcmp(q_method, "mvn_pdf") == 0)
    {
        E_step_mvn(X, W, V, buf, false, q_params);
    }
    else if (q_method != NULL && strcmp(q_method, "mvn_cdf") == 0)
    {
        E_step_mvn(X, W, V, buf, true, q_params);
    }
    else
    {
        E_step_mult(X, W, V, buf);
    }
}

void E_step(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf)
{
    E_step_method(X, W, V, buf, "mult", NULL);
}

double objective_function(Matrix *W, Matrix *V, EMBuffers *buf, const Matrix *alpha_eval, const Matrix *beta_eval,
                          const bool prob_valid)
{

    int B = buf->B, G = buf->G, C = buf->C;

    double loss = 0.0;

    // --- Get probabilities
    if (!prob_valid)
        getProbability(buf, V, alpha_eval, beta_eval);

    // --- Get the dot product
    for (int b = 0; b < B; b++)
    { // --- For each ballot box
        for (int g = 0; g < G; g++)
        { // --- For each group
            // Must be a continuos pointer, hence, the macro can't be used
            double dot = 0.0;
            for (int c = 0; c < C; c++)
            {
                double q = MATRIX_AT(buf->q_bgc[b], g, c);
                double p = MATRIX_AT(buf->prob[b], g, c);
                dot += q * log(fmax(p, 1e-12));
            }
            loss -= MATRIX_AT_PTR(W, b, g) * dot; // Check if it is to sum or to substract
        }
    }
    return loss;
}

void compute_gradients(const Matrix *W, Matrix *V, EMBuffers *buf, const Matrix *alpha_eval, const Matrix *beta_eval,
                       const bool prob_valid)
{

    int B = buf->B, G = buf->G, Cminus1 = buf->C - 1, A = buf->A;

    // --- Get probabilities
    if (!prob_valid)
        getProbability(buf, V, alpha_eval, beta_eval);

    for (int g = 0; g < G; g++)
    {
        for (int c = 0; c < Cminus1; c++)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int b = 0; b < B; b++)
            {
                double w = MATRIX_AT_PTR(W, b, g);
                sum1 += w * MATRIX_AT(buf->q_bgc[b], g, c);
                sum2 += w * MATRIX_AT(buf->prob[b], g, c);
            }
            MATRIX_AT(buf->grad_beta, g, c) = sum1 - sum2;
        }
    }

    for (int c = 0; c < Cminus1; c++)
    {
        for (int a = 0; a < A; a++)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int b = 0; b < B; b++)
            {
                for (int g = 0; g < G; g++)
                {
                    double w = MATRIX_AT_PTR(W, b, g);
                    double q = MATRIX_AT(buf->q_bgc[b], g, c);
                    double v = MATRIX_AT_PTR(V, b, a);
                    double p = MATRIX_AT(buf->prob[b], g, c);
                    sum1 += w * q * v;
                    sum2 += w * p * v;
                }
            }
            MATRIX_AT(buf->grad_alpha, c, a) = sum1 - sum2;
        }
    }
}

void compute_hessian(const Matrix *W, // B \times G
                     Matrix *V,       // B \times A
                     EMBuffers *buf,
                     const Matrix *alpha_eval, // (C-1) \times A
                     const Matrix *beta_eval,  // G \times (C-1)
                     const bool prob_valid)
{

    int B = buf->B, G = buf->G, Cm = buf->C - 1, A = buf->A, D = buf->D;
    int d_alpha = Cm * A;
    int d_beta = G * Cm;
    /* 1) zero out the entire Hessian */
    memset(buf->H.data, 0, D * D * sizeof(double));

    /* 2) compute p_bgc once */
    if (!prob_valid)
        getProbability(buf, V, alpha_eval, beta_eval);

    size_t n_iters = (size_t)B * G * Cm * Cm;
    // #ifdef _OPENMP
    // #pragma omp parallel for collapse(4) if (n_iters > 500) schedule(static)
    // #endif
    for (int b = 0; b < B; b++)
    {
        for (int g = 0; g < G; g++)
        {
            double w = MATRIX_AT_PTR(W, b, g);

            for (int ci = 0; ci < Cm; ci++)
            {
                double p_ci = MATRIX_AT(buf->prob[b], g, ci);

                for (int cj = 0; cj < Cm; cj++)
                {
                    double p_cj = MATRIX_AT(buf->prob[b], g, cj);
                    // negative of second‐derivative of log‑prob
                    double base = p_ci * p_cj + (ci == cj ? -p_ci : 0.0);
                    double factor = -w * base;

                    // alpha–alpha block
                    for (int ai = 0; ai < A; ai++)
                    {
                        double v_ai = MATRIX_AT_PTR(V, b, ai);
                        int row_alpha = ai * Cm + ci;
                        for (int aj = 0; aj < A; aj++)
                        {
                            double v_aj = MATRIX_AT_PTR(V, b, aj);
                            int col_alpha = aj * Cm + cj;
                            MATRIX_AT(buf->H, row_alpha, col_alpha) += v_ai * v_aj * factor;
                        }
                    }

                    // alpha–beta and beta–alpha blocks
                    for (int ai = 0; ai < A; ai++)
                    {
                        double v_ai = MATRIX_AT_PTR(V, b, ai);
                        int row_alpha = ai * Cm + ci;
                        int col_beta = d_alpha + g * Cm + cj;
                        double upd = v_ai * factor;
                        MATRIX_AT(buf->H, row_alpha, col_beta) += upd;
                        MATRIX_AT(buf->H, col_beta, row_alpha) += upd;
                    }

                    // beta–beta block
                    {
                        int row_beta = d_alpha + g * Cm + ci;
                        int col_beta = d_alpha + g * Cm + cj;
                        MATRIX_AT(buf->H, row_beta, col_beta) += factor;
                    }
                }
            }
        }
    }
}

// ----- HELPER FUNCTION ----- //
// Packs grad_alpha (C–1 x A, column‐major) followed by grad_beta (G x C–1, row‐major)
// into the flat vector g[0..D-1] to match the Hessian indexing.
static void pack_gradients(EMBuffers *buf)
{
    int Cminus1 = buf->grad_alpha.rows;
    int A = buf->grad_alpha.cols;
    int G = buf->grad_beta.rows;
    int idx = 0;

    // \alpha block
    for (int a = 0; a < A; a++)
    { // --- For each attribute
        for (int c = 0; c < Cminus1; c++)
        { // --- For each candidate
            buf->gvec[idx++] = MATRIX_AT(buf->grad_alpha, c, a);
        }
    }

    // \beta block (g-major / row-major)
    for (int gi = 0; gi < G; gi++)
    { // --- For each group
        for (int c = 0; c < Cminus1; c++)
        { // --- For each candidate
            buf->gvec[idx++] = MATRIX_AT(buf->grad_beta, gi, c);
        }
    }
}

// Unpacks a flat D‐vector v[] back into two matrices, in the same order.
static void unpack_step(EMBuffers *buf)
{
    int Cminus1 = buf->grad_alpha.rows;
    int A = buf->grad_alpha.cols;
    int G = buf->grad_beta.rows;
    int idx = 0;

    // \alpha block
    for (int a = 0; a < A; ++a)
    { // --- For each attribute
        for (int c = 0; c < Cminus1; ++c)
        { // --- For each candidate
            MATRIX_AT(buf->grad_alpha, c, a) = buf->vvec[idx++];
        }
    }

    // \beta block (g-major / row-major)
    for (int gi = 0; gi < G; ++gi)
    { // --- For each group
        for (int c = 0; c < Cminus1; ++c)
        { // --- For each candidate
            MATRIX_AT(buf->grad_beta, gi, c) = buf->vvec[idx++];
        }
    }
}

int Newton_damped(Matrix *W, // B \times G weights
                  Matrix *V, // B \times A covariates
                  EMBuffers *buf, double tol, int max_iter, double alpha_bs, double beta_bs)
{
    int B = V->rows;
    int A = V->cols;
    int Cminus1 = buf->alpha.rows;
    int C = Cminus1 + 1;
    int G = buf->beta.rows;
    int D = Cminus1 * A + G * Cminus1;

    // ---- Clone initial parameters
    Matrix *alpha = copMatrixPtr(&buf->alpha);
    Matrix *beta = copMatrixPtr(&buf->beta);

    // ---- Allocate temporaries
    // double *gvec = (double *)Calloc(D, double);
    // double *vvec = (double *)Calloc(D, double);
    // Matrix dalpha = createMatrix(Cminus1, A);
    // Matrix dbeta = createMatrix(G, Cminus1);
    // Matrix H = createMatrix(D, D);

    int iter;
    for (iter = 0; iter < max_iter; iter++)
    {
        // ---- Evaluate loss, gradient, Hessian at the current points
        // Loss
        getProbability(buf, V, alpha, beta);
        double f0 = objective_function(W, V, buf, alpha, beta, true);

        // Gradients
        compute_gradients(W, V, buf, alpha, beta, true);
        pack_gradients(buf);

        // Hessian (prezero H)
        compute_hessian(W, V, buf, alpha, beta, true);

        // Hessian damping: H <- (1-eta) H + eta * max_diag * I
        const double eta = tol;
        double max_diag = MATRIX_AT(buf->H, 0, 0);
        for (int i = 1; i < D; ++i)
        {
            max_diag = fmax(max_diag, MATRIX_AT(buf->H, i, i));
        }
        for (int i = 0; i < D; ++i)
        {
            for (int j = 0; j < D; ++j)
            {
                double hij = MATRIX_AT(buf->H, i, j);
                if (i == j)
                {
                    MATRIX_AT(buf->H, i, j) = (1.0 - eta) * hij + eta * max_diag;
                }
                else
                {
                    MATRIX_AT(buf->H, i, j) = (1.0 - eta) * hij;
                }
            }
        }

        // Solve H v = -g, for approximating it with Taylor expansion
        solve_linear_system(D, buf->H.data, buf->gvec, buf->vvec);
        for (int i = 0; i < D; i++)
        {
            buf->vvec[i] = -buf->vvec[i];
        }
        unpack_step(buf);

        // Convergence check
        double g_inf = 0;
        for (int i = 0; i < D; i++)
        {
            g_inf = fmax(g_inf, fabs(buf->gvec[i]));
        }
        if (g_inf < tol)
            break;

        // Armijo backtracking line search
        Matrix *alpha_t = copMatrixPtr(alpha);
        Matrix *beta_t = copMatrixPtr(beta);
        size_t alpha_elems = alpha->rows * alpha->cols;
        size_t beta_elems = beta->rows * beta->cols;
        BLAS_INT inc = 1;
        BLAS_INT alpha_n = (BLAS_INT)alpha_elems;
        BLAS_INT beta_n = (BLAS_INT)beta_elems;
        double t = 1.0; // We start with t = 1
        const double t_min = 1e-10;
        // compute grad(loss) * v (gvec stores -grad(loss))
        double gv = 0;
        for (int i = 0; i < D; i++)
            gv += -buf->gvec[i] * buf->vvec[i];

        while (1)
        {
            // Trial parameters
            F77_CALL(dcopy)(&alpha_n, alpha->data, &inc, alpha_t->data, &inc);
            F77_CALL(daxpy)(&alpha_n, &t, buf->grad_alpha.data, &inc, alpha_t->data, &inc);
            F77_CALL(dcopy)(&beta_n, beta->data, &inc, beta_t->data, &inc);
            F77_CALL(daxpy)(&beta_n, &t, buf->grad_beta.data, &inc, beta_t->data, &inc);

            double f_trial = objective_function(W, V, buf, alpha_t, beta_t, false);
            if (f_trial <= f0 + alpha_bs * t * gv)
            {
                break;
            }
            t *= beta_bs;
            if (t < t_min)
                break;
        }
        if (t < t_min)
            t = t_min;
        freeMatrix(alpha_t);
        freeMatrix(beta_t);

        // Update parameters. i.e, alpha and beta
        for (int i = 0; i < Cminus1; i++)
            for (int j = 0; j < A; j++)
                MATRIX_AT_PTR(alpha, i, j) += t * MATRIX_AT(buf->grad_alpha, i, j);
        for (int i = 0; i < G; i++)
            for (int j = 0; j < Cminus1; j++)
                MATRIX_AT_PTR(beta, i, j) += t * MATRIX_AT(buf->grad_beta, i, j);

    } // --- Newton iteration finishes

    // Copy results out
    size_t alpha_elems = alpha->rows * alpha->cols;
    size_t alpha_bytes = alpha_elems * sizeof(double);
    memcpy(buf->alpha.data, alpha->data, alpha_bytes);

    size_t beta_elems = beta->rows * beta->cols;
    size_t beta_bytes = beta_elems * sizeof(double);
    memcpy(buf->beta.data, beta->data, beta_bytes);

    // Cleanup
    freeMatrix(alpha);
    freeMatrix(beta);

    return iter + 1;
}

static double compute_ll_multinomial_log(const Matrix *X, // BxC
                                         const Matrix *W, // BxG
                                         Matrix *V,       // BxA
                                         EMBuffers *buf)
{
    int B = X->rows;
    int C = X->cols;
    int A = V->cols;
    int G = W->cols;

    double total_ll = 0.0;
    const double log_eps = log(1e-12);

    getProbability(buf, V, &buf->alpha, &buf->beta);

    for (int b = 0; b < B; ++b)
    {
        // compute denominator = sum_g w_bg[b,g]
        double wsum = 0;
        for (int g = 0; g < G; ++g)
            wsum += MATRIX_AT_PTR(W, b, g);
        double denom = wsum + 1e-12;
        double log_denom = log(denom);

        // multinomial factorial term
        double xb = 0;
        for (int c = 0; c < C; ++c)
        {
            double x = MATRIX_AT_PTR(X, b, c);
            double lg = lgamma(x + 1.0);
            total_ll -= lg;
            xb += x;
        }
        total_ll += lgamma(xb + 1.0);

        // data term \sum x_bc \cdot log(p_bc), computed in log-space
        for (int c = 0; c < C; ++c)
        {
            double max_log = -INFINITY;
            bool any = false;
            for (int g = 0; g < G; ++g)
            {
                double w = MATRIX_AT_PTR(W, b, g);
                double p = MATRIX_AT(buf->prob[b], g, c);
                if (w > 0.0 && p > 0.0)
                {
                    double log_term = log(w) + log(p);
                    if (!any || log_term > max_log)
                        max_log = log_term;
                    any = true;
                }
            }
            double log_marg = -INFINITY;
            if (any)
            {
                double sum_exp = 0.0;
                for (int g = 0; g < G; ++g)
                {
                    double w = MATRIX_AT_PTR(W, b, g);
                    double p = MATRIX_AT(buf->prob[b], g, c);
                    if (w > 0.0 && p > 0.0)
                    {
                        double log_term = log(w) + log(p);
                        sum_exp += exp(log_term - max_log);
                    }
                }
                log_marg = max_log + log(sum_exp);
            }
            double log_pbc = log_marg - log_denom;
            if (!isfinite(log_pbc) || log_pbc < log_eps)
                log_pbc = log_eps;
            total_ll += MATRIX_AT_PTR(X, b, c) * log_pbc;
        }
    }

    return total_ll;
}

void M_step(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf, const double tol, const int maxnewton, const bool verbose)
{
    int newton_iterations = Newton_damped(W, V, buf, tol, maxnewton, 0.5, 0.5);

    // if (verbose)
    // {
    //     Rprintf("The newton algorithm was made in %d iterations\n", newton_iterations - 1);
    // }
}

static void average_expected_outcomes_update_q(const Matrix *X, const Matrix *W, EMBuffers *forward, EMBuffers *reverse)
{
    const int B = forward->B;
    const int G = forward->G;
    const int C = forward->C;
    const size_t z_size = (size_t)B * (size_t)G * (size_t)C;
    double *z_avg = (double *)Calloc(z_size, double);

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            const double w_bg = MATRIX_AT_PTR(W, b, g);
            for (int c = 0; c < C; ++c)
            {
                const double x_bc = MATRIX_AT_PTR(X, b, c);
                const double z_forward = w_bg * MATRIX_AT(forward->q_bgc[b], g, c);
                const double z_reverse = x_bc * MATRIX_AT(reverse->q_bgc[b], c, g);
                Q_3D(z_avg, b, g, c, G, C) = 0.5 * (z_forward + z_reverse);
            }
        }
    }

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            const double den = MATRIX_AT_PTR(W, b, g);
            double sum = 0.0;
            for (int c = 0; c < C; ++c)
            {
                double value = den > 0.0 ? Q_3D(z_avg, b, g, c, G, C) / den : 1.0 / (double)C;
                if (!isfinite(value) || value < 0.0)
                    value = 0.0;
                MATRIX_AT(forward->q_bgc[b], g, c) = value;
                sum += value;
            }
            if (!isfinite(sum) || sum <= 0.0)
            {
                const double uniform = 1.0 / (double)C;
                for (int c = 0; c < C; ++c)
                    MATRIX_AT(forward->q_bgc[b], g, c) = uniform;
            }
            else
            {
                for (int c = 0; c < C; ++c)
                    MATRIX_AT(forward->q_bgc[b], g, c) /= sum;
            }
        }
    }

    for (int b = 0; b < B; ++b)
    {
        for (int c = 0; c < C; ++c)
        {
            const double den = MATRIX_AT_PTR(X, b, c);
            double sum = 0.0;
            for (int g = 0; g < G; ++g)
            {
                double value = den > 0.0 ? Q_3D(z_avg, b, g, c, G, C) / den : 1.0 / (double)G;
                if (!isfinite(value) || value < 0.0)
                    value = 0.0;
                MATRIX_AT(reverse->q_bgc[b], c, g) = value;
                sum += value;
            }
            if (!isfinite(sum) || sum <= 0.0)
            {
                const double uniform = 1.0 / (double)G;
                for (int g = 0; g < G; ++g)
                    MATRIX_AT(reverse->q_bgc[b], c, g) = uniform;
            }
            else
            {
                for (int g = 0; g < G; ++g)
                    MATRIX_AT(reverse->q_bgc[b], c, g) /= sum;
            }
        }
    }

    Free(z_avg);
}

static double *flatten_q_bgc(const EMBuffers *buf)
{
    const size_t ballot_size = (size_t)buf->G * (size_t)buf->C;
    double *q_flat = (double *)Calloc((size_t)buf->B * ballot_size, double);
    for (int b = 0; b < buf->B; ++b)
        memcpy(q_flat + (size_t)b * ballot_size, buf->q_bgc[b].data, ballot_size * sizeof(double));
    return q_flat;
}

static void copy_q_ballot_from_flat(EMBuffers *buf, const double *q_flat, int b)
{
    const size_t ballot_size = (size_t)buf->G * (size_t)buf->C;
    memcpy(buf->q_bgc[b].data, q_flat + (size_t)b * ballot_size, ballot_size * sizeof(double));
}

static Matrix transpose_ballots_candidates(const Matrix *X)
{
    Matrix Xt = createMatrix(X->cols, X->rows);
    for (int b = 0; b < X->rows; ++b)
        for (int c = 0; c < X->cols; ++c)
            MATRIX_AT(Xt, c, b) = MATRIX_AT_PTR(X, b, c);
    return Xt;
}

static void apply_joint_or_separate_lp(Matrix *X, Matrix *W, EMBuffers *forward, EMBuffers *reverse)
{
    Matrix X_ctx = transpose_ballots_candidates(X);
    double *q_forward = flatten_q_bgc(forward);
    double *q_reverse = flatten_q_bgc(reverse);

    EMContext ctx_forward = {0};
    ctx_forward.X = X_ctx;
    ctx_forward.W = *W;
    ctx_forward.q = q_forward;
    ctx_forward.B = (uint32_t)forward->B;
    ctx_forward.G = (uint16_t)forward->G;
    ctx_forward.C = (uint16_t)forward->C;

    EMContext ctx_reverse = {0};
    ctx_reverse.q = q_reverse;
    ctx_reverse.B = (uint32_t)reverse->B;
    ctx_reverse.G = (uint16_t)reverse->G;
    ctx_reverse.C = (uint16_t)reverse->C;

    for (int b = 0; b < forward->B; ++b)
    {
        int status = LPW_joint_symmetric_ctx(&ctx_forward, &ctx_reverse, b);
        if (status == 0)
        {
            copy_q_ballot_from_flat(forward, q_forward, b);
            copy_q_ballot_from_flat(reverse, q_reverse, b);
        }
        else
        {
            LPW(X, W, forward->q_bgc, b);
            LPW(W, X, reverse->q_bgc, b);
        }
    }

    freeMatrix(&X_ctx);
    Free(q_forward);
    Free(q_reverse);
}

Matrix *EM_Algorithm_Method(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                            const double maxtime, const double ll_threshold, const int maxnewton, const bool verbose,
                            double *out_elapsed, int *total_iterations, double *logLikelihood,
                            Matrix **out_q, Matrix **out_expected,
                            const char *q_method, QMethodInput *q_params,
                            const char *adjust_prob_cond_method, bool adjust_prob_cond_every)
{
    int B = V->rows;
    int A = V->cols;
    int Cm = alpha->rows;
    int G = beta->rows;

    // Initialize buffers
    EMBuffers buf;
    init_EMBuffers(&buf, B, G, Cm, A);
    buf.alpha = copMatrix(alpha);
    buf.beta = copMatrix(beta);

    bool use_project_lp = (adjust_prob_cond_method != NULL && strcmp(adjust_prob_cond_method, "project_lp") == 0);
    bool use_lp = (adjust_prob_cond_method != NULL && strcmp(adjust_prob_cond_method, "lp") == 0);
    double *scale_factors = NULL;
    Matrix norm = (Matrix){0};
    if (use_project_lp)
    {
        scale_factors = (double *)Calloc(B, double);
        for (int b = 0; b < B; b++)
            scale_factors[b] = 1.0;
        if (hasMismatch(X, W))
            precomputeScaleFactors(scale_factors, X, W);
        norm = precomputeNorm(scale_factors, W);
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t0); // Start timer
    double current_ll = -DBL_MAX;            // log-likelihood
    double new_ll = -DBL_MAX;                // log-likelihood
    double tol = 1.0;
    for (int iter = 0; iter < maxiter; iter++)
    {
        *total_iterations += 1;
        tol = 1.0 / (iter + 1);
        E_step_method(X, W, V, &buf, q_method, q_params);
        if (adjust_prob_cond_every)
        {
            if (use_project_lp)
                projectQ(X, W, &buf, &norm, scale_factors);
            else if (use_lp)
                for (int b = 0; b < B; b++)
                    LPW(X, W, buf.q_bgc, b);
        }
        M_step(X, W, V, &buf, tol, maxnewton, verbose);
        new_ll = compute_ll_multinomial_log(X, W, V, &buf);

        // Check if the user want to interrupt the process
        if (iter % 5 == 0)
            R_CheckUserInterrupt();

        if (verbose)
        {
            Rprintf("Iteration %d: log-likelihood = %.4f\n", iter + 1, new_ll);
        }

        // Check for convergence
        if (current_ll >= new_ll && verbose)
            Rprintf("Log-likelihood did not increase: %.6f -> %.6f\n", current_ll, new_ll);

        if (fabs(new_ll - current_ll) <= ll_threshold || current_ll >= new_ll)
        {
            if (verbose)
            {
                Rprintf("Converged after %d iterations.\n", iter + 1);
            }
            break;
        }
        current_ll = new_ll;
    }
    E_step_method(X, W, V, &buf, q_method, q_params);
    if (use_project_lp)
    {
        projectQ(X, W, &buf, &norm, scale_factors);
        M_step(X, W, V, &buf, tol, maxnewton, verbose);
        new_ll = compute_ll_multinomial_log(X, W, V, &buf);
    }
    else if (use_lp)
    {
        for (int b = 0; b < B; b++)
            LPW(X, W, buf.q_bgc, b);
        M_step(X, W, V, &buf, tol, maxnewton, verbose);
        new_ll = compute_ll_multinomial_log(X, W, V, &buf);
    }

    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

    if (out_q != NULL)
    {
        Matrix *q_out = alloc_matrix_array(B, G, buf.C);
        copy_matrix_array(q_out, buf.q_bgc, B);
        *out_q = q_out;
    }
    if (out_expected != NULL)
    {
        Matrix *expected_out = alloc_matrix_array(B, G, buf.C);
        compute_expected_outcome(W, buf.q_bgc, expected_out, B, G, buf.C);
        *out_expected = expected_out;
    }

    // Compute elapsed seconds
    double sec = (double)(t1.tv_sec - t0.tv_sec);
    double nsec = (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
    *out_elapsed = sec + nsec;
    *logLikelihood = new_ll;

    size_t na = buf.alpha.rows * buf.alpha.cols;
    memcpy(alpha->data, buf.alpha.data, na * sizeof(double));
    size_t nb = buf.beta.rows * buf.beta.cols;
    memcpy(beta->data, buf.beta.data, nb * sizeof(double));
    Matrix *finalProb = buf.prob;
    // detach buf.prob so we don't free it:
    buf.prob = NULL;
    free_EMBuffers(&buf);
    if (scale_factors != NULL)
        Free(scale_factors);
    if (norm.data != NULL)
        freeMatrix(&norm);

    // Matrix *finalProbability = getProbability(V, beta, alpha);
    return finalProb; // Return the final probabilities
}

Matrix *EM_Algorithm_Symmetric_Method(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                                      const double maxtime, const double ll_threshold, const int maxnewton,
                                      const bool verbose, double *out_elapsed, int *total_iterations,
                                      double *logLikelihood, Matrix **out_q, Matrix **out_expected,
                                      const char *q_method, QMethodInput *q_params,
                                      const char *adjust_prob_cond_method, bool adjust_prob_cond_every)
{
    int B = V->rows;
    int A = V->cols;
    int C = X->cols;
    int G = W->cols;
    int Cm = C - 1;
    int Gm = G - 1;

    if (C < 2 || G < 2)
        error("Parametric symmetric EM requires at least two candidates and two groups.");

    EMBuffers forward;
    EMBuffers reverse;
    init_EMBuffers(&forward, B, G, Cm, A);
    init_EMBuffers(&reverse, B, C, Gm, A);
    forward.alpha = copMatrix(alpha);
    forward.beta = copMatrix(beta);
    reverse.alpha = createMatrix(Gm, A);
    reverse.beta = createMatrix(C, Gm);

    bool use_project_lp = (adjust_prob_cond_method != NULL && strcmp(adjust_prob_cond_method, "project_lp") == 0);
    bool use_lp = (adjust_prob_cond_method != NULL && strcmp(adjust_prob_cond_method, "lp") == 0);

    double *scale_forward = NULL;
    double *scale_reverse = NULL;
    Matrix norm_forward = (Matrix){0};
    Matrix norm_reverse = (Matrix){0};
    if (use_project_lp)
    {
        scale_forward = (double *)Calloc(B, double);
        scale_reverse = (double *)Calloc(B, double);
        for (int b = 0; b < B; ++b)
        {
            scale_forward[b] = 1.0;
            scale_reverse[b] = 1.0;
        }
        if (hasMismatch(X, W))
            precomputeScaleFactors(scale_forward, X, W);
        if (hasMismatch(W, X))
            precomputeScaleFactors(scale_reverse, W, X);
        norm_forward = precomputeNorm(scale_forward, W);
        norm_reverse = precomputeNorm(scale_reverse, X);
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t0);
    double old_ll_forward = -DBL_MAX;
    double old_ll_reverse = -DBL_MAX;
    double new_ll_forward = -DBL_MAX;
    double new_ll_reverse = -DBL_MAX;
    double tol = 1.0;

    for (int iter = 0; iter < maxiter; ++iter)
    {
        *total_iterations += 1;
        tol = 1.0 / (iter + 1);

        E_step_method(X, W, V, &forward, q_method, q_params);
        if (adjust_prob_cond_every && use_project_lp)
            projectQ(X, W, &forward, &norm_forward, scale_forward);

        E_step_method(W, X, V, &reverse, q_method, q_params);
        if (adjust_prob_cond_every && use_project_lp)
            projectQ(W, X, &reverse, &norm_reverse, scale_reverse);

        if (adjust_prob_cond_every && use_lp)
            apply_joint_or_separate_lp(X, W, &forward, &reverse);

        average_expected_outcomes_update_q(X, W, &forward, &reverse);

        M_step(X, W, V, &forward, tol, maxnewton, verbose);
        M_step(W, X, V, &reverse, tol, maxnewton, verbose);

        new_ll_forward = compute_ll_multinomial_log(X, W, V, &forward);
        new_ll_reverse = compute_ll_multinomial_log(W, X, V, &reverse);

        if (iter % 5 == 0)
            R_CheckUserInterrupt();

        if (verbose)
        {
            Rprintf("Iteration %d: forward log-likelihood = %.4f, reverse log-likelihood = %.4f\n", iter + 1,
                    new_ll_forward, new_ll_reverse);
        }

        bool converged_forward =
            fabs(new_ll_forward - old_ll_forward) <= ll_threshold || old_ll_forward >= new_ll_forward;
        bool converged_reverse =
            fabs(new_ll_reverse - old_ll_reverse) <= ll_threshold || old_ll_reverse >= new_ll_reverse;
        if (converged_forward && converged_reverse)
        {
            if (verbose)
                Rprintf("Parametric symmetric EM converged after %d iterations.\n", iter + 1);
            break;
        }

        old_ll_forward = new_ll_forward;
        old_ll_reverse = new_ll_reverse;
    }

    E_step_method(X, W, V, &forward, q_method, q_params);
    E_step_method(W, X, V, &reverse, q_method, q_params);
    average_expected_outcomes_update_q(X, W, &forward, &reverse);
    if (use_project_lp)
    {
        projectQ(X, W, &forward, &norm_forward, scale_forward);
        projectQ(W, X, &reverse, &norm_reverse, scale_reverse);
    }
    else if (use_lp)
    {
        apply_joint_or_separate_lp(X, W, &forward, &reverse);
    }

    M_step(X, W, V, &forward, tol, maxnewton, verbose);
    M_step(W, X, V, &reverse, tol, maxnewton, verbose);
    new_ll_forward = compute_ll_multinomial_log(X, W, V, &forward);
    new_ll_reverse = compute_ll_multinomial_log(W, X, V, &reverse);

    clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

    if (out_q != NULL)
    {
        Matrix *q_out = alloc_matrix_array(B, G, forward.C);
        copy_matrix_array(q_out, forward.q_bgc, B);
        *out_q = q_out;
    }
    if (out_expected != NULL)
    {
        Matrix *expected_out = alloc_matrix_array(B, G, forward.C);
        compute_expected_outcome(W, forward.q_bgc, expected_out, B, G, forward.C);
        *out_expected = expected_out;
    }

    double sec = (double)(t1.tv_sec - t0.tv_sec);
    double nsec = (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
    *out_elapsed = sec + nsec;
    *logLikelihood = 0.5 * (new_ll_forward + new_ll_reverse);

    size_t na = forward.alpha.rows * forward.alpha.cols;
    memcpy(alpha->data, forward.alpha.data, na * sizeof(double));
    size_t nb = forward.beta.rows * forward.beta.cols;
    memcpy(beta->data, forward.beta.data, nb * sizeof(double));

    Matrix *finalProb = forward.prob;
    forward.prob = NULL;

    freeMatrix(&forward.alpha);
    freeMatrix(&forward.beta);
    freeMatrix(&reverse.alpha);
    freeMatrix(&reverse.beta);
    free_EMBuffers(&forward);
    free_EMBuffers(&reverse);
    if (scale_forward != NULL)
        Free(scale_forward);
    if (scale_reverse != NULL)
        Free(scale_reverse);
    if (norm_forward.data != NULL)
        freeMatrix(&norm_forward);
    if (norm_reverse.data != NULL)
        freeMatrix(&norm_reverse);

    return finalProb;
}

Matrix *EM_Algorithm(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                     const double maxtime, const double ll_threshold, const int maxnewton, const bool verbose,
                     double *out_elapsed, int *total_iterations, double *logLikelihood,
                     Matrix **out_q, Matrix **out_expected,
                     const char *adjust_prob_cond_method, bool adjust_prob_cond_every)
{
    return EM_Algorithm_Method(X, W, V, beta, alpha, maxiter, maxtime, ll_threshold, maxnewton, verbose, out_elapsed,
                               total_iterations, logLikelihood, out_q, out_expected, "mult", NULL,
                               adjust_prob_cond_method, adjust_prob_cond_every);
}

static double update_matrix_parametric_responsibilities(const Matrix *score, const Matrix *phi,
                                                        Matrix *responsibilities, double *row_scores)
{
    const int B = score->rows;
    const int K = score->cols;
    const double eps = 1e-12;
    double ll = 0.0;

    for (int b = 0; b < B; ++b)
    {
        for (int k = 0; k < K; ++k)
            row_scores[k] = MATRIX_AT_PTR(score, b, k);

        double maxv = -INFINITY;
        for (int k = 0; k < K; ++k)
        {
            double ph = MATRIX_AT_PTR(phi, b, k);
            ph = (isfinite(ph) && ph > eps) ? ph : eps;
            double v = row_scores[k] + log(ph);
            if (isfinite(v) && v > maxv)
                maxv = v;
        }
        if (!isfinite(maxv))
            maxv = log(eps);

        double sum_exp = 0.0;
        for (int k = 0; k < K; ++k)
        {
            double ph = MATRIX_AT_PTR(phi, b, k);
            ph = (isfinite(ph) && ph > eps) ? ph : eps;
            double v = row_scores[k] + log(ph);
            if (isfinite(v))
                sum_exp += exp(v - maxv);
        }
        double ll_b = (isfinite(sum_exp) && sum_exp > 0.0) ? maxv + log(sum_exp) : log(eps);
        ll += isfinite(ll_b) ? ll_b : log(eps);

        double den = 0.0;
        for (int k = 0; k < K; ++k)
        {
            double ph = MATRIX_AT_PTR(phi, b, k);
            ph = (isfinite(ph) && ph > eps) ? ph : eps;
            double val = isfinite(row_scores[k]) ? exp(row_scores[k] + log(ph) - ll_b) : 0.0;
            MATRIX_AT_PTR(responsibilities, b, k) = val;
            den += val;
        }

        if (!isfinite(den) || den <= 0.0)
        {
            double uniform = 1.0 / (double)K;
            for (int k = 0; k < K; ++k)
                MATRIX_AT_PTR(responsibilities, b, k) = uniform;
        }
        else
        {
            for (int k = 0; k < K; ++k)
                MATRIX_AT_PTR(responsibilities, b, k) /= den;
        }
    }

    return ll;
}

static double matrix_max_delta(const Matrix *a, const Matrix *b)
{
    double delta = 0.0;
    size_t n = (size_t)a->rows * (size_t)a->cols;
    for (size_t i = 0; i < n; ++i)
    {
        double d = fabs(a->data[i] - b->data[i]);
        if (isfinite(d) && d > delta)
            delta = d;
    }
    return delta;
}

static void normalize_matrix_mixture_component(Matrix *P)
{
    const double eps = 1e-12;
    for (int g = 0; g < P->rows; ++g)
    {
        double rowsum = 0.0;
        for (int c = 0; c < P->cols; ++c)
        {
            double v = MATRIX_AT_PTR(P, g, c);
            if (!isfinite(v) || v < eps)
                v = eps;
            MATRIX_AT_PTR(P, g, c) = v;
            rowsum += v;
        }
        if (!isfinite(rowsum) || rowsum <= 0.0)
            rowsum = 1.0;
        for (int c = 0; c < P->cols; ++c)
            MATRIX_AT_PTR(P, g, c) /= rowsum;
    }
}

static void compute_membership_phi(const Matrix *V, const Matrix *beta, Matrix *phi)
{
    const int B = V->rows;
    const int A = V->cols;
    const int K = phi->cols;

    for (int b = 0; b < B; ++b)
    {
        double max_eta = 0.0;
        for (int k = 0; k < K - 1; ++k)
        {
            double eta = 0.0;
            for (int a = 0; a < A; ++a)
                eta += MATRIX_AT_PTR(V, b, a) * MATRIX_AT_PTR(beta, a, k);
            MATRIX_AT_PTR(phi, b, k) = eta;
            if (eta > max_eta)
                max_eta = eta;
        }

        double denom = exp(-max_eta);
        for (int k = 0; k < K - 1; ++k)
            denom += exp(MATRIX_AT_PTR(phi, b, k) - max_eta);

        for (int k = 0; k < K - 1; ++k)
            MATRIX_AT_PTR(phi, b, k) = exp(MATRIX_AT_PTR(phi, b, k) - max_eta) / denom;
        MATRIX_AT_PTR(phi, b, K - 1) = exp(-max_eta) / denom;
    }
}

static double membership_negative_loglik(const Matrix *V, const Matrix *responsibilities, const Matrix *beta)
{
    const int B = V->rows;
    const int A = V->cols;
    const int K = responsibilities->cols;
    const double eps = 1e-12;
    double nll = 0.0;

    Matrix phi = createMatrix(B, K);
    compute_membership_phi(V, beta, &phi);
    for (int b = 0; b < B; ++b)
    {
        for (int k = 0; k < K; ++k)
        {
            double r = MATRIX_AT_PTR(responsibilities, b, k);
            double ph = MATRIX_AT(phi, b, k);
            ph = (isfinite(ph) && ph > eps) ? ph : eps;
            nll -= r * log(ph);
        }
    }
    freeMatrix(&phi);
    return nll;
}

static void newton_membership_beta(const Matrix *V, const Matrix *responsibilities, Matrix *beta, int max_iter,
                                   double ridge)
{
    const int B = V->rows;
    const int A = V->cols;
    const int K = responsibilities->cols;
    const int H = K - 1;
    const int D = A * H;
    if (D <= 0)
        return;

    Matrix phi = createMatrix(B, K);
    double *grad = (double *)Calloc((size_t)D, double);
    double *hess = (double *)Calloc((size_t)D * (size_t)D, double);
    double *step = (double *)Calloc((size_t)D, double);
    Matrix trial = createMatrix(A, H);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        compute_membership_phi(V, beta, &phi);
        memset(grad, 0, (size_t)D * sizeof(double));
        memset(hess, 0, (size_t)D * (size_t)D * sizeof(double));

        for (int b = 0; b < B; ++b)
        {
            for (int h = 0; h < H; ++h)
            {
                double diff = MATRIX_AT(phi, b, h) - MATRIX_AT_PTR(responsibilities, b, h);
                for (int a = 0; a < A; ++a)
                {
                    int idx = h * A + a;
                    double va = MATRIX_AT_PTR(V, b, a);
                    grad[idx] += diff * va;
                }
            }

            for (int h1 = 0; h1 < H; ++h1)
            {
                double p1 = MATRIX_AT(phi, b, h1);
                for (int h2 = 0; h2 < H; ++h2)
                {
                    double p2 = MATRIX_AT(phi, b, h2);
                    double factor = p1 * ((h1 == h2 ? 1.0 : 0.0) - p2);
                    for (int a1 = 0; a1 < A; ++a1)
                    {
                        int row = h1 * A + a1;
                        double v1 = MATRIX_AT_PTR(V, b, a1);
                        for (int a2 = 0; a2 < A; ++a2)
                        {
                            int col = h2 * A + a2;
                            double v2 = MATRIX_AT_PTR(V, b, a2);
                            hess[(size_t)col * (size_t)D + (size_t)row] += factor * v1 * v2;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < D; ++i)
            hess[(size_t)i * (size_t)D + (size_t)i] += ridge;

        solve_linear_system(D, hess, grad, step);
        double gdot = 0.0;
        double step_inf = 0.0;
        for (int i = 0; i < D; ++i)
        {
            gdot += grad[i] * step[i];
            step_inf = fmax(step_inf, fabs(step[i]));
        }
        if (step_inf < 1e-8)
            break;

        double f0 = membership_negative_loglik(V, responsibilities, beta);
        double t = 1.0;
        const double c1 = 1e-4;
        bool accepted = false;
        while (t >= 1e-8)
        {
            memcpy(trial.data, beta->data, (size_t)A * (size_t)H * sizeof(double));
            for (int h = 0; h < H; ++h)
                for (int a = 0; a < A; ++a)
                    MATRIX_AT(trial, a, h) += t * step[h * A + a];

            double ft = membership_negative_loglik(V, responsibilities, &trial);
            if (isfinite(ft) && ft <= f0 + c1 * t * gdot)
            {
                memcpy(beta->data, trial.data, (size_t)A * (size_t)H * sizeof(double));
                accepted = true;
                break;
            }
            t *= 0.5;
        }
        if (!accepted)
            break;
    }

    freeMatrix(&phi);
    freeMatrix(&trial);
    Free(grad);
    Free(hess);
    Free(step);
}

static void compute_matrix_component_mult(Matrix *X, Matrix *W, Matrix *P, double *q, double *ballot_loglik)
{
    const int B = X->rows;
    const int C = X->cols;
    const int G = W->cols;
    const double eps = 1e-12;

    for (int b = 0; b < B; ++b)
    {
        double wsum = 0.0;
        for (int g = 0; g < G; ++g)
            wsum += MATRIX_AT_PTR(W, b, g);
        if (!isfinite(wsum) || wsum <= eps)
            wsum = 1.0;

        double ll_b = 0.0;
        for (int c = 0; c < C; ++c)
        {
            double s = 0.0;
            for (int g = 0; g < G; ++g)
                s += MATRIX_AT_PTR(W, b, g) * MATRIX_AT_PTR(P, g, c);
            double marginal = s / wsum;
            if (!isfinite(marginal) || marginal < eps)
                marginal = eps;
            ll_b += MATRIX_AT_PTR(X, b, c) * log(marginal);

            for (int g = 0; g < G; ++g)
            {
                double num = MATRIX_AT_PTR(P, g, c) * MATRIX_AT_PTR(X, b, c);
                double den = s - MATRIX_AT_PTR(P, g, c);
                double val = (isfinite(den) && fabs(den) > eps) ? num / den : 0.0;
                if (!isfinite(val) || val < 0.0)
                    val = 0.0;
                Q_3D(q, b, g, c, G, C) = val;
            }
        }
        ballot_loglik[b] = ll_b;

        for (int g = 0; g < G; ++g)
        {
            double qsum = 0.0;
            for (int c = 0; c < C; ++c)
                qsum += Q_3D(q, b, g, c, G, C);
            if (!isfinite(qsum) || qsum <= eps)
            {
                for (int c = 0; c < C; ++c)
                    Q_3D(q, b, g, c, G, C) = MATRIX_AT_PTR(P, g, c);
            }
            else
            {
                for (int c = 0; c < C; ++c)
                    Q_3D(q, b, g, c, G, C) /= qsum;
            }
        }
    }
}

static void compute_matrix_component_mvn(Matrix *X, Matrix *W, Matrix *P, double *q, double *ballot_loglik,
                                         bool use_cdf, QMethodInput *q_params)
{
    const int B = X->rows;
    const int C = X->cols;
    const int G = W->cols;

    Matrix X_ctx = createMatrix(C, B);
    for (int b = 0; b < B; ++b)
        for (int c = 0; c < C; ++c)
            MATRIX_AT(X_ctx, c, b) = MATRIX_AT_PTR(X, b, c);

    Matrix *prob_by_ballot = (Matrix *)Calloc((size_t)B, Matrix);
    for (int b = 0; b < B; ++b)
        prob_by_ballot[b] = *P;

    EMContext ctx = {0};
    ctx.X = X_ctx;
    ctx.W = *W;
    ctx.probabilities = *P;
    ctx.q = q;
    ctx.ballot_loglik = ballot_loglik;
    ctx.B = (uint32_t)B;
    ctx.G = (uint16_t)G;
    ctx.C = (uint16_t)C;

    uint32_t old_total_ballots = TOTAL_BALLOTS;
    uint16_t old_total_candidates = TOTAL_CANDIDATES;
    uint16_t old_total_groups = TOTAL_GROUPS;
    TOTAL_BALLOTS = (uint32_t)B;
    TOTAL_CANDIDATES = (uint16_t)C;
    TOTAL_GROUPS = (uint16_t)G;

    QMethodInput params = {0};
    if (q_params != NULL)
        params = *q_params;
    params.computeLL = true;
    if (use_cdf)
    {
        if (params.monteCarloIter <= 0)
            params.monteCarloIter = 1000;
        if (!(params.errorThreshold > 0.0) || !isfinite(params.errorThreshold))
            params.errorThreshold = 1e-6;
        if (params.simulationMethod == NULL)
            params.simulationMethod = "genz";
        allocateSeed(&ctx, params);
    }

    double ll = 0.0;
    if (use_cdf)
        computeQMultivariateCDFByBallot(&ctx, prob_by_ballot, params, &ll);
    else
        computeQMultivariatePDFByBallot(&ctx, prob_by_ballot, params, &ll);

    if (ctx.cdf_seeds != NULL)
        Free(ctx.cdf_seeds);
    TOTAL_BALLOTS = old_total_ballots;
    TOTAL_CANDIDATES = old_total_candidates;
    TOTAL_GROUPS = old_total_groups;
    Free(prob_by_ballot);
    freeMatrix(&X_ctx);
}

static void compute_matrix_component_estep(Matrix *X, Matrix *W, Matrix *P, double *q, double *ballot_loglik,
                                           const char *q_method, QMethodInput *q_params)
{
    if (q_method != NULL && strcmp(q_method, "mvn_pdf") == 0)
        compute_matrix_component_mvn(X, W, P, q, ballot_loglik, false, q_params);
    else if (q_method != NULL && strcmp(q_method, "mvn_cdf") == 0)
        compute_matrix_component_mvn(X, W, P, q, ballot_loglik, true, q_params);
    else
        compute_matrix_component_mult(X, W, P, q, ballot_loglik);
}

EMParametricMixtureResult *EM_Algorithm_Parametric_Mixture(
    Matrix *X, Matrix *W, Matrix *V, Matrix *component_prob_init, Matrix *membership_beta_init, int mixture_h,
    const int maxiter, const int miniter, const double maxtime, const double convergence, const double ll_threshold,
    const int maxnewton, const bool verbose, const char *q_method, QMethodInput *q_params,
    const char *adjust_prob_cond_method, bool adjust_prob_cond_every)
{
    (void)adjust_prob_cond_method;
    (void)adjust_prob_cond_every;
    const double eps = 1e-12;
    if (mixture_h < 1)
        error("run_em: Invalid 'mixture'. Must be a positive integer.");

    int B = V->rows;
    int A = V->cols;
    int G = W->cols;
    int C = X->cols;
    int K = mixture_h;

    Matrix *components = alloc_matrix_array(K, G, C);
    Matrix beta = copMatrix(membership_beta_init);
    Matrix prev_beta = createMatrix(A, K - 1);
    double **q_components = (double **)Calloc((size_t)K, double *);
    double **ballot_scores = (double **)Calloc((size_t)K, double *);
    Matrix score = createMatrix(B, mixture_h);
    Matrix phi = createMatrix(B, mixture_h);
    Matrix responsibilities = createMatrix(B, mixture_h);
    double *row_counts = (double *)Calloc((size_t)C, double);
    double *row_scores = (double *)Calloc((size_t)mixture_h, double);
    double *prev_component_prob = (double *)Calloc((size_t)G * (size_t)C * (size_t)K, double);

    for (int k = 0; k < K; ++k)
    {
        memcpy(components[k].data, component_prob_init[k].data, (size_t)G * (size_t)C * sizeof(double));
        normalize_matrix_mixture_component(&components[k]);
        q_components[k] = (double *)Calloc((size_t)B * (size_t)G * (size_t)C, double);
        ballot_scores[k] = (double *)Calloc((size_t)B, double);
    }
    compute_membership_phi(V, &beta, &phi);
    for (int b = 0; b < B; ++b)
        for (int k = 0; k < K; ++k)
            MATRIX_AT(responsibilities, b, k) = MATRIX_AT(phi, b, k);

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double elapsed_total = 0.0;
    double final_ll = -INFINITY;
    double prev_ll = -INFINITY;
    int iter_done = 0;
    int finish_id = 2;

    for (int iter = 1; iter <= maxiter; ++iter)
    {
        iter_done = iter;

        memcpy(prev_beta.data, beta.data, (size_t)A * (size_t)(K - 1) * sizeof(double));
        for (int k = 0; k < K; ++k)
        {
            for (int c = 0; c < C; ++c)
                for (int g = 0; g < G; ++g)
                    prev_component_prob[(size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k)] =
                        MATRIX_AT(components[k], g, c);
        }

        compute_membership_phi(V, &beta, &phi);
        for (int k = 0; k < K; ++k)
        {
            compute_matrix_component_estep(X, W, &components[k], q_components[k], ballot_scores[k], q_method,
                                           q_params);
            for (int b = 0; b < B; ++b)
                MATRIX_AT(score, b, k) = ballot_scores[k][b];
        }

        final_ll = update_matrix_parametric_responsibilities(&score, &phi, &responsibilities, row_scores);

        for (int k = 0; k < K; ++k)
        {
            for (int g = 0; g < G; ++g)
            {
                memset(row_counts, 0, (size_t)C * sizeof(double));
                double den = 0.0;
                for (int b = 0; b < B; ++b)
                {
                    double w = MATRIX_AT(responsibilities, b, k) * MATRIX_AT_PTR(W, b, g);
                    if (w <= 0.0 || !isfinite(w))
                        continue;
                    den += w;
                    for (int c = 0; c < C; ++c)
                        row_counts[c] += w * Q_3D(q_components[k], b, g, c, G, C);
                }
                if (den <= eps || !isfinite(den))
                    continue;
                double rowsum = 0.0;
                for (int c = 0; c < C; ++c)
                {
                    double v = row_counts[c] / den;
                    if (!isfinite(v) || v < eps)
                        v = eps;
                    row_counts[c] = v;
                    rowsum += v;
                }
                if (!isfinite(rowsum) || rowsum <= 0.0)
                    rowsum = 1.0;
                for (int c = 0; c < C; ++c)
                    MATRIX_AT(components[k], g, c) = row_counts[c] / rowsum;
            }
        }

        newton_membership_beta(V, &responsibilities, &beta, maxnewton, 1e-8);

        double delta = 0.0;
        delta = fmax(delta, matrix_max_delta(&beta, &prev_beta));
        for (int k = 0; k < K; ++k)
        {
            for (int c = 0; c < C; ++c)
            {
                for (int g = 0; g < G; ++g)
                {
                    double prev = prev_component_prob[(size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k)];
                    double d = fabs(MATRIX_AT(components[k], g, c) - prev);
                    if (isfinite(d) && d > delta)
                        delta = d;
                }
            }
        }

        if (verbose)
        {
            Rprintf("Iteration %d: parametric mixture log-likelihood = %.4f, delta = %.10f\n", iter, final_ll,
                    delta);
        }

        if (iter % 5 == 0)
            R_CheckUserInterrupt();

        clock_gettime(CLOCK_MONOTONIC_RAW, &now);
        elapsed_total = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
        if (elapsed_total >= maxtime)
        {
            finish_id = 1;
            break;
        }
        if (iter > 1 && iter >= miniter &&
            (delta < convergence || (isfinite(prev_ll) && fabs(final_ll - prev_ll) <= ll_threshold)))
        {
            finish_id = 0;
            break;
        }
        prev_ll = final_ll;
    }

    compute_membership_phi(V, &beta, &phi);
    for (int k = 0; k < K; ++k)
    {
        compute_matrix_component_estep(X, W, &components[k], q_components[k], ballot_scores[k], q_method, q_params);
        for (int b = 0; b < B; ++b)
            MATRIX_AT(score, b, k) = ballot_scores[k][b];
    }
    final_ll = update_matrix_parametric_responsibilities(&score, &phi, &responsibilities, row_scores);
    clock_gettime(CLOCK_MONOTONIC_RAW, &now);
    elapsed_total = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;

    EMParametricMixtureResult *res = (EMParametricMixtureResult *)Calloc(1, EMParametricMixtureResult);
    res->probabilities = alloc_matrix_array(B, G, C);
    res->q = alloc_matrix_array(B, G, C);
    res->expected = alloc_matrix_array(B, G, C);
    res->beta = copMatrix(&beta);
    res->phi = copMatrix(&phi);
    res->component_prob = (double *)Calloc((size_t)G * (size_t)C * (size_t)K, double);

    for (int k = 0; k < K; ++k)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int g = 0; g < G; ++g)
            {
                size_t idx = (size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k);
                res->component_prob[idx] = MATRIX_AT(components[k], g, c);
            }
        }
    }

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int c = 0; c < C; ++c)
            {
                double p = 0.0;
                double q = 0.0;
                for (int k = 0; k < mixture_h; ++k)
                {
                    double r = MATRIX_AT(responsibilities, b, k);
                    p += r * MATRIX_AT(components[k], g, c);
                    q += r * Q_3D(q_components[k], b, g, c, G, C);
                }
                MATRIX_AT(res->probabilities[b], g, c) = p;
                MATRIX_AT(res->q[b], g, c) = q;
                MATRIX_AT(res->expected[b], g, c) = MATRIX_AT_PTR(W, b, g) * q;
            }
        }
    }

    res->responsibilities = responsibilities;
    res->mixture_h = mixture_h;
    res->total_iterations = iter_done;
    res->total_time = elapsed_total;
    res->log_likelihood = final_ll;
    res->finish_id = finish_id;

    for (int k = 0; k < K; ++k)
    {
        freeMatrix(&components[k]);
        Free(q_components[k]);
        Free(ballot_scores[k]);
    }
    Free(components);
    Free(q_components);
    Free(ballot_scores);
    freeMatrix(&beta);
    freeMatrix(&prev_beta);
    freeMatrix(&phi);
    freeMatrix(&score);
    Free(row_counts);
    Free(row_scores);
    Free(prev_component_prob);

    return res;
}

void cleanupParametricMixtureResult(EMParametricMixtureResult *res)
{
    if (res == NULL)
        return;

    if (res->probabilities != NULL)
    {
        int B = res->responsibilities.rows;
        for (int b = 0; b < B; ++b)
            freeMatrix(&res->probabilities[b]);
        Free(res->probabilities);
    }
    if (res->q != NULL)
    {
        int B = res->responsibilities.rows;
        for (int b = 0; b < B; ++b)
            freeMatrix(&res->q[b]);
        Free(res->q);
    }
    if (res->expected != NULL)
    {
        int B = res->responsibilities.rows;
        for (int b = 0; b < B; ++b)
            freeMatrix(&res->expected[b]);
        Free(res->expected);
    }
    if (res->beta.data != NULL)
        freeMatrix(&res->beta);
    if (res->phi.data != NULL)
        freeMatrix(&res->phi);
    if (res->responsibilities.data != NULL)
        freeMatrix(&res->responsibilities);
    if (res->component_prob != NULL)
        Free(res->component_prob);
    Free(res);
}
