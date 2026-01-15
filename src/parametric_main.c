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
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <R_ext/Utils.h> // for R_CheckUserInterrupt
#include <Rinternals.h>
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

void E_step(Matrix *X, Matrix *W, Matrix *V, EMBuffers *buf)
{

    int B = buf->B, G = buf->G, C = buf->C;

    // ---- Get the probabilities
    getProbability(buf, V, &buf->alpha, &buf->beta);

    for (int b = 0; b < B; ++b)
    {
        for (int c = 0; c < C; ++c)
        {
            double acc = 0.0;
            for (int g = 0; g < G; ++g)
            {
                acc += MATRIX_AT_PTR(W, b, g) * MATRIX_AT(buf->prob[b], g, c);
            }
            MATRIX_AT(buf->S_bc, b, c) = acc;
        }
    }

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
                double v = n / d;
                MATRIX_AT(buf->q_bgc[b], g, c) = v;
                denom += v;
            }
            // --- normalize
            for (int c = 0; c < C; ++c)
            {
                MATRIX_AT(buf->q_bgc[b], g, c) /= denom;
            }
        }
    }
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

Matrix *EM_Algorithm(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                     const double maxtime, const double ll_threshold, const int maxnewton, const bool verbose,
                     double *out_elapsed, int *total_iterations, double *logLikelihood,
                     Matrix **out_q, Matrix **out_expected,
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
        E_step(X, W, V, &buf);
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
    E_step(X, W, V, &buf);
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
