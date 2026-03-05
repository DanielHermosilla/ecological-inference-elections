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

#include "saddlepoint.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define SP_MAX_NEWTON 120
#define SP_LINESEARCH_STEPS 16
#define SP_TOL 1e-8
#define SP_FINAL_RES_TOL 1e-6
#define SP_T_CLAMP 40.0
#define SP_RIDGE_START 1e-10
#define SP_PROB_EPS 1e-15
#define SP_SECOND_ORDER_EPS 1e-12

typedef struct
{
    int d;
    int *cand_idx;
    double *t;
    double *t_trial;
    double *grad;
    double *grad_trial;
    double *rhs;
    double *step;
    double *hess;
    double *mat;
    double *q;
} SaddleScratch;

static inline double clampDouble(double x, double lo, double hi)
{
    if (x < lo)
        return lo;
    if (x > hi)
        return hi;
    return x;
}

static SaddleScratch initSaddleScratch(int d)
{
    SaddleScratch S = {0};
    S.d = d;
    if (d <= 0)
        return S;

    S.cand_idx = (int *)Calloc((size_t)d, int);
    S.t = (double *)Calloc((size_t)d, double);
    S.t_trial = (double *)Calloc((size_t)d, double);
    S.grad = (double *)Calloc((size_t)d, double);
    S.grad_trial = (double *)Calloc((size_t)d, double);
    S.rhs = (double *)Calloc((size_t)d, double);
    S.step = (double *)Calloc((size_t)d, double);
    S.hess = (double *)Calloc((size_t)d * (size_t)d, double);
    S.mat = (double *)Calloc((size_t)d * (size_t)d, double);
    S.q = (double *)Calloc((size_t)d, double);
    return S;
}

static void freeSaddleScratch(SaddleScratch *S)
{
    if (S == NULL)
        return;
    if (S->t)
        Free(S->t);
    if (S->cand_idx)
        Free(S->cand_idx);
    if (S->t_trial)
        Free(S->t_trial);
    if (S->grad)
        Free(S->grad);
    if (S->grad_trial)
        Free(S->grad_trial);
    if (S->rhs)
        Free(S->rhs);
    if (S->step)
        Free(S->step);
    if (S->hess)
        Free(S->hess);
    if (S->mat)
        Free(S->mat);
    if (S->q)
        Free(S->q);
    memset(S, 0, sizeof(*S));
}

static bool choleskyInPlace(double *A, int n, double *logdet)
{
    if (logdet)
        *logdet = 0.0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double sum = A[(size_t)i * (size_t)n + (size_t)j];
            for (int k = 0; k < j; ++k)
            {
                sum -= A[(size_t)i * (size_t)n + (size_t)k] * A[(size_t)j * (size_t)n + (size_t)k];
            }

            if (i == j)
            {
                if (!(sum > 0.0) || !isfinite(sum))
                    return false;
                A[(size_t)i * (size_t)n + (size_t)i] = sqrt(sum);
            }
            else
            {
                double den = A[(size_t)j * (size_t)n + (size_t)j];
                if (!(den > 0.0) || !isfinite(den))
                    return false;
                A[(size_t)i * (size_t)n + (size_t)j] = sum / den;
            }
        }

        for (int j = i + 1; j < n; ++j)
        {
            A[(size_t)i * (size_t)n + (size_t)j] = 0.0;
        }
    }

    if (logdet)
    {
        for (int i = 0; i < n; ++i)
        {
            double lii = A[(size_t)i * (size_t)n + (size_t)i];
            if (!(lii > 0.0) || !isfinite(lii))
                return false;
            *logdet += 2.0 * log(lii);
        }
    }
    return true;
}

static void solveFromCholesky(const double *L, const double *b, int n, double *x)
{
    // Forward solve: L y = b
    for (int i = 0; i < n; ++i)
    {
        double sum = b[i];
        for (int k = 0; k < i; ++k)
        {
            sum -= L[(size_t)i * (size_t)n + (size_t)k] * x[k];
        }
        x[i] = sum / L[(size_t)i * (size_t)n + (size_t)i];
    }

    // Backward solve: L' x = y
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = x[i];
        for (int k = i + 1; k < n; ++k)
        {
            sum -= L[(size_t)k * (size_t)n + (size_t)i] * x[k];
        }
        x[i] = sum / L[(size_t)i * (size_t)n + (size_t)i];
    }
}

static bool solveSPD(const double *A, const double *rhs, int n, double *x, double *work)
{
    if (n == 1)
    {
        double den = A[0];
        if (!(den > 0.0) || !isfinite(den))
            return false;
        x[0] = rhs[0] / den;
        return isfinite(x[0]);
    }

    double ridge = SP_RIDGE_START;
    for (int attempt = 0; attempt < 7; ++attempt)
    {
        memcpy(work, A, (size_t)n * (size_t)n * sizeof(double));
        for (int i = 0; i < n; ++i)
        {
            work[(size_t)i * (size_t)n + (size_t)i] += ridge;
        }

        if (choleskyInPlace(work, n, NULL))
        {
            solveFromCholesky(work, rhs, n, x);
            return true;
        }
        ridge *= 10.0;
    }
    return false;
}

static bool logdetSPD(const double *A, int n, double *logdet, double *work)
{
    if (n == 1)
    {
        if (!(A[0] > 0.0) || !isfinite(A[0]))
            return false;
        *logdet = log(A[0]);
        return true;
    }

    double ridge = SP_RIDGE_START;
    for (int attempt = 0; attempt < 7; ++attempt)
    {
        memcpy(work, A, (size_t)n * (size_t)n * sizeof(double));
        for (int i = 0; i < n; ++i)
        {
            work[(size_t)i * (size_t)n + (size_t)i] += ridge;
        }

        if (choleskyInPlace(work, n, logdet))
        {
            return true;
        }
        ridge *= 10.0;
    }
    return false;
}

static bool buildGradientAndHessian(const Matrix *probabilities, int G, int C, const int *n_by_group, const double *t,
                                    const int *cand_idx, int ref_c, double *grad, double *hess, double *qtmp)
{
    const int d = C - 1;

    for (int g = 0; g < G; ++g)
    {
        int ng = n_by_group[g];
        if (ng <= 0)
            continue;

        double m = MATRIX_AT_PTR(probabilities, g, ref_c);
        if (!isfinite(m) || m < 0.0)
            m = 0.0;

        for (int j = 0; j < d; ++j)
        {
            int c = cand_idx[j];
            double p = MATRIX_AT_PTR(probabilities, g, c);
            if (!isfinite(p) || p < 0.0)
                p = 0.0;

            double et = exp(t[j]);
            qtmp[j] = p * et;
            m += qtmp[j];
        }

        if (!(m > 0.0) || !isfinite(m))
            return false;

        for (int j = 0; j < d; ++j)
        {
            qtmp[j] /= m;
            grad[j] += (double)ng * qtmp[j];
        }

        if (hess != NULL)
        {
            for (int j = 0; j < d; ++j)
            {
                double qj = qtmp[j];
                hess[(size_t)j * (size_t)d + (size_t)j] += (double)ng * qj * (1.0 - qj);
                for (int k = 0; k < j; ++k)
                {
                    double v = -(double)ng * qj * qtmp[k];
                    hess[(size_t)j * (size_t)d + (size_t)k] += v;
                    hess[(size_t)k * (size_t)d + (size_t)j] += v;
                }
            }
        }
    }

    return true;
}

// Second-order correction on the saddlepoint log-PMF.
// This is a Lugannani-Rice-style scalar correction based on cumulants of the projection v'X
// with v aligned to the saddlepoint direction t.
static double secondOrderLogCorrection(const Matrix *probabilities, int G, int C, const int *n_by_group, const double *t,
                                       const int *cand_idx, int ref_c, double *work)
{
    const int d = C - 1;
    if (d <= 0)
        return 0.0;

    double t_norm2 = 0.0;
    for (int j = 0; j < d; ++j)
        t_norm2 += t[j] * t[j];
    if (!(t_norm2 > SP_SECOND_ORDER_EPS) || !isfinite(t_norm2))
        return 0.0;

    double inv_norm = 1.0 / sqrt(t_norm2);
    for (int j = 0; j < d; ++j)
        work[j] = t[j] * inv_norm; // projection direction v

    double k2 = 0.0, k3 = 0.0, k4 = 0.0;
    for (int g = 0; g < G; ++g)
    {
        int ng = n_by_group[g];
        if (ng <= 0)
            continue;

        double pref = MATRIX_AT_PTR(probabilities, g, ref_c);
        if (!isfinite(pref) || pref < 0.0)
            pref = 0.0;

        double m = pref;
        for (int j = 0; j < d; ++j)
        {
            int c = cand_idx[j];
            double p = MATRIX_AT_PTR(probabilities, g, c);
            if (!isfinite(p) || p < 0.0)
                p = 0.0;
            m += p * exp(t[j]);
        }
        if (!(m > 0.0) || !isfinite(m))
            return 0.0;

        double m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0;
        for (int j = 0; j < d; ++j)
        {
            int c = cand_idx[j];
            double p = MATRIX_AT_PTR(probabilities, g, c);
            if (!isfinite(p) || p < 0.0)
                p = 0.0;
            double pi = (p * exp(t[j])) / m; // tilted category probability
            double z = work[j];
            double z2 = z * z;
            m1 += pi * z;
            m2 += pi * z2;
            m3 += pi * z2 * z;
            m4 += pi * z2 * z2;
        }

        double mu2 = m2 - m1 * m1;
        if (!(mu2 > SP_SECOND_ORDER_EPS) || !isfinite(mu2))
            continue;
        double mu3 = m3 - 3.0 * m2 * m1 + 2.0 * m1 * m1 * m1;
        double mu4 = m4 - 4.0 * m3 * m1 + 6.0 * m2 * m1 * m1 - 3.0 * m1 * m1 * m1 * m1;
        double k4_one = mu4 - 3.0 * mu2 * mu2;

        k2 += (double)ng * mu2;
        k3 += (double)ng * mu3;
        k4 += (double)ng * k4_one;
    }

    if (!(k2 > SP_SECOND_ORDER_EPS) || !isfinite(k2))
        return 0.0;

    double delta = k4 / (8.0 * k2 * k2) - 5.0 * k3 * k3 / (24.0 * k2 * k2 * k2);
    if (!isfinite(delta))
        return 0.0;

    double corr = 1.0 + delta;
    if (!isfinite(corr) || corr <= SP_SECOND_ORDER_EPS)
        corr = SP_SECOND_ORDER_EPS;

    return log(corr);
}

static void setGroupToPriorOrUniform(const Matrix *probabilities, double *q, int b, int g, int G, int C)
{
    double prior_sum = 0.0;
    for (int c = 0; c < C; ++c)
    {
        double pgc = MATRIX_AT_PTR(probabilities, g, c);
        if (isfinite(pgc) && pgc > 0.0)
            prior_sum += pgc;
    }

    if (prior_sum > 0.0)
    {
        for (int c = 0; c < C; ++c)
        {
            double pgc = MATRIX_AT_PTR(probabilities, g, c);
            double qgc = (isfinite(pgc) && pgc > 0.0) ? (pgc / prior_sum) : 0.0;
            Q_3D(q, b, g, c, G, C) = qgc;
        }
        return;
    }

    double u = 1.0 / (double)C;
    for (int c = 0; c < C; ++c)
        Q_3D(q, b, g, c, G, C) = u;
}

static bool setGroupTiltedProbabilities(const Matrix *probabilities, const SaddleScratch *S, const int *idx_by_c, int b,
                                        int g, int G, int C, int ref_c, double *q)
{
    const int d = C - 1;
    if (S == NULL || idx_by_c == NULL || d <= 0)
        return false;

    double denom = MATRIX_AT_PTR(probabilities, g, ref_c);
    if (!isfinite(denom) || denom < 0.0)
        denom = 0.0;

    for (int j = 0; j < d; ++j)
    {
        int c = S->cand_idx[j];
        double p = MATRIX_AT_PTR(probabilities, g, c);
        if (!isfinite(p) || p < 0.0)
            p = 0.0;
        denom += p * exp(S->t[j]);
    }
    if (!(denom > 0.0) || !isfinite(denom))
        return false;

    double qsum = 0.0;
    for (int c = 0; c < C; ++c)
    {
        double p = MATRIX_AT_PTR(probabilities, g, c);
        if (!isfinite(p) || p < 0.0)
            p = 0.0;

        double qgc = 0.0;
        if (c == ref_c)
        {
            qgc = p / denom;
        }
        else
        {
            int j = idx_by_c[c];
            if (j < 0 || j >= d)
                return false;
            qgc = p * exp(S->t[j]) / denom;
        }

        if (!isfinite(qgc) || qgc < 0.0)
            qgc = 0.0;
        Q_3D(q, b, g, c, G, C) = qgc;
        qsum += qgc;
    }

    if (!(qsum > 0.0) || !isfinite(qsum))
        return false;

    for (int c = 0; c < C; ++c)
        Q_3D(q, b, g, c, G, C) /= qsum;
    return true;
}

static int chooseReferenceCandidate(const int *x, int C)
{
    int ref = 0;
    int best = x[0];
    for (int c = 1; c < C; ++c)
    {
        if (x[c] > best)
        {
            best = x[c];
            ref = c;
        }
    }
    return ref;
}

static double weightedAverageProb(const Matrix *probabilities, const int *n_by_group, int G, int c, int n_total)
{
    if (n_total <= 0)
        return 0.0;

    double num = 0.0;
    for (int g = 0; g < G; ++g)
    {
        int ng = n_by_group[g];
        if (ng <= 0)
            continue;
        double p = MATRIX_AT_PTR(probabilities, g, c);
        if (!isfinite(p) || p < 0.0)
            p = 0.0;
        num += (double)ng * p;
    }
    return num / (double)n_total;
}

static void initializeSaddlepointStart(const Matrix *probabilities, int G, int C, const int *n_by_group, const int *x,
                                       const int *cand_idx, int ref_c, double *t)
{
    const int d = C - 1;
    int n_total = 0;
    for (int g = 0; g < G; ++g)
        n_total += n_by_group[g];

    double x_ref = (double)x[ref_c] + 0.5;
    if (!(x_ref > 0.0) || !isfinite(x_ref))
        x_ref = 0.5;

    double p_ref = weightedAverageProb(probabilities, n_by_group, G, ref_c, n_total);
    if (!(p_ref > SP_PROB_EPS) || !isfinite(p_ref))
        p_ref = SP_PROB_EPS;

    for (int j = 0; j < d; ++j)
    {
        int c = cand_idx[j];
        double x_c = (double)x[c] + 0.5;
        if (!(x_c > 0.0) || !isfinite(x_c))
            x_c = 0.5;

        double p_c = weightedAverageProb(probabilities, n_by_group, G, c, n_total);
        if (!(p_c > SP_PROB_EPS) || !isfinite(p_c))
            p_c = SP_PROB_EPS;

        double tj = log(x_c / x_ref) - log(p_c / p_ref);
        t[j] = clampDouble(tj, -SP_T_CLAMP, SP_T_CLAMP);
    }
}

static int prepareScaledGroupCounts(const EMContext *ctx, int b, int *counts, int x_total)
{
    const int G = (int)ctx->G;
    double scale = 1.0;
    if (ctx->scale_factors != NULL && isfinite(ctx->scale_factors[b]) && ctx->scale_factors[b] > 0.0)
    {
        scale = ctx->scale_factors[b];
    }

    int total = 0;
    int max_idx = 0;
    double max_raw = -1.0;

    for (int g = 0; g < G; ++g)
    {
        double raw = MATRIX_AT(ctx->W, b, g) * scale;
        if (!isfinite(raw) || raw < 0.0)
            raw = 0.0;

        int ng = (int)llround(raw);
        if (ng < 0)
            ng = 0;

        counts[g] = ng;
        total += ng;

        if (raw > max_raw)
        {
            max_raw = raw;
            max_idx = g;
        }
    }

    int diff = x_total - total;
    if (diff > 0)
    {
        counts[max_idx] += diff;
    }
    else if (diff < 0)
    {
        int remaining = -diff;
        while (remaining > 0)
        {
            int best = -1;
            int best_count = 0;
            for (int g = 0; g < G; ++g)
            {
                if (counts[g] > best_count)
                {
                    best_count = counts[g];
                    best = g;
                }
            }
            if (best < 0 || best_count <= 0)
                break;

            int take = (best_count < remaining) ? best_count : remaining;
            counts[best] -= take;
            remaining -= take;
        }
    }

    int final_total = 0;
    for (int g = 0; g < G; ++g)
        final_total += counts[g];
    return final_total;
}

static double saddlepointLogPMF(const Matrix *probabilities, int G, int C, const int *n_by_group, const int *x,
                                int ref_c, const double *t_start,
                                SaddleScratch *S)
{
    const int d = C - 1;
    if (d <= 0 || S == NULL || S->d != d)
        return -INFINITY;
    if (ref_c < 0 || ref_c >= C)
        return -INFINITY;

    int n_total = 0;
    int x_total = 0;
    for (int g = 0; g < G; ++g)
    {
        if (n_by_group[g] < 0)
            return -INFINITY;
        n_total += n_by_group[g];
    }
    for (int c = 0; c < C; ++c)
    {
        if (x[c] < 0)
            return -INFINITY;
        x_total += x[c];
    }
    if (x_total != n_total)
        return -INFINITY;

    int idx = 0;
    for (int c = 0; c < C; ++c)
    {
        if (c == ref_c)
            continue;
        S->cand_idx[idx++] = c;
    }
    if (idx != d)
        return -INFINITY;

    int sum_rest = 0;
    for (int j = 0; j < d; ++j)
        sum_rest += x[S->cand_idx[j]];
    if (x[ref_c] != (n_total - sum_rest))
        return -INFINITY;

    if (t_start != NULL)
    {
        memcpy(S->t, t_start, (size_t)d * sizeof(double));
    }
    else
    {
        initializeSaddlepointStart(probabilities, G, C, n_by_group, x, S->cand_idx, ref_c, S->t);
    }

    for (int iter = 0; iter < SP_MAX_NEWTON; ++iter)
    {
        memset(S->grad, 0, (size_t)d * sizeof(double));
        memset(S->hess, 0, (size_t)d * (size_t)d * sizeof(double));

        if (!buildGradientAndHessian(probabilities, G, C, n_by_group, S->t, S->cand_idx, ref_c, S->grad, S->hess,
                                     S->q))
            return -INFINITY;

        double max_abs_res = 0.0;
        for (int j = 0; j < d; ++j)
        {
            S->rhs[j] = S->grad[j] - (double)x[S->cand_idx[j]];
            double a = fabs(S->rhs[j]);
            if (a > max_abs_res)
                max_abs_res = a;
        }

        if (max_abs_res < SP_TOL)
            break;

        if (!solveSPD(S->hess, S->rhs, d, S->step, S->mat))
            return -INFINITY;

        bool accepted = false;
        double alpha = 1.0;
        for (int ls = 0; ls < SP_LINESEARCH_STEPS; ++ls)
        {
            for (int j = 0; j < d; ++j)
            {
                double tj = S->t[j] - alpha * S->step[j];
                S->t_trial[j] = clampDouble(tj, -SP_T_CLAMP, SP_T_CLAMP);
            }

            memset(S->grad_trial, 0, (size_t)d * sizeof(double));
                if (!buildGradientAndHessian(probabilities, G, C, n_by_group, S->t_trial, S->cand_idx, ref_c,
                                             S->grad_trial, NULL, S->q))
                {
                    alpha *= 0.5;
                    continue;
            }

            double trial_max_abs = 0.0;
            for (int j = 0; j < d; ++j)
            {
                double a = fabs(S->grad_trial[j] - (double)x[S->cand_idx[j]]);
                if (a > trial_max_abs)
                    trial_max_abs = a;
            }

            if (trial_max_abs < max_abs_res)
            {
                memcpy(S->t, S->t_trial, (size_t)d * sizeof(double));
                accepted = true;
                break;
            }

            alpha *= 0.5;
        }

        if (!accepted)
        {
            for (int j = 0; j < d; ++j)
            {
                double tj = S->t[j] - 0.1 * S->step[j];
                S->t[j] = clampDouble(tj, -SP_T_CLAMP, SP_T_CLAMP);
            }
        }
    }

    memset(S->grad, 0, (size_t)d * sizeof(double));
    memset(S->hess, 0, (size_t)d * (size_t)d * sizeof(double));
    if (!buildGradientAndHessian(probabilities, G, C, n_by_group, S->t, S->cand_idx, ref_c, S->grad, S->hess, S->q))
        return -INFINITY;

    double final_res = 0.0;
    for (int j = 0; j < d; ++j)
    {
        double a = fabs(S->grad[j] - (double)x[S->cand_idx[j]]);
        if (a > final_res)
            final_res = a;
    }
    if (!(final_res <= SP_FINAL_RES_TOL) || !isfinite(final_res))
        return -INFINITY;

    double K = 0.0;
    for (int g = 0; g < G; ++g)
    {
        int ng = n_by_group[g];
        if (ng <= 0)
            continue;

        double m = MATRIX_AT_PTR(probabilities, g, ref_c);
        if (!isfinite(m) || m < 0.0)
            m = 0.0;

        for (int j = 0; j < d; ++j)
        {
            int c = S->cand_idx[j];
            double p = MATRIX_AT_PTR(probabilities, g, c);
            if (!isfinite(p) || p < 0.0)
                p = 0.0;
            m += p * exp(S->t[j]);
        }

        if (!(m > 0.0) || !isfinite(m))
            return -INFINITY;

        K += (double)ng * log(m);
    }

    double logdet = 0.0;
    if (!logdetSPD(S->hess, d, &logdet, S->mat))
        return -INFINITY;

    double t_dot_x = 0.0;
    for (int j = 0; j < d; ++j)
    {
        t_dot_x += S->t[j] * (double)x[S->cand_idx[j]];
    }

    double out = K - t_dot_x - 0.5 * ((double)d * log(2.0 * M_PI) + logdet);
    out += secondOrderLogCorrection(probabilities, G, C, n_by_group, S->t, S->cand_idx, ref_c, S->t_trial);
    if (!isfinite(out))
        return -INFINITY;
    return out;
}

void computeQSaddlepoint(EMContext *ctx, QMethodInput params, double *ll)
{
    if (ll != NULL)
        *ll = 0.0;
    if (ctx->ballot_loglik != NULL)
    {
        memset(ctx->ballot_loglik, 0, (size_t)ctx->B * sizeof(double));
    }

    Matrix *probabilities = &ctx->probabilities;
    const int B = (int)ctx->B;
    const int G = (int)ctx->G;
    const int C = (int)ctx->C;

    if (C < 2 || G < 1 || B < 1)
    {
        return;
    }

    SaddleScratch S = initSaddleScratch(C - 1);
    int *x = (int *)Calloc((size_t)C, int);
    int *n = (int *)Calloc((size_t)G, int);
    int *idx_by_c = (int *)Calloc((size_t)C, int);

    for (int b = 0; b < B; ++b)
    {
        int x_total = 0;
        for (int c = 0; c < C; ++c)
        {
            int v = MATRIX_AT(ctx->intX, c, b);
            if (v < 0)
                v = 0;
            x[c] = v;
            x_total += v;
        }

        prepareScaledGroupCounts(ctx, b, n, x_total);
        int ref_b = chooseReferenceCandidate(x, C);
        double logp = saddlepointLogPMF(probabilities, G, C, n, x, ref_b, NULL, &S);
        bool valid_saddle = isfinite(logp);

        if (params.computeLL && ll != NULL)
        {
            if (ctx->ballot_loglik != NULL)
            {
                ctx->ballot_loglik[b] = isfinite(logp) ? logp : 0.0;
            }
            if (isfinite(logp))
                *ll += logp;
        }

        for (int c = 0; c < C; ++c)
        {
            idx_by_c[c] = -1;
        }
        if (valid_saddle)
        {
            for (int j = 0; j < C - 1; ++j)
            {
                int c = S.cand_idx[j];
                if (c >= 0 && c < C)
                {
                    idx_by_c[c] = j;
                }
            }
        }

        for (int g = 0; g < G; ++g)
        {
            if (n[g] <= 0 || !valid_saddle ||
                !setGroupTiltedProbabilities(probabilities, &S, idx_by_c, b, g, G, C, ref_b, ctx->q))
            {
                setGroupToPriorOrUniform(probabilities, ctx->q, b, g, G, C);
            }
        }
    }

    if (params.computeLL && ll != NULL && (!isfinite(*ll) || isnan(*ll)))
    {
        *ll = 0.0;
    }

    Free(x);
    Free(n);
    Free(idx_by_c);
    freeSaddleScratch(&S);
}
