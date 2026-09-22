// Joint KL projection for symmetric EM, computed with iterative proportional fitting.

#include "KL.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <math.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#define KL_EPS 1e-12
#define KL_TOL 1e-10
#define KL_MAX_ITER 10000

typedef struct
{
    const Matrix *X;
    const Matrix *W;
    Matrix *q_forward_bgc;
    Matrix *q_reverse_bcg;
    double *q_forward;
    double *q_reverse;
    int G;
    int C;
    bool x_is_cb;
} KLJointInput;

static inline double kl_get_x(const KLJointInput *input, int b, int c)
{
    if (input->x_is_cb)
        return MATRIX_AT_PTR(input->X, c, b);
    return MATRIX_AT_PTR(input->X, b, c);
}

static inline double kl_get_w(const KLJointInput *input, int b, int g)
{
    return MATRIX_AT_PTR(input->W, b, g);
}

static inline double kl_get_q_forward(const KLJointInput *input, int b, int g, int c)
{
    if (input->q_forward_bgc != NULL)
        return MATRIX_AT(input->q_forward_bgc[b], g, c);
    return Q_3D(input->q_forward, b, g, c, input->G, input->C);
}

static inline double kl_get_q_reverse(const KLJointInput *input, int b, int c, int g)
{
    if (input->q_reverse_bcg != NULL)
        return MATRIX_AT(input->q_reverse_bcg[b], c, g);
    return Q_3D(input->q_reverse, b, c, g, input->C, input->G);
}

static inline void kl_set_q_forward(const KLJointInput *input, int b, int g, int c, double value)
{
    if (input->q_forward_bgc != NULL)
        MATRIX_AT(input->q_forward_bgc[b], g, c) = value;
    else
        Q_3D(input->q_forward, b, g, c, input->G, input->C) = value;
}

static inline void kl_set_q_reverse(const KLJointInput *input, int b, int c, int g, double value)
{
    if (input->q_reverse_bcg != NULL)
        MATRIX_AT(input->q_reverse_bcg[b], c, g) = value;
    else
        Q_3D(input->q_reverse, b, c, g, input->C, input->G) = value;
}

static int project_ballot_kl(const KLJointInput *input, int b)
{
    const int G = input->G;
    const int C = input->C;
    const bool use_reverse = input->q_reverse_bcg != NULL || input->q_reverse != NULL;
    double sum_w = 0.0;
    double sum_x = 0.0;

    for (int g = 0; g < G; ++g)
    {
        const double value = kl_get_w(input, b, g);
        if (!isfinite(value) || value < 0.0)
            return -100;
        sum_w += value;
    }
    for (int c = 0; c < C; ++c)
    {
        const double value = kl_get_x(input, b, c);
        if (!isfinite(value) || value < 0.0)
            return -100;
        sum_x += value;
    }

    if (!isfinite(sum_w) || !isfinite(sum_x))
        return -100;

    const double total = fmax(fabs(sum_w), fabs(sum_x));
    const double rel = total > 0.0 ? fabs(sum_w - sum_x) / total : 0.0;
    const double alpha = rel > 1e-12 && sum_w > 0.0 ? sum_x / sum_w : 1.0;
    const double tol = KL_TOL * fmax(1.0, fmax(sum_w, sum_x));
    double *w = Calloc(G, double);
    double *x = Calloc(C, double);
    double *z = Calloc((size_t)G * (size_t)C, double);

    for (int g = 0; g < G; ++g)
        w[g] = kl_get_w(input, b, g) * alpha;
    for (int c = 0; c < C; ++c)
        x[c] = kl_get_x(input, b, c);

    if (sum_x == 0.0)
    {
        Free(w);
        Free(x);
        Free(z);
        return 0;
    }

    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            const int gc = g * C + c;
            if (w[g] <= 0.0 || x[c] <= 0.0)
            {
                z[gc] = 0.0;
                continue;
            }

            double z_forward = w[g] * kl_get_q_forward(input, b, g, c);
            if (!isfinite(z_forward) || z_forward < KL_EPS)
                z_forward = KL_EPS;
            z[gc] = z_forward;
            if (use_reverse)
            {
                double z_reverse = x[c] * kl_get_q_reverse(input, b, c, g);
                if (!isfinite(z_reverse) || z_reverse < KL_EPS)
                    z_reverse = KL_EPS;
                // The equal-weight joint KL objective is an I-projection from this
                // geometric mean of the forward and reverse estimated counts.
                z[gc] = exp(0.5 * (log(z_forward) + log(z_reverse)));
            }
        }
    }

    bool converged = false;
    for (int iter = 0; iter < KL_MAX_ITER; ++iter)
    {
        for (int g = 0; g < G; ++g)
        {
            double row_sum = 0.0;
            for (int c = 0; c < C; ++c)
                row_sum += z[g * C + c];

            if (w[g] > 0.0 && (!isfinite(row_sum) || row_sum <= 0.0))
                goto cleanup;

            const double scale = row_sum > 0.0 ? w[g] / row_sum : 0.0;
            for (int c = 0; c < C; ++c)
                z[g * C + c] *= scale;
        }

        for (int c = 0; c < C; ++c)
        {
            double col_sum = 0.0;
            for (int g = 0; g < G; ++g)
                col_sum += z[g * C + c];

            if (x[c] > 0.0 && (!isfinite(col_sum) || col_sum <= 0.0))
                goto cleanup;

            const double scale = col_sum > 0.0 ? x[c] / col_sum : 0.0;
            for (int g = 0; g < G; ++g)
                z[g * C + c] *= scale;
        }

        double error = 0.0;
        for (int g = 0; g < G; ++g)
        {
            double row_sum = 0.0;
            for (int c = 0; c < C; ++c)
                row_sum += z[g * C + c];
            error = fmax(error, fabs(row_sum - w[g]));
        }
        for (int c = 0; c < C; ++c)
        {
            double col_sum = 0.0;
            for (int g = 0; g < G; ++g)
                col_sum += z[g * C + c];
            error = fmax(error, fabs(col_sum - x[c]));
        }

        if (error <= tol)
        {
            converged = true;
            break;
        }
    }

    if (!converged)
        goto cleanup;

    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            const double value = z[g * C + c];
            if (!isfinite(value) || value < 0.0)
                goto cleanup;
        }
    }

    for (int g = 0; g < G; ++g)
    {
        for (int c = 0; c < C; ++c)
        {
            const double value = z[g * C + c];
            if (w[g] > 0.0)
                kl_set_q_forward(input, b, g, c, value / w[g]);
            if (use_reverse && x[c] > 0.0)
                kl_set_q_reverse(input, b, c, g, value / x[c]);
        }
    }

    Free(w);
    Free(x);
    Free(z);
    return 0;

cleanup:
    Free(w);
    Free(x);
    Free(z);
    return -100;
}

int KL_project(const Matrix *X, const Matrix *W, Matrix *q_forward, int b)
{
    if (X == NULL || W == NULL || q_forward == NULL)
        return -1;
    if (X->rows != W->rows || b < 0 || b >= X->rows)
        return -1;

    KLJointInput input = {0};
    input.X = X;
    input.W = W;
    input.q_forward_bgc = q_forward;
    input.G = W->cols;
    input.C = X->cols;
    input.x_is_cb = false;
    return project_ballot_kl(&input, b);
}

int KL_project_ctx(EMContext *ctx, int b)
{
    if (ctx == NULL || b < 0 || b >= (int)ctx->B)
        return -1;

    KLJointInput input = {0};
    input.X = &ctx->X;
    input.W = &ctx->W;
    input.q_forward = ctx->q;
    input.G = (int)ctx->G;
    input.C = (int)ctx->C;
    input.x_is_cb = true;
    return project_ballot_kl(&input, b);
}

int KL_joint_symmetric(const Matrix *X, const Matrix *W, Matrix *q_forward, Matrix *q_reverse, int b)
{
    if (X == NULL || W == NULL || q_forward == NULL || q_reverse == NULL)
        return -1;
    if (X->rows != W->rows || b < 0 || b >= X->rows)
        return -1;

    KLJointInput input = {0};
    input.X = X;
    input.W = W;
    input.q_forward_bgc = q_forward;
    input.q_reverse_bcg = q_reverse;
    input.G = W->cols;
    input.C = X->cols;
    input.x_is_cb = false;
    return project_ballot_kl(&input, b);
}

int KL_joint_symmetric_ctx(EMContext *ctx_forward, EMContext *ctx_reverse, int b)
{
    if (ctx_forward == NULL || ctx_reverse == NULL)
        return -1;
    if (ctx_forward->B != ctx_reverse->B || ctx_forward->G != ctx_reverse->C || ctx_forward->C != ctx_reverse->G)
        return -1;
    if (b < 0 || b >= (int)ctx_forward->B)
        return -1;

    KLJointInput input = {0};
    input.X = &ctx_forward->X;
    input.W = &ctx_forward->W;
    input.q_forward = ctx_forward->q;
    input.q_reverse = ctx_reverse->q;
    input.G = (int)ctx_forward->G;
    input.C = (int)ctx_forward->C;
    input.x_is_cb = true;
    return project_ballot_kl(&input, b);
}
