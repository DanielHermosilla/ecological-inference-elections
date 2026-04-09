#include "main_symmetric.h"
#include "globals.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

static void setGlobalsFromCtx(const EMContext *ctx)
{
    TOTAL_BALLOTS = ctx->B;
    TOTAL_CANDIDATES = ctx->C;
    TOTAL_GROUPS = ctx->G;
    TOTAL_VOTES = (uint32_t)ctx->total_votes;
}

static void applyProbabilityCondition(EMContext *ctx, QMethodInput inputParams, bool force_every)
{
    if (((!force_every) && !inputParams.prob_cond_every) || inputParams.prob_cond == NULL ||
        strlen(inputParams.prob_cond) == 0)
        return;

    if (strcmp(inputParams.prob_cond, "project_lp") == 0)
    {
        projectQ(ctx, inputParams);
    }
    else if (strcmp(inputParams.prob_cond, "lp") == 0)
    {
        for (int b = 0; b < (int)ctx->B; ++b)
            LPW_ctx(ctx, b);
    }

    // Keep q as a proper conditional probability after any adjustment method.
    for (int b = 0; b < (int)ctx->B; ++b)
    {
        for (int g = 0; g < (int)ctx->G; ++g)
        {
            double sum = 0.0;
            for (int c = 0; c < (int)ctx->C; ++c)
            {
                double v = Q_3D(ctx->q, b, g, c, ctx->G, ctx->C);
                if (!isfinite(v) || v < 0.0)
                    v = 0.0;
                Q_3D(ctx->q, b, g, c, ctx->G, ctx->C) = v;
                sum += v;
            }
            if (!isfinite(sum) || sum <= 0.0)
            {
                const double uniform = 1.0 / (double)ctx->C;
                for (int c = 0; c < (int)ctx->C; ++c)
                    Q_3D(ctx->q, b, g, c, ctx->G, ctx->C) = uniform;
            }
            else
            {
                for (int c = 0; c < (int)ctx->C; ++c)
                    Q_3D(ctx->q, b, g, c, ctx->G, ctx->C) /= sum;
            }
        }
    }
}

static void buildReverseMatrices(const EMContext *ctx_forward, Matrix *out_x_reverse, Matrix *out_w_reverse)
{
    const int B = (int)ctx_forward->B;
    const int G = (int)ctx_forward->G;
    const int C = (int)ctx_forward->C;

    *out_x_reverse = createMatrix(G, B); // X_rev = t(W): G x B
    *out_w_reverse = createMatrix(B, C); // W_rev = t(X): B x C

    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
            MATRIX_AT_PTR(out_x_reverse, g, b) = MATRIX_AT(ctx_forward->W, b, g);
        for (int c = 0; c < C; ++c)
            MATRIX_AT_PTR(out_w_reverse, b, c) = MATRIX_AT(ctx_forward->X, c, b);
    }
}

static Matrix buildReverseCustomInitialProb(const Matrix *forward_prob, const EMContext *ctx_forward)
{
    const int G = (int)ctx_forward->G;
    const int C = (int)ctx_forward->C;
    Matrix reverse_prob = createMatrix(C, G); // reverse rows=original candidates, cols=original groups

    double *den = Calloc(G, double);

    for (int g = 0; g < G; ++g)
    {
        double row_sum = 0.0;
        for (int c = 0; c < C; ++c)
            row_sum += MATRIX_AT_PTR(forward_prob, g, c) * ctx_forward->candidates_votes[c];
        den[g] = row_sum;
    }

    for (int c = 0; c < C; ++c)
    {
        double row_sum = 0.0;
        for (int g = 0; g < G; ++g)
        {
            double num = MATRIX_AT_PTR(forward_prob, g, c) * ctx_forward->candidates_votes[c];
            double value = den[g] > 0.0 ? num / den[g] : 0.0;
            MATRIX_AT(reverse_prob, c, g) = value;
            row_sum += value;
        }

        if (!isfinite(row_sum) || row_sum <= 0.0)
        {
            const double uniform = 1.0 / (double)G;
            for (int g = 0; g < G; ++g)
                MATRIX_AT(reverse_prob, c, g) = uniform;
        }
        else
        {
            for (int g = 0; g < G; ++g)
                MATRIX_AT(reverse_prob, c, g) /= row_sum;
        }
    }

    Free(den);
    return reverse_prob;
}

static void averageEstimatedVotesAndUpdateQ(EMContext *ctx_forward, EMContext *ctx_reverse)
{
    const int B = (int)ctx_forward->B;
    const int G = (int)ctx_forward->G;
    const int C = (int)ctx_forward->C;
    const size_t z_size = (size_t)B * (size_t)G * (size_t)C;

    double *z_avg = Calloc(z_size, double);

    // ------------------------------------------------------------
    // 1) Compute symmetric average of estimated counts Z (NOT q)
    //    z_forward[b,g,c] = W[b,g] * q_fwd[b,g,c]
    //    z_reverse[b,g,c] = X[b,c] * q_rev[b,c,g]
    //    z_avg   [b,g,c] = 0.5*(z_forward + z_reverse)
    // ------------------------------------------------------------
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            const double w_bg = MATRIX_AT(ctx_forward->W, b, g);

            for (int c = 0; c < C; ++c)
            {
                const double x_bc = MATRIX_AT(ctx_forward->X, c, b);

                const double z_forward = w_bg * Q_3D(ctx_forward->q, b, g, c, G, C);
                const double z_reverse = x_bc * Q_3D(ctx_reverse->q, b, c, g, ctx_reverse->G, ctx_reverse->C);

                Q_3D(z_avg, b, g, c, G, C) = 0.5 * (z_forward + z_reverse);
            }
        }
    }

    // ------------------------------------------------------------
    // 2) Update forward q from z_avg WITHOUT renormalization
    //    q_fwd[b,g,c] <- z_avg[b,g,c] / W[b,g]  (if W[b,g] > 0)
    //
    //    NOTE: We keep the old safety clamps (finite, nonnegative),
    //          but we DO NOT force sum_c q_fwd[b,g,c] = 1.
    // ------------------------------------------------------------
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            const double den = MATRIX_AT(ctx_forward->W, b, g);

            double sum = 0.0; // (was used only for renormalization)

            for (int c = 0; c < C; ++c)
            {
                double value = den > 0.0 ? Q_3D(z_avg, b, g, c, G, C) / den : 1.0 / (double)C;

                if (!isfinite(value) || value < 0.0)
                    value = 0.0;

                Q_3D(ctx_forward->q, b, g, c, G, C) = value;

                sum += value; // (normalization accumulator)
            }

            // ------------------------------
            // ORIGINAL NORMALIZATION (kept commented)
            // ------------------------------
            if (!isfinite(sum) || sum <= 0.0)
            {
                const double uniform = 1.0 / (double)C;
                for (int c = 0; c < C; ++c)
                    Q_3D(ctx_forward->q, b, g, c, G, C) = uniform;
            }
            else
            {
                for (int c = 0; c < C; ++c)
                    Q_3D(ctx_forward->q, b, g, c, G, C) /= sum;
            }
        }
    }

    // ------------------------------------------------------------
    // 3) Update reverse q from z_avg WITHOUT renormalization
    //    q_rev[b,c,g] <- z_avg[b,g,c] / X[b,c]  (if X[b,c] > 0)
    //
    //    NOTE: We keep the old safety clamps (finite, nonnegative),
    //          but we DO NOT force sum_g q_rev[b,c,g] = 1.
    // ------------------------------------------------------------
    for (int b = 0; b < B; ++b)
    {
        for (int c = 0; c < C; ++c)
        {
            const double den = MATRIX_AT(ctx_forward->X, c, b);

            double sum = 0.0; // (was used only for renormalization)

            for (int g = 0; g < G; ++g)
            {
                double value = den > 0.0 ? Q_3D(z_avg, b, g, c, G, C) / den : 1.0 / (double)G;

                if (!isfinite(value) || value < 0.0)
                    value = 0.0;

                Q_3D(ctx_reverse->q, b, c, g, ctx_reverse->G, ctx_reverse->C) = value;

                sum += value; // (normalization accumulator)
            }

            // ------------------------------
            // ORIGINAL NORMALIZATION (kept commented)
            // ------------------------------
            if (!isfinite(sum) || sum <= 0.0)
            {
                const double uniform = 1.0 / (double)G;
                for (int g = 0; g < G; ++g)
                    Q_3D(ctx_reverse->q, b, c, g, ctx_reverse->G, ctx_reverse->C) = uniform;
            }
            else
            {
                for (int g = 0; g < G; ++g)
                    Q_3D(ctx_reverse->q, b, c, g, ctx_reverse->G, ctx_reverse->C) /= sum;
            }
        }
    }

    Free(z_avg);
}

static bool shouldRunFinalMStep(QMethodInput inputParams)
{
    if (inputParams.prob_cond == NULL)
        return false;
    return strcmp(inputParams.prob_cond, "project_lp") == 0 || strcmp(inputParams.prob_cond, "lp") == 0;
}

bool shouldRunSymmetricEMWeight(const QMethodInput *inputParams)
{
    return inputParams != NULL && inputParams->symmetric && inputParams->symmetric_weight_method != NULL &&
           strcmp(inputParams->symmetric_weight_method, "joint") == 0;
}

void runSymmetricEMWeight(EMContext *ctx_forward, const char *p_method, const char *q_method, const double convergence,
                          const double LLconvergence, const int maxIter, const double maxSeconds, const bool verbose,
                          double *time, int *iterTotal, double *logLLarr, int *finishing_reason, Matrix *probMatrix,
                          QMethodInput *inputParams, QMethodConfig config_forward)
{
    Matrix reverse_x = {NULL, 0, 0};
    Matrix reverse_w = {NULL, 0, 0};
    Matrix reverse_prob_custom = {NULL, 0, 0};
    Matrix old_forward_prob = createMatrix(ctx_forward->G, ctx_forward->C);
    Matrix old_reverse_prob = {NULL, 0, 0};

    QMethodInput reverse_params = *inputParams;
    QMethodConfig config_reverse = {0};
    EMContext *ctx_reverse = NULL;

    buildReverseMatrices(ctx_forward, &reverse_x, &reverse_w);
    ctx_reverse = createEMContext(&reverse_x, &reverse_w, q_method, reverse_params);

    freeMatrix(&reverse_x);
    freeMatrix(&reverse_w);

    if (strcmp(p_method, "custom") == 0)
    {
        reverse_prob_custom = buildReverseCustomInitialProb(probMatrix, ctx_forward);
        getInitialP(ctx_reverse, "custom", &reverse_prob_custom);
        freeMatrix(&reverse_prob_custom);
    }
    else
    {
        getInitialP(ctx_reverse, p_method, probMatrix);
    }

    config_reverse = getQMethodConfig(q_method, reverse_params);
    old_reverse_prob = createMatrix(ctx_reverse->G, ctx_reverse->C);

    if (ctx_forward->B != ctx_reverse->B || ctx_forward->G != ctx_reverse->C || ctx_forward->C != ctx_reverse->G)
    {
        cleanup(ctx_reverse);
        setGlobalsFromCtx(ctx_forward);
        freeMatrix(&old_forward_prob);
        freeMatrix(&old_reverse_prob);
        error("Symmetric EM weight: incompatible forward/reverse dimensions.");
    }

    struct timespec iter_start, iter_end;
    double elapsed_total = 0.0;
    double oldLL_forward = -DBL_MAX;
    double oldLL_reverse = -DBL_MAX;
    double newLL_forward = 0.0;
    double newLL_reverse = 0.0;
    bool converged = false;
    bool timeout_reached = false;

    for (int i = 0; i < maxIter; ++i)
    {
        clock_gettime(CLOCK_MONOTONIC, &iter_start);
        *iterTotal = i;
        ctx_forward->iteration = i;
        ctx_reverse->iteration = i;

        const bool has_prob_cond = inputParams->prob_cond != NULL && strlen(inputParams->prob_cond) > 0;
        const bool run_prob_cond_each_iter = has_prob_cond && inputParams->prob_cond_every;
        const bool is_lp = has_prob_cond && strcmp(inputParams->prob_cond, "lp") == 0;
        const bool is_project_lp = has_prob_cond && strcmp(inputParams->prob_cond, "project_lp") == 0;

        setGlobalsFromCtx(ctx_forward);
        config_forward.computeQ(ctx_forward, config_forward.params, &newLL_forward);
        if (run_prob_cond_each_iter && is_project_lp)
            projectQ(ctx_forward, *inputParams);

        setGlobalsFromCtx(ctx_reverse);
        config_reverse.computeQ(ctx_reverse, config_reverse.params, &newLL_reverse);
        if (run_prob_cond_each_iter && is_project_lp)
            projectQ(ctx_reverse, *inputParams);

        if (run_prob_cond_each_iter && is_lp)
        {
            for (int b = 0; b < (int)ctx_forward->B; ++b)
            {
                int status = LPW_joint_symmetric_ctx(ctx_forward, ctx_reverse, b);
                if (status != 0)
                {
                    // Safety fallback: keep the original behavior if the joint LP fails.
                    LPW_ctx(ctx_forward, b);
                    LPW_ctx(ctx_reverse, b);
                }
            }
        }

        averageEstimatedVotesAndUpdateQ(ctx_forward, ctx_reverse);

        memcpy(old_forward_prob.data, ctx_forward->probabilities.data,
               sizeof(double) * old_forward_prob.rows * old_forward_prob.cols);
        memcpy(old_reverse_prob.data, ctx_reverse->probabilities.data,
               sizeof(double) * old_reverse_prob.rows * old_reverse_prob.cols);

        setGlobalsFromCtx(ctx_forward);
        getP(ctx_forward);
        setGlobalsFromCtx(ctx_reverse);
        getP(ctx_reverse);

        *logLLarr = 0.5 * (newLL_forward + newLL_reverse);

        if (verbose)
        {
            Rprintf("\n----------\nIteration: %d\nProbability matrix:\n", i + 1);
            setGlobalsFromCtx(ctx_forward);
            printMatrix(&ctx_forward->probabilities);
            Rprintf("Forward log-likelihood: %f\n", newLL_forward);
            Rprintf("Reverse log-likelihood: %f\n", newLL_reverse);
            if (i != 0)
            {
                Rprintf("Delta forward log-likelihood: %f\n", fabs(newLL_forward - oldLL_forward));
                Rprintf("Delta reverse log-likelihood: %f\n", fabs(newLL_reverse - oldLL_reverse));
            }
        }

        bool decreasing_forward = oldLL_forward > newLL_forward;
        bool decreasing_reverse = oldLL_reverse > newLL_reverse;
        bool early_stop_forward = decreasing_forward && strcmp(q_method, "exact") == 0;
        bool early_stop_reverse = decreasing_reverse && strcmp(q_method, "exact") == 0;

        bool converged_forward = fabs(newLL_forward - oldLL_forward) < LLconvergence ||
                                 convergeMatrix(&old_forward_prob, &ctx_forward->probabilities, convergence) ||
                                 early_stop_forward;
        bool converged_reverse = fabs(newLL_reverse - oldLL_reverse) < LLconvergence ||
                                 convergeMatrix(&old_reverse_prob, &ctx_reverse->probabilities, convergence) ||
                                 early_stop_reverse;

        if (i >= 1 && i >= inputParams->miniter && converged_forward && converged_reverse)
        {
            converged = true;
            *finishing_reason = 0;
            break;
        }

        clock_gettime(CLOCK_MONOTONIC, &iter_end);
        elapsed_total += (iter_end.tv_sec - iter_start.tv_sec) + (iter_end.tv_nsec - iter_start.tv_nsec) / 1e9;
        R_CheckUserInterrupt();

        if (verbose)
            Rprintf("Elapsed time: %f\n----------\n", elapsed_total);

        if (elapsed_total >= maxSeconds)
        {
            timeout_reached = true;
            *finishing_reason = 1;
            break;
        }

        oldLL_forward = newLL_forward;
        oldLL_reverse = newLL_reverse;
    }

    if (!converged && !timeout_reached)
    {
        if (verbose)
            Rprintf("Maximum number of iterations reached without convergence.\n");
        *finishing_reason = 2;
    }

    setGlobalsFromCtx(ctx_forward);
    config_forward.computeQ(ctx_forward, config_forward.params, &newLL_forward);

    setGlobalsFromCtx(ctx_reverse);
    config_reverse.computeQ(ctx_reverse, config_reverse.params, &newLL_reverse);

    averageEstimatedVotesAndUpdateQ(ctx_forward, ctx_reverse);

    if (shouldRunFinalMStep(*inputParams))
    {
        if (strcmp(inputParams->prob_cond, "lp") == 0)
        {
            for (int b = 0; b < (int)ctx_forward->B; ++b)
            {
                int status = LPW_joint_symmetric_ctx(ctx_forward, ctx_reverse, b);
                if (status != 0)
                {
                    // Safety fallback: keep the original behavior if the joint LP fails.
                    LPW_ctx(ctx_forward, b);
                    LPW_ctx(ctx_reverse, b);
                }
            }
        }
        else
        {
            applyProbabilityCondition(ctx_forward, *inputParams, true);
            applyProbabilityCondition(ctx_reverse, *inputParams, true);
        }
    }

    setGlobalsFromCtx(ctx_forward);
    getP(ctx_forward);
    setGlobalsFromCtx(ctx_reverse);
    getP(ctx_reverse);

    setGlobalsFromCtx(ctx_forward);
    getPredictedVotes(ctx_forward);

    *logLLarr = 0.5 * (newLL_forward + newLL_reverse);
    *time = elapsed_total;

    cleanup(ctx_reverse);
    setGlobalsFromCtx(ctx_forward);

    freeMatrix(&old_forward_prob);
    freeMatrix(&old_reverse_prob);
}
