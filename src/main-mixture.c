#include "globals.h"
#include "main.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h>
#include <Rinternals.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#ifndef CLOCK_MONOTONIC_RAW
#define CLOCK_MONOTONIC_RAW 4
#endif

static double log_sum_exp_weighted_phi(const double *d_b, const double *phi, int K)
{
    // Numerically stable version of:
    //   log(sum_k exp(d_{bk}) * phi_k),
    // where d_{bk} = log P(X_b | V_{bk}=1; p).
    const double eps = 1e-300;
    double max_term = -INFINITY;

    for (int k = 0; k < K; ++k)
    {
        if (!isfinite(d_b[k]))
            continue;
        double phi_k = (isfinite(phi[k]) && phi[k] > eps) ? phi[k] : eps;
        double term = d_b[k] + log(phi_k);
        if (term > max_term)
            max_term = term;
    }

    if (!isfinite(max_term))
        return -INFINITY;

    double sum_exp = 0.0;
    for (int k = 0; k < K; ++k)
    {
        if (!isfinite(d_b[k]))
            continue;
        double phi_k = (isfinite(phi[k]) && phi[k] > eps) ? phi[k] : eps;
        double term = d_b[k] + log(phi_k);
        sum_exp += exp(term - max_term);
    }

    if (!(sum_exp > 0.0) || !isfinite(sum_exp))
        return -INFINITY;

    return max_term + log(sum_exp);
}

static void print_mixture_state(int iter, double loglik, EMContext **components, const double *phi, int K)
{
    Rprintf("\n----------\nIteration: %d\nLog-likelihood: %.10f\n", iter, loglik);
    for (int k = 0; k < K; ++k)
    {
        Rprintf("Component %d probability matrix (GxC):\n", k + 1);
        printMatrix(&components[k]->probabilities);
    }

    Rprintf("Phi vector (K=%d): [", K);
    for (int k = 0; k < K; ++k)
    {
        if (k < K - 1)
        {
            Rprintf("%4.3f, ", phi[k]);
        }
        else
        {
            Rprintf("%4.3f", phi[k]);
        }
    }
    Rprintf("]\n");
}

static int compare_components_lex(const EMContext *a, const EMContext *b, double phi_a, double phi_b)
{
    const double tol = 1e-12;
    // for (int g = 0; g < a->G; ++g)
    // {
    //     for (int c = 0; c < a->C; ++c)
    //     {
    //         double da = MATRIX_AT(a->probabilities, g, c);
    //         double db = MATRIX_AT(b->probabilities, g, c);
    //         double d = da - db;
    //         if (fabs(d) > tol)
    //         {
    //             return (d < 0.0) ? -1 : 1;
    //         }
    //     }
    // }
    if (fabs(phi_a - phi_b) > tol)
    {
        return (phi_a < phi_b) ? -1 : 1;
    }
    return 0;
}

static void canonicalize_mixture_labels(EMContext **components, double *phi, int K)
{
    if (K <= 1)
    {
        return;
    }

    int *order = (int *)Calloc((size_t)K, int);
    EMContext **components_sorted = (EMContext **)Calloc((size_t)K, EMContext *);
    double *phi_sorted = (double *)Calloc((size_t)K, double);

    for (int i = 0; i < K; ++i)
    {
        order[i] = i;
    }

    // Stable insertion sort by a deterministic lexicographic rule on P_k and phi_k.
    for (int i = 1; i < K; ++i)
    {
        int key = order[i];
        int j = i - 1;
        while (j >= 0 && compare_components_lex(components[key], components[order[j]], phi[key], phi[order[j]]) < 0)
        {
            order[j + 1] = order[j];
            --j;
        }
        order[j + 1] = key;
    }

    for (int i = 0; i < K; ++i)
    {
        components_sorted[i] = components[order[i]];
        phi_sorted[i] = phi[order[i]];
    }
    for (int i = 0; i < K; ++i)
    {
        components[i] = components_sorted[i];
        phi[i] = phi_sorted[i];
    }

    Free(order);
    Free(components_sorted);
    Free(phi_sorted);
}

EMMixtureResult *EMAlgoritmMixture(Matrix *X, Matrix *W, const char *p_method, const char *q_method,
                                   const double convergence, const double LLconvergence, const int maxIter,
                                   const double maxSeconds, const bool verbose, Matrix *probMatrix,
                                   QMethodInput *inputParams, int mixture_h)
{
    const double eps = 1e-12;
    if (mixture_h < 1)
    {
        error("run_em: Invalid 'mixture'. Must be a positive integer.");
    }

    if (mixture_h == 1)
    {
        double t = 0.0, ll = 0.0;
        int it = 0, finish = 2;
        EMContext *ctx = EMAlgoritm(X, W, p_method, q_method, convergence, LLconvergence, maxIter, maxSeconds, verbose,
                                    &t, &it, &ll, &finish, probMatrix, inputParams);

        EMMixtureResult *res = Calloc(1, EMMixtureResult);
        res->probabilities = copMatrix(&ctx->probabilities);
        size_t N = (size_t)ctx->B * (size_t)ctx->G * (size_t)ctx->C;
        res->q = (double *)Calloc(N, double);
        memcpy(res->q, ctx->q, N * sizeof(double));
        res->predicted_votes = (double *)Calloc(N, double);
        memcpy(res->predicted_votes, ctx->predicted_votes, N * sizeof(double));
        res->phi = (double *)Calloc(1, double);
        res->phi[0] = 1.0;
        res->responsibilities = createMatrix(ctx->B, 1);
        for (int b = 0; b < ctx->B; ++b)
        {
            MATRIX_AT(res->responsibilities, b, 0) = 1.0;
        }
        res->component_prob = (double *)Calloc((size_t)ctx->G * (size_t)ctx->C, double);
        memcpy(res->component_prob, ctx->probabilities.data, (size_t)ctx->G * (size_t)ctx->C * sizeof(double));
        res->mixture_h = 1;
        res->total_iterations = it;
        res->total_time = t;
        res->log_likelihood = ll;
        res->finish_id = finish;

        cleanup(ctx);
        return res;
    }
    //
    // ---- Timer ---- //
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double elapsed_total = 0.0;

    // ---- Allocate per-component contexts ---- //
    QMethodInput ctx_params = *inputParams;
    ctx_params.computeLL = true;

    EMContext **components = (EMContext **)Calloc((size_t)mixture_h, EMContext *);
    if (components == NULL)
    {
        error("Failed to allocate mixture contexts.");
    }

    for (int k = 0; k < mixture_h; ++k)
    {
        components[k] = createEMContext(X, W, q_method, ctx_params);
    }

    // ---- Initialize probabilities for each component ---- //
    getInitialP(components[0], p_method, probMatrix);
    normalizeProbabilityRows(&components[0]->probabilities);
    for (int k = 1; k < mixture_h; ++k)
    {
        memcpy(components[k]->probabilities.data, components[0]->probabilities.data,
               (size_t)components[0]->G * (size_t)components[0]->C * sizeof(double));
    }
    GetRNGstate();
    for (int k = 1; k < mixture_h; ++k)
    {
        for (int g = 0; g < components[k]->G; ++g)
        {
            for (int c = 0; c < components[k]->C; ++c)
            {
                double jit = 0.9 + 0.2 * unif_rand();
                MATRIX_AT(components[k]->probabilities, g, c) *= jit;
            }
        }
        normalizeProbabilityRows(&components[k]->probabilities);
    }
    PutRNGstate();

    const int B = components[0]->B;
    const int G = components[0]->G;
    const int C = components[0]->C;

    double *phi = (double *)Calloc((size_t)mixture_h, double);
    double *phi_new = (double *)Calloc((size_t)mixture_h, double);
    Matrix responsibilities = createMatrix(B, mixture_h); // B x K
    Matrix score = createMatrix(B, mixture_h);            // B x K
    double *row_counts = (double *)Calloc((size_t)C, double);
    double *d_row = (double *)Calloc((size_t)mixture_h, double);
    double *prev_component_prob = (double *)Calloc((size_t)G * (size_t)C * (size_t)mixture_h, double);
    double *prev_phi_params = (double *)Calloc((size_t)mixture_h, double);

    for (int k = 0; k < mixture_h; ++k)
    {
        phi[k] = 1.0 / (double)mixture_h;
    }

    QMethodConfig config = getQMethodConfig(q_method, *inputParams);
    QMethodInput mix_params = config.params;
    mix_params.computeLL = true;

    int finish_id = 2;
    int iter_done = 0;
    double final_ll = NA_REAL;
    double prev_ll = -INFINITY;

    // ---- Mixture EM loop ---- //
    for (int iter = 1; iter <= maxIter; ++iter)
    {
        iter_done = iter;
        int fallback_underflow_count = 0;
        int fallback_phi_prior_count = 0;
        int fallback_uniform_count = 0;

        // Snapshot previous canonical parameters for convergence checking.
        for (int k = 0; k < mixture_h; ++k)
        {
            prev_phi_params[k] = phi[k];
            for (int c = 0; c < C; ++c)
            {
                for (int g = 0; g < G; ++g)
                {
                    prev_component_prob[(size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k)] =
                        MATRIX_AT(components[k]->probabilities, g, c);
                }
            }
        }

        // E-step (per component): compute q_{bkgc} and d_{bk} = log(a_{bk})
        for (int k = 0; k < mixture_h; ++k)
        {
            EMContext *ctxk = components[k];
            double ll_dummy = 0.0;
            config.computeQ(ctxk, mix_params, &ll_dummy); // Aproximation
            if (inputParams->prob_cond_every)
            {
                if (strcmp(inputParams->prob_cond, "project_lp") == 0)
                {
                    projectQ(ctxk, *inputParams);
                }
                else if (strcmp(inputParams->prob_cond, "lp") == 0)
                {
                    for (int b = 0; b < B; ++b)
                    {
                        LPW_ctx(ctxk, b);
                    }
                }
            }

            for (int b = 0; b < B; ++b)
            {
                double d_bk = ctxk->ballot_loglik[b];
                if (!isfinite(d_bk))
                {
                    d_bk = log(eps);
                }
                MATRIX_AT(score, b, k) = d_bk;
            }
        }

        // Responsibilities r_{bk}
        final_ll = 0.0;
        for (int b = 0; b < B; ++b)
        {
            for (int k = 0; k < mixture_h; ++k)
            {
                d_row[k] = MATRIX_AT(score, b, k);
            }

            double ll_b = log_sum_exp_weighted_phi(d_row, phi, mixture_h);
            final_ll += isfinite(ll_b) ? ll_b : log(eps);

            double den = 0.0;
            for (int k = 0; k < mixture_h; ++k)
            {
                double phi_k = (isfinite(phi[k]) && phi[k] > eps) ? phi[k] : eps;
                double e = 0.0;
                if (isfinite(ll_b) && isfinite(d_row[k]))
                {
                    e = exp(d_row[k] + log(phi_k) - ll_b); // ver esto
                }
                MATRIX_AT(responsibilities, b, k) = e;
                den += e;
            }
            if (!isfinite(den) || den <= 0)
            {
                fallback_underflow_count++;
                double phi_sum = 0.0;
                for (int k = 0; k < mixture_h; ++k)
                {
                    if (isfinite(phi[k]) && phi[k] > 0.0)
                        phi_sum += phi[k];
                }
                if (phi_sum > 0.0 && isfinite(phi_sum))
                {
                    fallback_phi_prior_count++;
                    for (int k = 0; k < mixture_h; ++k)
                    {
                        double v = (isfinite(phi[k]) && phi[k] > 0.0) ? (phi[k] / phi_sum) : 0.0;
                        MATRIX_AT(responsibilities, b, k) = v;
                    }
                }
                else
                {
                    fallback_uniform_count++;
                    den = (double)mixture_h;
                    for (int k = 0; k < mixture_h; ++k)
                    {
                        MATRIX_AT(responsibilities, b, k) = 1.0 / den;
                    }
                }
            }
            else
            {
                for (int k = 0; k < mixture_h; ++k)
                {
                    MATRIX_AT(responsibilities, b, k) /= den;
                }
            }
        }
        if (verbose && fallback_underflow_count > 0)
        {
            Rprintf("[debug] Iteration %d: responsibility fallback on %d ballots (phi-prior=%d, uniform=%d)\n", iter,
                    fallback_underflow_count, fallback_phi_prior_count, fallback_uniform_count);
        }

        // M-step for each component p_{kgc}
        double delta = 0.0;
        for (int k = 0; k < mixture_h; ++k)
        {
            EMContext *ctxk = components[k];
            for (int g = 0; g < G; ++g)
            {
                memset(row_counts, 0, (size_t)C * sizeof(double));
                double den = 0.0;
                for (int b = 0; b < B; ++b)
                {
                    double w = MATRIX_AT(responsibilities, b, k) * MATRIX_AT(ctxk->W, b, g);
                    if (w <= 0)
                    {
                        continue;
                    }
                    den += w;
                    for (int c = 0; c < C; ++c)
                    {
                        row_counts[c] += w * Q_3D(ctxk->q, b, g, c, G, C);
                    }
                }
                if (den <= eps)
                {
                    continue;
                }

                double row_sum = 0.0;
                for (int c = 0; c < C; ++c)
                {
                    double v = row_counts[c] / den;
                    if (!isfinite(v) || v < eps)
                    {
                        v = eps;
                    }
                    row_counts[c] = v;
                    row_sum += v;
                }
                if (!isfinite(row_sum) || row_sum <= 0.0)
                {
                    row_sum = (double)C;
                    for (int c = 0; c < C; ++c)
                    {
                        row_counts[c] = 1.0 / row_sum;
                    }
                }
                else
                {
                    for (int c = 0; c < C; ++c)
                    {
                        row_counts[c] /= row_sum;
                    }
                }

                for (int c = 0; c < C; ++c)
                {
                    double oldv = MATRIX_AT(ctxk->probabilities, g, c);
                    double newv = row_counts[c];
                    double d = fabs(newv - oldv);
                    if (d > delta)
                    {
                        delta = d;
                    }
                    MATRIX_AT(ctxk->probabilities, g, c) = newv;
                }
            }
        }

        // M-step for phi_k
        double phi_sum = 0.0;
        for (int k = 0; k < mixture_h; ++k)
        {
            double s = 0.0;
            for (int b = 0; b < B; ++b)
            {
                s += MATRIX_AT(responsibilities, b, k);
            }
            phi_new[k] = s / (double)B;
            if (!isfinite(phi_new[k]) || phi_new[k] < eps)
            {
                phi_new[k] = eps;
            }
            phi_sum += phi_new[k];
        }
        if (!isfinite(phi_sum) || phi_sum <= 0)
        {
            phi_sum = (double)mixture_h;
            for (int k = 0; k < mixture_h; ++k)
            {
                phi_new[k] = 1.0 / phi_sum;
            }
        }
        else
        {
            for (int k = 0; k < mixture_h; ++k)
            {
                phi_new[k] /= phi_sum;
            }
        }
        for (int k = 0; k < mixture_h; ++k)
        {
            double d = fabs(phi_new[k] - phi[k]);
            if (d > delta)
            {
                delta = d;
            }
            phi[k] = phi_new[k];
        }

        // Break label invariance with a deterministic component order.
        canonicalize_mixture_labels(components, phi, mixture_h);

        // Convergence metric in canonical order (not affected by label permutations).
        delta = 0.0;
        for (int k = 0; k < mixture_h; ++k)
        {
            for (int c = 0; c < C; ++c)
            {
                for (int g = 0; g < G; ++g)
                {
                    double oldv = prev_component_prob[(size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k)];
                    double newv = MATRIX_AT(components[k]->probabilities, g, c);
                    double d = fabs(newv - oldv);
                    if (d > delta)
                    {
                        delta = d;
                    }
                }
            }
            double dphi = fabs(phi[k] - prev_phi_params[k]);
            if (dphi > delta)
            {
                delta = dphi;
            }
        }

        if (verbose && iter % 1 == 0)
        {
            print_mixture_state(iter, final_ll, components, phi, mixture_h);
        }

        // Early stop safeguard: after 50 iterations, stop if log-likelihood decreases.
        if (iter > 10000000 && isfinite(prev_ll) && isfinite(final_ll) && final_ll < prev_ll)
        {
            if (verbose)
            {
                Rprintf("Stopping: log-likelihood decreased after iteration 50.\n");
            }
            finish_id = 0;
            break;
        }

        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC_RAW, &now);
        elapsed_total = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
        if (elapsed_total >= maxSeconds)
        {
            finish_id = 1;
            break;
        }
        if (iter >= inputParams->miniter && delta < convergence)
        {
            finish_id = 0;
            break;
        }

        prev_ll = final_ll;
    }

    // ---- Final E-step with updated parameters ---- //
    for (int k = 0; k < mixture_h; ++k)
    {
        EMContext *ctxk = components[k];
        double ll_dummy = 0.0;
        config.computeQ(ctxk, mix_params, &ll_dummy);
        if (inputParams->prob_cond_every)
        {
            if (strcmp(inputParams->prob_cond, "project_lp") == 0)
            {
                projectQ(ctxk, *inputParams);
            }
            else if (strcmp(inputParams->prob_cond, "lp") == 0)
            {
                for (int b = 0; b < B; ++b)
                {
                    LPW_ctx(ctxk, b);
                }
            }
        }

        for (int b = 0; b < B; ++b)
        {
            double d_bk = ctxk->ballot_loglik[b];
            if (!isfinite(d_bk))
            {
                d_bk = log(eps);
            }
            MATRIX_AT(score, b, k) = d_bk;
        }
    }
    final_ll = 0.0;
    for (int b = 0; b < B; ++b)
    {
        for (int k = 0; k < mixture_h; ++k)
        {
            d_row[k] = MATRIX_AT(score, b, k);
        }

        double ll_b = log_sum_exp_weighted_phi(d_row, phi, mixture_h);
        final_ll += isfinite(ll_b) ? ll_b : log(eps);

        double den = 0.0;
        for (int k = 0; k < mixture_h; ++k)
        {
            double phi_k = (isfinite(phi[k]) && phi[k] > eps) ? phi[k] : eps;
            double e = 0.0;
            if (isfinite(ll_b) && isfinite(d_row[k]))
            {
                e = exp(d_row[k] + log(phi_k) - ll_b);
            }
            MATRIX_AT(responsibilities, b, k) = e;
            den += e;
        }
        if (!isfinite(den) || den <= 0)
        {
            double phi_sum = 0.0;
            for (int k = 0; k < mixture_h; ++k)
            {
                if (isfinite(phi[k]) && phi[k] > 0.0)
                    phi_sum += phi[k];
            }
            if (phi_sum > 0.0 && isfinite(phi_sum))
            {
                for (int k = 0; k < mixture_h; ++k)
                {
                    double v = (isfinite(phi[k]) && phi[k] > 0.0) ? (phi[k] / phi_sum) : 0.0;
                    MATRIX_AT(responsibilities, b, k) = v;
                }
            }
            else
            {
                den = (double)mixture_h;
                for (int k = 0; k < mixture_h; ++k)
                {
                    MATRIX_AT(responsibilities, b, k) = 1.0 / den;
                }
            }
        }
        else
        {
            for (int k = 0; k < mixture_h; ++k)
            {
                MATRIX_AT(responsibilities, b, k) /= den;
            }
        }
    }

    // ---- Aggregate q and expected outcome ---- //
    size_t N = (size_t)B * (size_t)G * (size_t)C;
    double *q_mix = (double *)Calloc(N, double);
    double *expected_mix = (double *)Calloc(N, double);
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            double wbg = MATRIX_AT(components[0]->W, b, g);
            for (int c = 0; c < C; ++c)
            {
                double qv = 0.0;
                for (int k = 0; k < mixture_h; ++k)
                {
                    qv += MATRIX_AT(responsibilities, b, k) * Q_3D(components[k]->q, b, g, c, G, C);
                }
                Q_3D(q_mix, b, g, c, G, C) = qv;
                Q_3D(expected_mix, b, g, c, G, C) = wbg * qv;
            }
        }
    }

    Matrix final_prob = createMatrix(G, C);
    for (int g = 0; g < G; ++g)
    {
        double den = 0.0;
        for (int b = 0; b < B; ++b)
        {
            den += MATRIX_AT(components[0]->W, b, g);
        }
        if (!isfinite(den) || den <= 0)
        {
            den = 1.0;
        }
        for (int c = 0; c < C; ++c)
        {
            double num = 0.0;
            for (int b = 0; b < B; ++b)
            {
                num += Q_3D(expected_mix, b, g, c, G, C);
            }
            MATRIX_AT(final_prob, g, c) = num / den;
        }
    }

    double *component_prob = (double *)Calloc((size_t)G * (size_t)C * (size_t)mixture_h, double);
    for (int k = 0; k < mixture_h; ++k)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int g = 0; g < G; ++g)
            {
                component_prob[(size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k)] =
                    MATRIX_AT(components[k]->probabilities, g, c);
            }
        }
    }

    EMMixtureResult *res = Calloc(1, EMMixtureResult);
    res->probabilities = final_prob;
    res->q = q_mix;
    res->predicted_votes = expected_mix;
    res->phi = (double *)Calloc((size_t)mixture_h, double);
    memcpy(res->phi, phi, (size_t)mixture_h * sizeof(double));
    res->responsibilities = responsibilities;
    res->component_prob = component_prob;
    res->mixture_h = mixture_h;
    res->total_iterations = iter_done;
    res->total_time = elapsed_total;
    res->log_likelihood = final_ll;
    res->finish_id = finish_id;

    // ---- Cleanup temporaries ---- //
    freeMatrix(&score);
    Free(phi);
    Free(phi_new);
    Free(row_counts);
    Free(d_row);
    Free(prev_component_prob);
    Free(prev_phi_params);
    for (int k = 0; k < mixture_h; ++k)
    {
        cleanup(components[k]);
    }
    Free(components);

    return res;
}

static int pow_int_safe(int base, int exp, int *out)
{
    long double value = 1.0;
    for (int i = 0; i < exp; ++i)
    {
        value *= (long double)base;
        if (!(value <= (long double)INT_MAX))
        {
            return 0;
        }
    }
    *out = (int)value;
    return 1;
}

static void decode_assignment(int index, int H, int G, int *out)
{
    for (int g = 0; g < G; ++g)
    {
        out[g] = index % H;
        index /= H;
    }
}

static inline size_t comp_idx(int g, int c, int k, int G, int C)
{
    return (size_t)g + (size_t)G * ((size_t)c + (size_t)C * (size_t)k);
}

static inline double log_sum_exp_terms(const double *terms, int n)
{
    double max_term = -INFINITY;
    for (int i = 0; i < n; ++i)
    {
        if (isfinite(terms[i]) && terms[i] > max_term)
        {
            max_term = terms[i];
        }
    }
    if (!isfinite(max_term))
    {
        return -INFINITY;
    }
    double sum_exp = 0.0;
    for (int i = 0; i < n; ++i)
    {
        if (isfinite(terms[i]))
        {
            sum_exp += exp(terms[i] - max_term);
        }
    }
    if (!(sum_exp > 0.0) || !isfinite(sum_exp))
    {
        return -INFINITY;
    }
    return max_term + log(sum_exp);
}

static void run_row_mixture_estep(EMContext *eval_ctx, const double *component_prob, const int *assignments, int S,
                                  int H, int B, int G, int C, size_t N, QMethodConfig config, QMethodInput ll_params,
                                  QMethodInput *inputParams, Matrix *score, double *q_store, double eps)
{
    (void)H;
    // E-step for row-level mixture:
    // for each assignment h in {1,...,H}^G, build P_h(g,c)=p_{h_g,g,c},
    // then compute q_{b,h,g,c} and d_{b,h}=log P(X_b | h; p).
    for (int s = 0; s < S; ++s)
    {
        const int *assign = assignments + (size_t)s * (size_t)G;
        for (int g = 0; g < G; ++g)
        {
            const int kg = assign[g];
            for (int c = 0; c < C; ++c)
            {
                MATRIX_AT(eval_ctx->probabilities, g, c) = component_prob[comp_idx(g, c, kg, G, C)];
            }
        }

        double ll_dummy = 0.0;
        config.computeQ(eval_ctx, ll_params, &ll_dummy);
        if (inputParams->prob_cond_every)
        {
            if (strcmp(inputParams->prob_cond, "project_lp") == 0)
            {
                projectQ(eval_ctx, *inputParams);
            }
            else if (strcmp(inputParams->prob_cond, "lp") == 0)
            {
                for (int b = 0; b < B; ++b)
                {
                    LPW_ctx(eval_ctx, b);
                }
            }
        }

        for (int b = 0; b < B; ++b)
        {
            double d_bh = eval_ctx->ballot_loglik[b];
            if (!isfinite(d_bh))
            {
                d_bh = log(eps);
            }
            MATRIX_AT_PTR(score, b, s) = d_bh;
        }

        memcpy(q_store + (size_t)s * N, eval_ctx->q, N * sizeof(double));
    }
}

static double update_row_mixture_responsibilities(const Matrix *score, Matrix *responsibilities, const double *phi_gk,
                                                  const int *assignments, int B, int G, int H, int S, double *terms,
                                                  double *log_prior, double eps)
{
    // log_prior[h] = log P(h | phi) = sum_g log(phi_{g,h_g}).
    for (int s = 0; s < S; ++s)
    {
        const int *assign = assignments + (size_t)s * (size_t)G;
        double lp = 0.0;
        for (int g = 0; g < G; ++g)
        {
            double v = phi_gk[(size_t)g * (size_t)H + (size_t)assign[g]];
            if (!isfinite(v) || v < eps)
            {
                v = eps;
            }
            lp += log(v);
        }
        log_prior[s] = lp;
    }

    double final_ll = 0.0;
    for (int b = 0; b < B; ++b)
    {
        for (int s = 0; s < S; ++s)
        {
            // terms[h] = d_{b,h} + log P(h | phi).
            terms[s] = MATRIX_AT_PTR(score, b, s) + log_prior[s];
        }

        // ll_b = log(sum_h exp(d_{b,h}) * P(h | phi)).
        double ll_b = log_sum_exp_terms(terms, S);
        final_ll += isfinite(ll_b) ? ll_b : log(eps);

        double den = 0.0;
        for (int s = 0; s < S; ++s)
        {
            double e = 0.0;
            if (isfinite(ll_b) && isfinite(terms[s]))
            {
                e = exp(terms[s] - ll_b);
            }
            MATRIX_AT_PTR(responsibilities, b, s) = e;
            den += e;
        }

        if (!isfinite(den) || den <= 0.0)
        {
            const double uniform = 1.0 / (double)S;
            for (int s = 0; s < S; ++s)
            {
                MATRIX_AT_PTR(responsibilities, b, s) = uniform;
            }
        }
        else
        {
            for (int s = 0; s < S; ++s)
            {
                MATRIX_AT_PTR(responsibilities, b, s) /= den;
            }
        }
        // responsibilities[b,h] = r_{b,h} = P(h | X_b; tau).
    }
    return final_ll;
}

static void print_row_mixture_state(int iter, double loglik, const double *component_prob, const double *phi, int G,
                                    int C, int H)
{
    Rprintf("\n----------\nIteration: %d\nLog-likelihood: %.10f\n", iter, loglik);
    for (int k = 0; k < H; ++k)
    {
        Matrix Pk = {(double *)(component_prob + (size_t)k * (size_t)G * (size_t)C), G, C};
        Rprintf("Component %d probability matrix (GxC):\n", k + 1);
        printMatrix(&Pk);
    }

    Rprintf("Phi matrix (GxH):\nMatrix (%dx%d):\n", G, H);
    for (int g = 0; g < G; ++g)
    {
        Rprintf("| ");
        for (int k = 0; k < H; ++k)
        {
            double v = phi[(size_t)g * (size_t)H + (size_t)k];
            if (k < H - 1)
            {
                Rprintf("%4.3f\t", v);
            }
            else
            {
                Rprintf("%4.3f", v);
            }
        }
        Rprintf(" |\n");
    }
}

EMMixtureResult *EMAlgoritmRowMixture(Matrix *X, Matrix *W, const char *p_method, const char *q_method,
                                      const double convergence, const double LLconvergence, const int maxIter,
                                      const double maxSeconds, const bool verbose, Matrix *probMatrix,
                                      QMethodInput *inputParams, int row_mixture_h)
{
    (void)LLconvergence;

    const double eps = 1e-12;
    const int H = row_mixture_h;
    if (H < 1)
    {
        error("run_em: Invalid 'H'. Must be a positive integer.");
    }
    if (H == 1)
    {
        return EMAlgoritmMixture(X, W, p_method, q_method, convergence, LLconvergence, maxIter, maxSeconds, verbose,
                                 probMatrix, inputParams, 1);
    }

    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double elapsed_total = 0.0;

    QMethodInput ctx_params = *inputParams;
    ctx_params.computeLL = true;
    EMContext *eval_ctx = createEMContext(X, W, q_method, ctx_params);
    getInitialP(eval_ctx, p_method, probMatrix);
    normalizeProbabilityRows(&eval_ctx->probabilities);

    const int B = (int)eval_ctx->B;
    const int G = (int)eval_ctx->G;
    const int C = (int)eval_ctx->C;

    int S = 0;
    if (!pow_int_safe(H, G, &S) || S <= 0)
    {
        cleanup(eval_ctx);
        error("run_em: H^G is too large to handle.");
    }

    int *assignments = (int *)Calloc((size_t)S * (size_t)G, int);
    for (int s = 0; s < S; ++s)
    {
        // assignments[s,g] stores h_g for the s-th assignment h.
        decode_assignment(s, H, G, assignments + (size_t)s * (size_t)G);
    }

    double *component_prob = (double *)Calloc((size_t)G * (size_t)C * (size_t)H, double);
    for (int k = 0; k < H; ++k)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int g = 0; g < G; ++g)
            {
                component_prob[comp_idx(g, c, k, G, C)] = MATRIX_AT(eval_ctx->probabilities, g, c);
            }
        }
    }
    GetRNGstate();
    for (int k = 1; k < H; ++k)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int g = 0; g < G; ++g)
            {
                double jit = 0.9 + 0.2 * unif_rand();
                component_prob[comp_idx(g, c, k, G, C)] *= jit;
            }
        }
        for (int g = 0; g < G; ++g)
        {
            double row_sum = 0.0;
            for (int c = 0; c < C; ++c)
            {
                double v = component_prob[comp_idx(g, c, k, G, C)];
                if (!isfinite(v) || v < eps)
                {
                    v = eps;
                }
                component_prob[comp_idx(g, c, k, G, C)] = v;
                row_sum += v;
            }
            if (!(row_sum > 0.0) || !isfinite(row_sum))
            {
                row_sum = (double)C;
                for (int c = 0; c < C; ++c)
                {
                    component_prob[comp_idx(g, c, k, G, C)] = 1.0 / row_sum;
                }
            }
            else
            {
                for (int c = 0; c < C; ++c)
                {
                    component_prob[comp_idx(g, c, k, G, C)] /= row_sum;
                }
            }
        }
    }
    PutRNGstate();

    double *phi = (double *)Calloc((size_t)G * (size_t)H, double);
    double *phi_new = (double *)Calloc((size_t)G * (size_t)H, double);
    for (int g = 0; g < G; ++g)
    {
        for (int k = 0; k < H; ++k)
        {
            phi[(size_t)g * (size_t)H + (size_t)k] = 1.0 / (double)H;
        }
    }

    Matrix score = createMatrix(B, S);
    Matrix responsibilities = createMatrix(B, S);
    size_t N = (size_t)B * (size_t)G * (size_t)C;
    double *q_store = (double *)Calloc((size_t)S * N, double);
    double *terms = (double *)Calloc((size_t)S, double);
    double *log_prior = (double *)Calloc((size_t)S, double);
    double *acc_num = (double *)Calloc((size_t)H * (size_t)G * (size_t)C, double);
    double *acc_den = (double *)Calloc((size_t)H * (size_t)G, double);

    QMethodConfig config = getQMethodConfig(q_method, *inputParams);
    QMethodInput mix_params = config.params;
    mix_params.computeLL = true;

    int finish_id = 2;
    int iter_done = 0;
    double final_ll = NA_REAL;
    double prev_ll = -INFINITY;

    for (int iter = 1; iter <= maxIter; ++iter)
    {
        iter_done = iter;

        // E-step: compute q_{b,h,g,c} and d_{b,h} for all assignments h.
        run_row_mixture_estep(eval_ctx, component_prob, assignments, S, H, B, G, C, N, config, mix_params, inputParams,
                              &score, q_store, eps);

        // E-step posterior over assignments: r_{b,h} = P(h | X_b; tau).
        final_ll = update_row_mixture_responsibilities(&score, &responsibilities, phi, assignments, B, G, H, S, terms,
                                                       log_prior, eps);

        memset(acc_num, 0, (size_t)H * (size_t)G * (size_t)C * sizeof(double));
        memset(acc_den, 0, (size_t)H * (size_t)G * sizeof(double));

        for (int s = 0; s < S; ++s)
        {
            const int *assign = assignments + (size_t)s * (size_t)G;
            const double *q_s = q_store + (size_t)s * N;
            for (int b = 0; b < B; ++b)
            {
                double r = MATRIX_AT(responsibilities, b, s);
                if (!(r > 0.0) || !isfinite(r))
                {
                    continue;
                }
                for (int g = 0; g < G; ++g)
                {
                    int k = assign[g];
                    double w = r * MATRIX_AT(eval_ctx->W, b, g);
                    if (!(w > 0.0) || !isfinite(w))
                    {
                        continue;
                    }
                    acc_den[(size_t)g + (size_t)G * (size_t)k] += w;
                    for (int c = 0; c < C; ++c)
                    {
                        double qv = q_s[(size_t)b * (size_t)G * (size_t)C + (size_t)c * (size_t)G + (size_t)g];
                        // Accumulates numerator for:
                        // p_{k,g,c} <- sum_b sum_{h in H_{gk}} r_{b,h} w_{b,g} q_{b,h,g,c}.
                        acc_num[comp_idx(g, c, k, G, C)] += w * qv;
                    }
                }
            }
        }

        double delta = 0.0;
        for (int k = 0; k < H; ++k)
        {
            for (int g = 0; g < G; ++g)
            {
                double row_sum = 0.0;
                for (int c = 0; c < C; ++c)
                {
                    row_sum += acc_num[comp_idx(g, c, k, G, C)];
                }
                if (!(row_sum > eps) || !isfinite(row_sum))
                {
                    continue;
                }
                for (int c = 0; c < C; ++c)
                {
                    double oldv = component_prob[comp_idx(g, c, k, G, C)];
                    // M-step closed form:
                    // p_{k,g,c} = num_{k,g,c} / sum_{c'} num_{k,g,c'}.
                    double newv = acc_num[comp_idx(g, c, k, G, C)] / row_sum;
                    double d = fabs(newv - oldv);
                    if (d > delta)
                    {
                        delta = d;
                    }
                    component_prob[comp_idx(g, c, k, G, C)] = newv;
                }
            }
        }

        memset(phi_new, 0, (size_t)G * (size_t)H * sizeof(double));
        for (int b = 0; b < B; ++b)
        {
            for (int s = 0; s < S; ++s)
            {
                double r = MATRIX_AT(responsibilities, b, s);
                if (!(r > 0.0) || !isfinite(r))
                {
                    continue;
                }
                const int *assign = assignments + (size_t)s * (size_t)G;
                for (int g = 0; g < G; ++g)
                {
                    // Accumulates sum_b sum_{h in H_{gk}} r_{b,h} for phi_{g,k}.
                    phi_new[(size_t)g * (size_t)H + (size_t)assign[g]] += r;
                }
            }
        }
        for (int g = 0; g < G; ++g)
        {
            double row_sum = 0.0;
            for (int k = 0; k < H; ++k)
            {
                // For fixed g, denominator is:
                // sum_{k'} sum_b sum_{h in H_{gk'}} r_{b,h} = B.
                double v = phi_new[(size_t)g * (size_t)H + (size_t)k] / (double)B;
                if (!isfinite(v) || v < eps)
                {
                    v = eps;
                }
                phi_new[(size_t)g * (size_t)H + (size_t)k] = v;
                row_sum += v;
            }
            if (!(row_sum > 0.0) || !isfinite(row_sum))
            {
                row_sum = (double)H;
                for (int k = 0; k < H; ++k)
                {
                    phi_new[(size_t)g * (size_t)H + (size_t)k] = 1.0 / row_sum;
                }
            }
            else
            {
                for (int k = 0; k < H; ++k)
                {
                    phi_new[(size_t)g * (size_t)H + (size_t)k] /= row_sum;
                }
            }
        }
        for (int g = 0; g < G; ++g)
        {
            for (int k = 0; k < H; ++k)
            {
                size_t idx = (size_t)g * (size_t)H + (size_t)k;
                double d = fabs(phi_new[idx] - phi[idx]);
                if (d > delta)
                {
                    delta = d;
                }
                phi[idx] = phi_new[idx];
            }
        }

        if (verbose && iter % 1 == 0)
        {
            print_row_mixture_state(iter, final_ll, component_prob, phi, G, C, H);
        }

        // Early stop safeguard: after 50 iterations, stop if log-likelihood decreases.
        // if (iter > 50 && isfinite(prev_ll) && isfinite(final_ll) && final_ll < prev_ll)
        // {
        //     if (verbose)
        //     {
        //         Rprintf("Stopping: log-likelihood decreased after iteration 50.\n");
        //     }
        //     finish_id = 0;
        //     break;
        // }

        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC_RAW, &now);
        elapsed_total = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
        if (elapsed_total >= maxSeconds)
        {
            finish_id = 1;
            break;
        }
        if (iter >= inputParams->miniter && delta < convergence)
        {
            finish_id = 0;
            break;
        }

        prev_ll = final_ll;
    }

    run_row_mixture_estep(eval_ctx, component_prob, assignments, S, H, B, G, C, N, config, mix_params, inputParams,
                          &score, q_store, eps);
    final_ll = update_row_mixture_responsibilities(&score, &responsibilities, phi, assignments, B, G, H, S, terms,
                                                   log_prior, eps);

    double *q_mix = (double *)Calloc(N, double);
    double *expected_mix = (double *)Calloc(N, double);
    for (int b = 0; b < B; ++b)
    {
        for (int g = 0; g < G; ++g)
        {
            double wbg = MATRIX_AT(eval_ctx->W, b, g);
            for (int c = 0; c < C; ++c)
            {
                double qv = 0.0;
                for (int s = 0; s < S; ++s)
                {
                    const double *q_s = q_store + (size_t)s * N;
                    // Final conditional mean:
                    // q_{b,g,c} = sum_h r_{b,h} q_{b,h,g,c}.
                    qv += MATRIX_AT(responsibilities, b, s) *
                          q_s[(size_t)b * (size_t)G * (size_t)C + (size_t)c * (size_t)G + (size_t)g];
                }
                Q_3D(q_mix, b, g, c, G, C) = qv;
                Q_3D(expected_mix, b, g, c, G, C) = wbg * qv;
            }
        }
    }

    Matrix final_prob = createMatrix(G, C);
    for (int g = 0; g < G; ++g)
    {
        double den = 0.0;
        for (int b = 0; b < B; ++b)
        {
            den += MATRIX_AT(eval_ctx->W, b, g);
        }
        if (!(den > 0.0) || !isfinite(den))
        {
            den = 1.0;
        }
        for (int c = 0; c < C; ++c)
        {
            double num = 0.0;
            for (int b = 0; b < B; ++b)
            {
                num += Q_3D(expected_mix, b, g, c, G, C);
            }
            MATRIX_AT(final_prob, g, c) = num / den;
        }
    }

    EMMixtureResult *res = Calloc(1, EMMixtureResult);
    res->probabilities = final_prob;
    res->q = q_mix;
    res->predicted_votes = expected_mix;
    res->phi = (double *)Calloc((size_t)G * (size_t)H, double);
    memcpy(res->phi, phi, (size_t)G * (size_t)H * sizeof(double));
    res->responsibilities = responsibilities;
    res->component_prob = component_prob;
    res->mixture_h = H;
    res->total_iterations = iter_done;
    res->total_time = elapsed_total;
    res->log_likelihood = final_ll;
    res->finish_id = finish_id;

    freeMatrix(&score);
    Free(phi);
    Free(phi_new);
    Free(q_store);
    Free(terms);
    Free(log_prior);
    Free(acc_num);
    Free(acc_den);
    Free(assignments);
    cleanup(eval_ctx);

    return res;
}

void cleanupMixtureResult(EMMixtureResult *res)
{
    if (res == NULL)
    {
        return;
    }
    if (res->probabilities.data != NULL)
    {
        freeMatrix(&res->probabilities);
        res->probabilities.data = NULL;
    }
    if (res->q != NULL)
    {
        Free(res->q);
        res->q = NULL;
    }
    if (res->predicted_votes != NULL)
    {
        Free(res->predicted_votes);
        res->predicted_votes = NULL;
    }
    if (res->phi != NULL)
    {
        Free(res->phi);
        res->phi = NULL;
    }
    if (res->responsibilities.data != NULL)
    {
        freeMatrix(&res->responsibilities);
        res->responsibilities.data = NULL;
    }
    if (res->component_prob != NULL)
    {
        Free(res->component_prob);
        res->component_prob = NULL;
    }
    Free(res);
}
