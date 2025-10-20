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

#include "main.h"
#include "globals.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <Rinternals.h>
#include <dirent.h>
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

// ---- Inititalize global variables ---- //
uint32_t TOTAL_VOTES = 0;
uint32_t TOTAL_BALLOTS = 0;
uint16_t TOTAL_CANDIDATES = 0;
uint16_t TOTAL_GROUPS = 0;
// ---...--- //

EMContext *createEMContext(Matrix *X, Matrix *W, const char *method, QMethodInput params)
{
    // Create the context object
    EMContext *ctx = Calloc(1, EMContext);

    // Generate the C, B, G parameters
    ctx->C = X->rows;
    ctx->B = X->cols;
    ctx->G = W->cols;
    TOTAL_BALLOTS = ctx->B;
    TOTAL_CANDIDATES = ctx->C;
    TOTAL_GROUPS = ctx->G;

    // Generate the matrices
    ctx->X = *copMatrixPtr(X);
    ctx->W = *copMatrixPtr(W);
    ctx->intX = copMatrixDI(X);
    ctx->intW = copMatrixDI(W);
    ctx->probabilities = createMatrix(ctx->G, ctx->C);
    ctx->q = (double *)Calloc(ctx->B * ctx->C * ctx->G, double);
    ctx->iteration = 0;

    ctx->ballots_votes = Calloc(ctx->B, uint16_t);
    ctx->inv_ballots_votes = Calloc(ctx->B, double);
    ctx->candidates_votes = Calloc(ctx->C, uint32_t);
    ctx->group_votes = Calloc(ctx->G, uint32_t);
    ctx->total_votes = 0;

    ctx->Wnorm = precomputeNorm(W);

    // Fill utility arrays
    for (uint32_t b = 0; b < ctx->B; ++b)
    {
        for (uint16_t c = 0; c < ctx->C; ++c)
        {
            uint32_t v = MATRIX_AT(ctx->intX, c, b);
            ctx->candidates_votes[c] += v;
            ctx->ballots_votes[b] += v;
            ctx->total_votes += v;
        }
        ctx->inv_ballots_votes[b] = 1.0 / (double)ctx->ballots_votes[b];
        for (uint16_t g = 0; g < ctx->G; ++g)
            ctx->group_votes[g] += MATRIX_AT(ctx->intW, b, g);
    }
    TOTAL_VOTES = ctx->total_votes;

    if (strcmp(method, "mult") == 0 && params.computeLL)
    {
        precomputeLogGammas(ctx);
    }
    if (strcmp(method, "mcmc") == 0)
    {
        precomputeLogGammas(ctx);
        generateOmegaSet(ctx, params.M, params.S, params.burnInSteps);
        encode(ctx);
        precomputeQConstant(ctx, params.S);
        preComputeMultinomial(ctx);
    }
    if (strcmp(method, "exact") == 0)
    {
        generateKSets(ctx);
        generateHSets(ctx);
    }
    if (strcmp(method, "mvn_cdf") == 0)
    {
        allocateSeed(ctx, params);
    }

    return ctx;
}

Matrix precomputeNorm(Matrix *W)
{
    Matrix returnMat = createMatrix(TOTAL_BALLOTS, TOTAL_GROUPS);
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        double sum = 0;
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            sum += pow(MATRIX_AT_PTR(W, b, g), 2);
        }
        // sum = sqrt(sum);
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            MATRIX_AT(returnMat, b, g) = MATRIX_AT_PTR(W, b, g) / sum;
        }
    }
    return returnMat;
}

/*
 * @brief Computes the predicted votes outcome for each ballot box
 */
void getPredictedVotes(EMContext *ctx)
{
    ctx->predicted_votes = (double *)Calloc(ctx->B * ctx->C * ctx->G, double);
    double *q = ctx->q;

    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            int W_bg = MATRIX_AT(ctx->intW, b, g);
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {
                Q_3D(ctx->predicted_votes, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) +=
                    W_bg * Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
            }
        }
    }
}

/**
 * @brief Computes the initial probability of the EM algoritm.
 *
 * Given the observables results, it computes a convenient initial "p" value for initiating the
 * algorithm. Currently it supports the "uniform", "group_proportional" and "proportional" methods.
 *
 * @param[in] p_method The method for calculating the initial parameter. Currently it supports "uniform",
 * "group_proportional" and "proportional" methods.
 *
 * @return Matrix of dimension (gxc) with the initial probability for each demographic group "g" voting for a given
 * candidate "c".
 * @note This should be used only that the first iteration of the EM-algorithm.
 * @warning
 * - Pointers shouldn't be NULL.
 * - `x` and `w` dimensions must be coherent.
 *
 */
void getInitialP(EMContext *ctx, const char *p_method, Matrix *probMatrix)
{

    // bool custom_matrix = MATRIX_AT_PTR(probMatrix, 0, 0) != -1 ? true : false;
    // ---- Validation: check the method input ----//
    if (strcmp(p_method, "custom") && strcmp(p_method, "uniform") != 0 && strcmp(p_method, "proportional") != 0 &&
        strcmp(p_method, "group_proportional") != 0 && strcmp(p_method, "random") != 0 &&
        strcmp(p_method, "mult") != 0 && strcmp(p_method, "mvn_cdf") != 0 && strcmp(p_method, "mvn_pdf") != 0 &&
        strcmp(p_method, "exact") != 0)
    {
        error("run_em: The method `%s` to calculate the initial probability doesn't exist.\nThe supported methods "
              "are: `uniform`, `proportional`, `random`, `group_proportional`, `mult`, `mvn_cdf`, `mvn_pdf` and "
              "`exact`.\n",
              p_method);
    }
    // ---...--- //
    // ---- Compute the random method ---- //
    else if (strcmp(p_method, "random") == 0)
    {
        // Integrate with R's RNG
        GetRNGstate();

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            double rowSum = 0.0;

            // Fill row with random values in [0,1)
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double r = unif_rand(); // R's uniform RNG
                MATRIX_AT(ctx->probabilities, g, c) = r;
                rowSum += r;
            }

            // Normalize row so that it sums to 1
            if (rowSum > 0.0)
            {
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    MATRIX_AT(ctx->probabilities, g, c) /= rowSum;
                }
            }
            else
            {
                // In the unlikely event rowSum is exactly 0.0, assign uniform
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    MATRIX_AT(ctx->probabilities, g, c) = 1.0 / (double)TOTAL_CANDIDATES;
                }
            }
        }
        PutRNGstate();
    }

    // ---- Compute the uniform method ---- //
    // ---- It assumes a uniform distribution among candidates ----
    else if (strcmp(p_method, "uniform") == 0)
    {
        fillMatrix(&ctx->probabilities, 1.0 / (double)TOTAL_CANDIDATES);
    }
    // ---...--- //

    // ---- Compute the proportional method ---- //
    // ---- It calculates the proportion of votes of each candidate, and assigns that probability to every demographic
    // group ----
    else if (strcmp(p_method, "proportional") == 0)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate
            double ratio = (double)ctx->candidates_votes[c] /
                           (double)TOTAL_VOTES; // Proportion of candidates votes per total votes.
            for (int g = 0; g < TOTAL_GROUPS; g++)
            { // --- For each group, given a candidate
                MATRIX_AT(ctx->probabilities, g, c) = ratio;
            } // --- End group loop
        } // --- End candidate loop
    }
    // ---...--- //

    // ---- Compute the group_proportional method ---- //
    // ---- Considers the proportion of candidates votes and demographic groups aswell ----
    else if (strcmp(p_method, "group_proportional") == 0)
    {
        // ---- Create a temporary matrix to store the first results ----
        Matrix ballotProbability = createMatrix(TOTAL_BALLOTS, TOTAL_CANDIDATES);
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot vote
            int den = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate, given a ballot box
                MATRIX_AT(ballotProbability, b, c) = MATRIX_AT(ctx->intX, c, b);
                den += MATRIX_AT(ballotProbability, b, c);
            }
            // ---- Handle border case ----
            if (den != 0)
            {
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // --- For each candidate, given a ballot box
                    MATRIX_AT(ballotProbability, b, c) /= (double)den;
                }
            }
        } // --- End ballot box loop

        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot box
            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            { // --- For each group given a ballot box
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // --- For each candidate, given a ballot box and a group
                    MATRIX_AT(ctx->probabilities, g, c) +=
                        MATRIX_AT(ballotProbability, b, c) * MATRIX_AT(ctx->intW, b, g);
                }
            }
        }

        // ---- Add the final values to the matrix
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group
                // ---- Handle border case ----
                if (ctx->group_votes[g] == 0)
                    MATRIX_AT(ctx->probabilities, g, c) = 0;
                else
                    MATRIX_AT(ctx->probabilities, g, c) /= (double)ctx->group_votes[g];
            }
        }
        freeMatrix(&ballotProbability);
    }
    else if (strcmp(p_method, "custom") == 0)
    {
        if (probMatrix->rows != ctx->G || probMatrix->cols != ctx->C)
        {
            error("Custom matrix dimensions do not match the expected size of (%d x %d).", ctx->G, ctx->C);
        }
        ctx->probabilities = *copMatrixPtr(probMatrix);
    }
    else
    {
        int iterTotal, finishing_reason;
        double time, logLLarr;
        QMethodInput inputParams = {0};
        EMContext *newCtx =
            EMAlgoritm(&ctx->X, &ctx->W, "group_proportional", p_method, 0.001, 0.0001, 1000, 1000, false, &time,
                       &iterTotal, &logLLarr, &finishing_reason, probMatrix, &inputParams);
        ctx->probabilities = createMatrix(newCtx->probabilities.rows, newCtx->probabilities.cols);
        ctx->q = newCtx->q;
        ctx->predicted_votes = newCtx->predicted_votes;
        // printMatrix(&newCtx->probabilities);
        size_t nel = (size_t)ctx->probabilities.rows * ctx->probabilities.cols;

        memcpy(ctx->probabilities.data, newCtx->probabilities.data, nel * sizeof *ctx->probabilities.data);
        // now compute number of elements in qMultinomial
        size_t nel2 = nel * (size_t)newCtx->W.rows;

        // allocate storage for q
        ctx->q = (double *)Calloc(nel2, double);
        ctx->predicted_votes = (double *)Calloc(nel2, double);

        if (!ctx->q)
        {
            error("Allocation error, submit a ticket in the Github repository.");
        }

        // copy qMultinomial
        memcpy(ctx->q, newCtx->q, nel2 * sizeof *ctx->q);
        memcpy(ctx->predicted_votes, newCtx->predicted_votes, nel2 * sizeof *ctx->predicted_votes);

        // cleanup(newCtx);
    }
    // ---...--- //
}

/*
 *
 * @brief Sets the configuration for the Q method, using a modularized approach.
 *
 * Given that different `Q` methods receive different parameters, a modularized approach is given towards each method
 *
 * @input[in] q_method A char with the q_method. Currently it supports "exact", "mcmc", "mult", "mvn_cdf", "metropolis",
 * and "mvn_pdf"
 * @input[in] inputParams A QMethodInput struct, that should be defined in a main function, with the parameters for the
 * distinct methods
 *
 * return A QMethodConfig struct that defines a function pointer towards the corresponding process for getting the `Q`
 * parameter according the method given.
 */
QMethodConfig getQMethodConfig(const char *q_method, QMethodInput inputParams)
{
    QMethodConfig config = {NULL}; // Initialize everything to NULL/0

    config.computeQ = computeQMultinomial;

    if (strcmp(q_method, "mult") == 0)
    {
        config.computeQ = computeQMultinomial;
    }
    else if (strcmp(q_method, "mcmc") == 0)
    {
        config.computeQ = computeQHitAndRun;
    }
    else if (strcmp(q_method, "mvn_pdf") == 0)
    {
        config.computeQ = computeQMultivariatePDF;
    }
    else if (strcmp(q_method, "exact") == 0)
    {
        config.computeQ = computeQExact;
    }
    else if (strcmp(q_method, "mvn_cdf") == 0)
    {
        config.computeQ = computeQMultivariateCDF;
    }
    else
    {
        error("Compute: An invalid method was provided: `%s`\nThe supported methods are: `exact`, `mcmc`"
              ", `mult`, `mvn_cdf` and `mvn_pdf`.\n",
              q_method);
    }

    // Directly store the input parameters
    config.params = inputParams;
    return config;
}

/*
 * @brief Computes the optimal solution for the `M` step
 *
 * Given the conditional probability and the votations per demographic group, it calculates the new probability for
 * the next iteration.
 *
 * @param[in] q Array of matrices of dimension (bxgxc) that represents the probability that a voter of group "g" in
 * ballot box "b" voted for candidate "c" conditional on the observed result.
 *
 * @return A matrix with the optimal probabilities according maximizing the Log-likelihood.
 *
 * @see getInitialP() for getting initial probabilities. This method is recommended to be used exclusively for the EM
 * Algorithm, unless there's a starting "q" to start with.
 *
 */
void getP(EMContext *ctx)
{
    // ---- Inititalize variables ---- //
    const double *q = ctx->q;

    // ---- Compute the dot products ---- //
    int stride = TOTAL_GROUPS * TOTAL_CANDIDATES;
    int tBal = TOTAL_BALLOTS;
    int newStride = 1;
    for (int g = 0; g < TOTAL_GROUPS; g++)
    { // --- For each group
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate given a group
            // Dot product over b=0..B-1 of W_{b,g} * Q_{b,g,c}
            const double *baseY = q + (c * TOTAL_GROUPS + g);

            double val;
            val = F77_CALL(ddot)(&tBal,
                                 &ctx->W.data[g * TOTAL_BALLOTS], // indexing W in column-major
                                 &newStride,                      // Column-major: stride is 1 for W
                                 baseY,                           // Column-major: index properly
                                 &stride                          // Stride: move down rows (1 step per row)
            );

            MATRIX_AT(ctx->probabilities, g, c) = val / ctx->group_votes[g];
        }
    }
    // ---...--- //
}

bool hasMismatch(EMContext *ctx)
{
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        double sumW = 0.0;
        double sumX = 0.0;
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            sumW += MATRIX_AT(ctx->W, b, g);
        }
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            sumX += MATRIX_AT(ctx->X, c, b);
        }
        if (fabs(sumW - sumX) > 1e-6)
        {
            return true;
        }
    }
    return false;
}

void projectQ(EMContext *ctx, QMethodInput inputParams)
{
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    Matrix *norm = &ctx->Wnorm;
    bool mismatch = hasMismatch(ctx);
    // getPredictedVotes(ctx); // Obtain WQ

    Matrix temp = createMatrix(TOTAL_BALLOTS, TOTAL_CANDIDATES);

    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double sum = 0.0;
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                // Here we should rescale if there's a mismatch between W and X totals
                if (!mismatch)
                {
                    sum += Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) * MATRIX_AT(ctx->W, b, g);
                }
                else
                {
                    double sum_x = ctx->ballots_votes[b];
                    double sum_w = 0.0;
                    for (int g = 0; g < TOTAL_GROUPS; g++)
                        sum_w += MATRIX_AT(ctx->W, b, g);
                    sum +=
                        Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) * MATRIX_AT(ctx->W, b, g) * sum_x / sum_w;
                }
            }
            MATRIX_AT(temp, b, c) = sum;
        }
    }
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {

        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double predictedVote = MATRIX_AT(temp, b, c);
                Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) =
                    Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) -
                    (predictedVote - MATRIX_AT_PTR(X, c, b)) * MATRIX_AT_PTR(norm, b, g);
            }
        }
    }

    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {
                if (Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) < 0 ||
                    Q_3D(ctx->q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) > 1)
                {
                    int status;
                    status = LPW(ctx, b);
                }
            }
        }
    }
}

int checkGroups(EMContext ctx)
{
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        if (ctx.group_votes[g] == 0)
        {
            ctx.total_votes -= ctx.group_votes[g];
            TOTAL_GROUPS--;
            return g;
        }
    }
    return -1;
}

/**
 * @brief Implements the whole EM algorithm.
 *
 * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
 * paper, it currently supports mcmc, mult, mvn_cdf and mvn_pdf methods.
 *
 * @param[in, out] X The candidate matrix, with the votes per candidate per ballot box.
 * @param[in, out] W The demographic matrix, with the votes per demographic group per ballot box.
 * @param[in] p_method Pointer to a string that indicates the method or calculating "p".
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
 * @param[in] LLconvergence Threshold for the log-likelihood convergence regarding its variation.
 * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
 * @param[in] verbose Wether to verbose useful outputs.
 *
 * @return Matrix: A matrix with the final probabilities. In case it doesn't converges, it returns the last
 * probability that was computed
 *
 * @note This is the main function that calls every other function for "q"
 *
 * @see getInitialP() for getting initial probabilities. group_proportional method is recommended.
 *
 * @warning
 * - Pointers shouldn't be NULL.
 * - `x` and `w` dimensions must be coherent.
 *
 */
EMContext *EMAlgoritm(Matrix *X, Matrix *W, const char *p_method, const char *q_method, const double convergence,
                      const double LLconvergence, const int maxIter, const double maxSeconds, const bool verbose,
                      double *time, int *iterTotal, double *logLLarr, int *finishing_reason, Matrix *probMatrix,
                      QMethodInput *inputParams)
{
    // ---- Error handling is done on getQMethodConfig!
    if (verbose)
    {
        Rprintf("Starting the EM algorithm.\n");
        Rprintf("Conditional probability will be estimated using the '%s' method with the following "
                "parameters:\n- Probability convergence threshold:\t%f\n- Log-likelihood convergence "
                "threshold:\t%f\n- Maximum number of iterations:\t%d\n",
                q_method, convergence, LLconvergence, maxIter);
    }

    // ---- Define the parameters for the main loop ---- //
    // ---- Start timer
    struct timespec start, end, iter_start, iter_end; // Declare timers for overall and per-iteration
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double elapsed_total = 0;

    // ---- Precomputations
    EMContext *ctx = createEMContext(X, W, q_method, *inputParams); // Allocate the important variables
    getInitialP(ctx, p_method, probMatrix);                         // Get the initial probabilities
    QMethodConfig config = getQMethodConfig(q_method, *inputParams);
    double newLL;
    double oldLL = -DBL_MAX;
    // ---...--- //

    // ---- Check the case where there's a group without voters ---- //
    int invalidGroup = checkGroups(*ctx);
    if (invalidGroup != -1)
    {
        removeColumn(&ctx->W, invalidGroup);
        EMContext *newCtx =
            EMAlgoritm(&ctx->X, &ctx->W, "group_proportional", q_method, convergence, LLconvergence, maxIter,
                       maxSeconds, verbose, time, iterTotal, logLLarr, finishing_reason, probMatrix, inputParams);
        cleanup(ctx);
        addColumnOfZeros(&newCtx->W, invalidGroup);
        // setParameters(X, W);
        addRowOfNaN(&newCtx->probabilities, invalidGroup);
        return newCtx;
        // TODO: Free the context
    }
    // ---...--- //
    Matrix oldProbabilities = createMatrix(ctx->G, ctx->C);
    // ---...--- //
    // ---- Execute the EM-iterations ---- //
    for (int i = 0; i < maxIter; i++)
    {
        // ---- Timer for the current iteration
        clock_gettime(CLOCK_MONOTONIC, &iter_start);
        *iterTotal = i;
        // config.params.iters = i; // Update the iteration number in the parameters
        ctx->iteration = i;
        config.computeQ(ctx, config.params, &newLL);
        // ---- Project Q if needed ---- //
        if (inputParams->prob_cond_every)
        {
            if (strcmp(inputParams->prob_cond, "project_lp") == 0)
                projectQ(ctx, *inputParams);
            else if (strcmp(inputParams->prob_cond, "lp") == 0)
                for (int b = 0; b < TOTAL_BALLOTS; b++)
                    LPW(ctx, b);
        }
        // ---...--- //

        memcpy(oldProbabilities.data, ctx->probabilities.data,
               sizeof(double) * oldProbabilities.rows * oldProbabilities.cols);
        *logLLarr = newLL;
        getP(ctx); // M-Step

        if (verbose)
        {
            Rprintf("\n----------\nIteration: %d\nProbability matrix:\n", i + 1);
            printMatrix(&ctx->probabilities);
            Rprintf("Log-likelihood: %f\n", newLL);
            if (i != 0)
                Rprintf("Delta log-likelihood: %f\n", fabs(newLL - oldLL));
        }
        // ---...--- //
        /*
         * For avoiding loops between same iterations (such as in the case of mvn_cdf), we impose that the
         * log-likelihood shouldn't decrease from the 50th iteration and on.
         */
        bool decreasing = oldLL > newLL ? true : false;
        bool early_stop = decreasing && strcmp(q_method, "exact") == 0 ? true : false;

        // ---- Check convergence ---- //
        if (i >= 1 && i >= config.params.miniter &&
            (fabs(newLL - oldLL) < LLconvergence ||
             convergeMatrix(&oldProbabilities, &ctx->probabilities, convergence) || early_stop))
        {
            // ---- End timer ----
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);

            if (verbose)
            {
                Rprintf("Converged after %d iterations (log-likelihood: %.4f) in %.5f "
                        "seconds.\n",
                        i + 1, newLL, elapsed_total);
            }
            *finishing_reason = 0;
            goto results;
        }
        // ---- Convergence wasn't found
        // ---- Stop the timer for verbose calls that aren't related to the algorithm
        clock_gettime(CLOCK_MONOTONIC, &iter_end);
        double elapsed_iter = (iter_end.tv_sec - iter_start.tv_sec) + (iter_end.tv_nsec - iter_start.tv_nsec) / 1e9;
        elapsed_total += elapsed_iter;
        R_CheckUserInterrupt();

        if (verbose)
            Rprintf("Elapsed time: %f\n----------\n", elapsed_total);

        // ---- The maximum time was reached
        if (elapsed_total >= maxSeconds)
        {
            if (verbose)
                Rprintf("Time limit reached.\n");
            *finishing_reason = 1;
            goto results;
        }

        oldLL = newLL;
    }
    // ---- Handle case where maxiter is achieved ----
    if (verbose)
        Rprintf("Maximum number of iterations reached without convergence.\n");

    *finishing_reason = 2;
    // ---...--- //
results:
    config.computeQ(ctx, config.params, &newLL);
    if (strcmp(inputParams->prob_cond, "project_lp") == 0)
    {
        projectQ(ctx, *inputParams);
        getP(ctx); // M-Step
    }
    else if (strcmp(inputParams->prob_cond, "lp") == 0)
    {
        for (int b = 0; b < TOTAL_BALLOTS; b++)
            LPW(ctx, b);
        getP(ctx); // M-Step
    }
    getPredictedVotes(ctx); // Compute the predicted votes for each ballot box

    *logLLarr = newLL;
    *time = elapsed_total;
    return ctx;
}

// ---- Clean all of the global variables ---- //
// __attribute__((destructor)) // Executes when the library is ready
void cleanup(EMContext *ctx)
{
    if (ctx == NULL)
        return;

    TOTAL_VOTES = 0;
    TOTAL_BALLOTS = 0;
    TOTAL_CANDIDATES = 0;
    TOTAL_GROUPS = 0;

    if (ctx->candidates_votes != NULL)
    {
        Free(ctx->candidates_votes);
        ctx->candidates_votes = NULL;
    }
    if (ctx->group_votes != NULL)
    {
        Free(ctx->group_votes);
        ctx->group_votes = NULL;
    }
    if (ctx->ballots_votes != NULL)
    {
        Free(ctx->ballots_votes);
        ctx->ballots_votes = NULL;
    }
    if (ctx->inv_ballots_votes != NULL)
    {
        Free(ctx->inv_ballots_votes);
        ctx->inv_ballots_votes = NULL;
    }

    if (ctx->X.data != NULL)
    { // Note that the columns and rows are usually stack.
        freeMatrix(&ctx->X);
        ctx->X.data = NULL;
    }
    if (ctx->W.data != NULL)
    {
        freeMatrix(&ctx->W);
        ctx->W.data = NULL;
    }
    if (ctx->Wnorm.data != NULL)
    { // <--- FALTABA
        freeMatrix(&ctx->Wnorm);
        ctx->Wnorm.data = NULL;
    }
    if (ctx->intW.data != NULL)
    {
        freeMatrixInt(&ctx->intW);
        ctx->intW.data = NULL;
    }
    if (ctx->intX.data != NULL)
    {
        freeMatrixInt(&ctx->intX);
        ctx->intX.data = NULL;
    }
    if (ctx->probabilities.data != NULL)
    {
        freeMatrix(&ctx->probabilities);
        ctx->probabilities.data = NULL;
    }

    if (ctx->q != NULL)
    {
        Free(ctx->q);
        ctx->q = NULL;
    }
    if (ctx->cdf_seeds != NULL)
    {
        Free(ctx->cdf_seeds);
        ctx->cdf_seeds = NULL;
    }
    if (ctx->predicted_votes != NULL)
    {
        Free(ctx->predicted_votes);
        ctx->predicted_votes = NULL;
    }

    // ---- OmegaSet **
    if (ctx->omegaset != NULL)
    {
        for (uint32_t b = 0; b < ctx->B; b++)
        {
            if (ctx->omegaset[b] != NULL)
            {
                OmegaSet *os = ctx->omegaset[b];
                if (os->data != NULL)
                {
                    for (size_t s = 0; s < os->size; s++)
                    {
                        freeMatrixInt(&os->data[s]); // deja os->data[s].data = NULL
                    }
                    Free(os->data);
                    os->data = NULL;
                }
                if (os->counts != NULL)
                { // <--- FALTABA
                    Free(os->counts);
                    os->counts = NULL;
                }
                Free(os);
                ctx->omegaset[b] = NULL;
            }
        }
        Free(ctx->omegaset);
        ctx->omegaset = NULL;
    }

    // ---- multinomial [b][s]
    if (ctx->multinomial != NULL)
    {
        for (uint32_t b = 0; b < ctx->B; b++)
        {
            if (ctx->multinomial[b] != NULL)
            {
                Free(ctx->multinomial[b]);
                ctx->multinomial[b] = NULL;
            }
        }
        Free(ctx->multinomial);
        ctx->multinomial = NULL;
    }

    if (ctx->logGamma != NULL)
    {
        Free(ctx->logGamma);
        ctx->logGamma = NULL;
    }

    // ---- Qconstant [b][s]
    if (ctx->Qconstant != NULL)
    {
        for (uint32_t b = 0; b < ctx->B; b++)
        {
            if (ctx->Qconstant[b] != NULL)
            {
                Free(ctx->Qconstant[b]);
                ctx->Qconstant[b] = NULL;
            }
        }
        Free(ctx->Qconstant);
        ctx->Qconstant = NULL;
    }

    // ---- hset: debe ser B*G si lo indexás como [b * G + g]
    if (ctx->hset != NULL)
    {
        size_t total = (size_t)ctx->B * (size_t)ctx->G; // asegurate que así se alocó
        for (size_t idx = 0; idx < total; ++idx)
        {
            Set *s = &ctx->hset[idx];
            if (s->data != NULL)
            {
                for (size_t i = 0; i < s->size; ++i)
                {
                    if (s->data[i] != NULL)
                    {
                        Free(s->data[i]);
                        s->data[i] = NULL;
                    }
                }
                Free(s->data);
                s->data = NULL;
                s->size = 0;
            }
        }
        Free(ctx->hset);
        ctx->hset = NULL;
    }

    // ---- kset: idem B*G
    if (ctx->kset != NULL)
    {
        size_t total = (size_t)ctx->B * (size_t)ctx->G; // asegurate que así se alocó
        for (size_t idx = 0; idx < total; ++idx)
        {
            Set *s = &ctx->kset[idx];
            if (s->data != NULL)
            {
                for (size_t i = 0; i < s->size; ++i)
                {
                    if (s->data[i] != NULL)
                    {
                        Free(s->data[i]);
                        s->data[i] = NULL;
                    }
                }
                Free(s->data);
                s->data = NULL;
                s->size = 0;
            }
        }
        Free(ctx->kset);
        ctx->kset = NULL;
    }

    // Free the context itself
    Free(ctx);
}
