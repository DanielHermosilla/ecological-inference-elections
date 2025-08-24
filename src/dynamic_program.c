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

#include "dynamic_program.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <dirent.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

double getSigmaForRange(const Matrix *xmat, const Matrix *wmat, int g1, int g2, double *ballotVotes)
{
    int ballotBoxes = wmat->rows;
    // --- Compute mean proportions across ballotBoxes --- //
    double *mean = (double *)Calloc(ballotBoxes, double);
    double mu = 0.0;
    double sum = 0.0;

    for (int b = 0; b < ballotBoxes; b++)
    {
        sum = 0.0;
        for (int g = g1; g <= g2; g++)
        {
            sum += MATRIX_AT_PTR(wmat, b, g);
        }
        mean[b] = (double)sum / (double)ballotVotes[b];
        mu += mean[b];
    }
    mu /= ballotBoxes;
    // --- Compute standard deviation of these proportions --- //
    double num = 0.0;
    for (int b = 0; b < ballotBoxes; b++)
    {
        double diff = mean[b] - mu;
        num += diff * diff;
    }
    Free(mean);

    double sqrted = R_pow(num / ballotBoxes, 0.5);
    return sqrted;
}
/*
 * Builds a [G x G] matrix with the rewards of incorporating a group determined within (g1, g2).
 *
 */
Matrix buildRewards(const Matrix *xmat, const Matrix *wmat, int setSize)
{
    // ---- Create table to store the rewards values ---- //
    Matrix table = createMatrix(setSize, setSize);
    // ---...--- //

    // ---- Store the amount of votes per ballot box ---- //
    double *ballotVotes = (double *)Calloc(wmat->rows, double);
    colSum(wmat, ballotVotes);
    // ---...--- //

    // ---- Calculate the value of closing a group from 'g1' to 'g2', only where g1 <= g2 ---- //
    for (int g1 = 0; g1 < setSize; g1++)
    {
        for (int g2 = g1; g2 < setSize; g2++)
        {
            MATRIX_AT(table, g1, g2) = getSigmaForRange(xmat, wmat, g1, g2, ballotVotes);
        }
    }
    // ---...--- //

    Free(ballotVotes);
    return table;
}

double dpReward(int s, int t, int u, int G, int A, const Matrix *lastReward, double *memo, bool *used, int *action)
{
    // Base case: if we've gone beyond the last index:
    if (t == G)
    {
        // This is valid, if and only if, all of the macrogroups have been closed
        // before.
        return (u == A) ? 0.0 : -DBL_MAX;
    }

    // ---- If we have already formed A groups but t < G => we can't form more. This would
    // be the case where the system "enforces" us to create more macrogroups.
    if (u == A)
    {
        // ---- no valid way to place the remaining single-ranges => -inf
        return -DBL_MAX;
    }
    // ---...--- //

    // Access memo/used
    double *memoRef = &Q_3D(memo, s, t, u, (G + 1), (A + 1));
    bool *usedRef = &Q_3D(used, s, t, u, (G + 1), (A + 1));
    int *actionRef = &Q_3D(action, s, t, u, (G + 1), (A + 1));

    // ---- If the value has been already used, return the memorizated value.
    if (*usedRef)
    {
        return *memoRef;
    }

    // ---- Base case: reaching the LAST group ---- //
    // Case where we are in the LAST group. We should close here, hopefully end up with A macrogroups
    // From here we should jump to t == G, with u+1 closes.
    if (t == G - 1)
    {
        // The objective value would be from the last non-closed index "s" to the last one "G-1 = t"
        double val = MATRIX_AT_PTR(lastReward, s, t);
        // Here we jump to t+1, we should know that if we had formed |A| macrogroups, then nextVal would be 0.
        double nextVal = dpReward(t + 1, t + 1, u + 1, G, A, lastReward, memo, used, action);
        // If nextVal is -infty then we had closed more groups -> it's and invalid combination
        double totalVal = (val == DBL_MAX || nextVal == -DBL_MAX) ? -DBL_MAX : (val + nextVal); // Vellman

        // Impose that this was already calculated
        *usedRef = true;
        // The optimal value when starting from 's', being in state 't' and having 'u' groups closed is 'totalVal'
        *memoRef = totalVal;
        // Closed at "t==G-1"
        *actionRef = 1; // '1' means "closed" at t.
        return totalVal;
    }
    // ---...--- //

    // ---- Choose between two possibilities ---- //
    //  (1) close the group at 't'
    //  (2) keep open (move on to t+1) if possible
    double bestVal = -DBL_MAX;
    int bestAct = -1; // 1=close, 0=open

    // ---- (1) "close" this macro-group at t => we collect reward for [s..t]. ---- //
    {
        // ---- 'valClose' would mean the sigma of the current Bellman combination, since at would be one
        double valClose = MATRIX_AT_PTR(lastReward, s, t);
        if (valClose > -DBL_MAX)
        {
            // ---- The next value would have a starting index 's' of the next state 't+1' and have formed 'u+1'
            // macrogroups.
            double valNext = dpReward(t + 1, t + 1, u + 1, G, A, lastReward, memo, used, action);
            // ---- This condition will check the 'future' value and later compare it to the (2) option (not closing)
            if (valNext > -DBL_MAX)
            {
                // ---- The candidate utility value; the current value and the next one
                double candidate = valClose + valNext;
                // ---- This will trivially be true, however, it would stay for clarity.
                if (candidate > bestVal)
                {
                    bestVal = candidate;
                    bestAct = 1;
                }
            }
        }
    }
    // ---...--- //

    // ----  (2) "open"/keep going => do NOT close, just move to t+1 in the same group ---- //
    // ---- Only do this if t+1 < G, but we already handled t == G-1 above, so we are good
    {
        // ---- The next value would have a starting index of 's' (since we're not closing, we maintain the index) and
        // next state of 't+1', maintaining the formed macrogroups.
        double valOpen = dpReward(s, t + 1, u, G, A, lastReward, memo, used, action);
        // ---- 'valOpen' would be the utility of not closing. If the utility is better than the best value (i.e,
        // closing), then the best action is to not close.
        if (valOpen > bestVal)
        {
            bestVal = valOpen;
            bestAct = 0;
        }
    }
    // ---...--- //

    // ---- Remember the computed values
    *usedRef = true;
    *memoRef = bestVal;
    *actionRef = bestAct;
    return bestVal;
}

/**
 * Reconstructs the solution path (which boundaries got closed).
 */
/**
 * collectCuts(s,t,u):
 *   Recursively follow the action[] decisions from (s,t,u),
 *   appending 't' to the cuts[] array whenever action=1 (i.e. close).
 *
 *   We'll pass in a pointer to an integer 'pos' that tracks
 *   how many closures have been recorded so far.
 *
 *   If dpReward was feasible for exactly A groups,
 *   we eventually get exactly A closure indices.
 */
void collectCuts(int s, int t, int u, int G, int A, const Matrix *lastReward, double *memo, bool *used, int *action,
                 int *cuts, // array of size at least A
                 int *pos   // how many closures so far
)
{
    // ---- If we pass t == G, we are done (or no more single-ranges).
    if (t == G)
    {
        return;
    }

    // ---- If t == G-1, we must have closed, it's a border condition:
    if (t == G - 1)
    {
        // ---- The DP forced a closure at t
        cuts[*pos] = t;
        (*pos)++;
        // next; we could finish here but it will remain this for clarity.
        collectCuts(t + 1, t + 1, u + 1, G, A, lastReward, memo, used, action, cuts, pos);
        return;
    }

    // ---- Check if the action was to close ---- //
    int act = Q_3D(action, s, t, u, (G + 1), (A + 1));
    if (act == 1)
    {
        // ---- Save the value
        cuts[*pos] = t;
        (*pos)++;

        // ---- Go to the next state => (t+1, t+1, u+1); start at 't+1' and group 't+1' with 'u+1' closes.
        collectCuts(t + 1, t + 1, u + 1, G, A, lastReward, memo, used, action, cuts, pos);
    }
    else
    {
        // ---- The action was to keep it open
        // ---- Go to next state => (s, t+1, u); start at 's' (same as before), on group 't+1' with 'u' closes.
        collectCuts(s, t + 1, u, G, A, lastReward, memo, used, action, cuts, pos);
    }
}

/**
 * solveDP(...) -> returns the array of cut indices of length A
 *   If no valid partition, returns NULL.
 *
 * 'lastReward' is GxG, where
 *   MATRIX_AT_PTR(lastReward, i, j) = reward for grouping [i..j].
 */
int *solveDP(int G, int A, const Matrix *lastReward,
             double *outBestValue // optional: to store the best total reward
)
{

    // ---- Create the 3D arrays ---- //
    int totalSize = (G + 1) * (G + 1) * (A + 1);
    // Table for remembering the past results
    double *memo = Calloc(totalSize, double);
    // Table for remembering if the value was used
    bool *used = Calloc(totalSize, bool);
    // Table for determining an "a_t" action.
    int *action = Calloc(totalSize, int);
    // ---...--- //

    // --- Initialize the arrays
    for (int i = 0; i < totalSize; i++)
    {
        memo[i] = 0.0;
        used[i] = false;
        action[i] = -1;
    }

    // ---- Compute best total reward from (s=0, t=0, u=0) ---- //
    double bestVal = dpReward(0, 0, 0, G, A, lastReward, memo, used, action);
    // ---...--- //

    // ---- For avoiding overflows, we look for -0.5 * DBL_MAX. Anyway, it would mean that
    // there are no valid partitions
    if (bestVal <= -0.5 * DBL_MAX)
    {
        // Means we got -DBL_MAX => no valid partition
        if (outBestValue)
        {
            *outBestValue = -DBL_MAX; // indicate invalid
        }
        Free(memo);
        Free(used);
        Free(action);
        return NULL;
    }

    // ---- Reconstruct the closure indices, i.e, the at actions ----
    int *cuts = Calloc(A, int); // We would expect A closings (including the last one)
    int pos = 0;

    collectCuts(0, 0, 0, G, A, lastReward, memo, used, action, cuts, &pos);

    // pos should be exactly A if the DP formed exactly A groups
    if (pos != A)
    {
        error("WARNING: we expected exactly %d closures, got %d. Something is off.\n", A, pos);
    }

    // 3) Store bestVal if desired
    if (outBestValue)
    {
        *outBestValue = bestVal;
    }

    // Cleanup DP
    Free(memo);
    Free(used);
    Free(action);

    // Return the array with A closure indices
    return cuts;
}

/*
 * Obtain the bootstrapping values of the group aggregations and the convergence value
 *
 */
Matrix testBootstrap(EMContext *ctx, double *quality, const char *set_method, const Matrix *xmat, const Matrix *wmat,
                     const int *boundaries, int A, int bootiter, const char *q_method, const char *p_method,
                     const double convergence, const double log_convergence, const int maxIter, const double maxSeconds,
                     QMethodInput inputParams)
{

    Matrix mergedMat, standardMat;
    // ---- Merge within macrogroups ---- //
    if (A != -1)
    {
        mergedMat = A == wmat->cols ? *wmat : mergeColumns(wmat, boundaries, A); // Boundaries is of length A
                                                                                 // ---...--- //
        // ---- Obtain the bootstrapped results ---- //
        GetRNGstate();
        standardMat = bootstrapA(xmat, &mergedMat, bootiter, q_method, p_method, convergence, log_convergence, maxIter,
                                 maxSeconds, false, &inputParams);
        PutRNGstate();
        // ---...--- //
    }
    else
    {
        mergedMat = createMatrix(wmat->rows, 1);
        for (int i = 0; i < (int)TOTAL_BALLOTS; i++)
        {
            MATRIX_AT(mergedMat, i, 0) = ctx->ballots_votes[i];
        }
        // printMatrix(&mergedMat);
        GetRNGstate();
        standardMat = bootSingleMat(xmat, &mergedMat, bootiter, false);
        PutRNGstate();
    }

    // ---- Maximum method ---- //
    if (strcmp(set_method, "maximum") == 0)
    {
        double maxval = maxElement(&standardMat);
        *quality = maxval;
    }
    // ---...--- //
    // ---- Mean method ---- //
    else
    {
        double mean = 0;
        for (int j = 0; j < standardMat.rows; j++)
        {
            for (int k = 0; k < standardMat.cols; k++)
            {
                mean += MATRIX_AT(standardMat, j, k);
            }
        }
        mean /= (double)(standardMat.rows * standardMat.cols);
        *quality = mean;
    }

    if (findNaN(&standardMat))
        *quality = INFINITY;

    return standardMat;
    // ---...--- //
}

/*
 * Main function to obtain the heuristical best group aggregation, using dynamic programming. Tries every combination
 * using a standard deviation approximate. Given the approximate, computes the bootstrapped standard deviation and
 * checks if it accomplishes the proposed threshold.
 *
 * @param[in] xmat The candidate (c x b) matrix.
 * @param[in] wmat The group (b x g) matrix.
 * @param[in, out results An array with the slicing indices.
 * @param[in] set_threshold The threshold of the proposed method
 * @param[in] set_method The method for evaluating the bootstrapping threshold.
 * @param[in] bootiter The amount of bootstrap iterations
 * @param[in] p_method The method for calculating the initial probability.
 * @param[in] q_method The method for calculating the EM algorithm of the boot samples.
 * @param[in] convergence The convegence threshold for the EM algorithm.
 * @param[in] maxIter The maximum amount of iterations to perform on the EM algorithm.
 * @param[in] maxSeconds The maximum amount of seconds to run the algorithm.
 * @param[in] verbose Boolean to whether verbose useful outputs.
 * @param[in] inputParams The parameters for specific methods.
 *
 * @return The heuristic optimal matrix with bootstrapped standard deviations.
 */
Matrix aggregateGroups(
    // ---- Matrices
    EMContext *ctx,

    // ---- Results
    int *results, // Array with cutting indices
    int *cuts,    // Amount of cuts
    bool *bestResult,

    // ---- EM and Bootstrap parameters
    double set_threshold, const char *set_method, bool feasible, int bootiter, const char *p_method,
    const char *q_method, const double convergence, const double log_convergence, const int maxIter, double maxSeconds,
    const bool verbose, QMethodInput *inputParams)
{
    Matrix *wmat = &ctx->W;
    Matrix *xmat = &ctx->X;

    // ---- Define initial parameters ---- //
    double bestValue = DBL_MAX;
    Matrix bestMatrix, bootstrapMatrix;
    Matrix lastReward = buildRewards(xmat, wmat, wmat->cols);
    int *boundaries;
    double quality;
    // ---...--- //

    // ---- Loop through all possible macrogroups, starting from |G| cuts to 2 ---- //
    for (int i = wmat->cols; i > 0; i--)
    { // --- For every macrogroup cut
        // --- Base case, try with |A| = |G|, basically, get the bootstrap of the whole matrix.
        if (verbose)
        {
            Rprintf("\n----------\nCalculating %d macro-groups\n", i);
        }
        double bestVal;
        // ---- Handle case where there's need to be a DP done ---- //
        if (i != 1)
        {
            boundaries = solveDP(wmat->cols, i, &lastReward, &bestVal);

            if (verbose)
            {
                Rprintf("Group aggregation:\t[");
                for (int k = 0; k < i - 1; k++)
                {
                    // Sum 1 to the index for using R's indexing
                    Rprintf("%d, ", boundaries[k] + 1);
                }
                Rprintf("%d]\n", boundaries[i - 1] + 1);
                Rprintf("Groups standard deviation:\t%f\n", bestVal);
            }
            // ---- Calculate the bootstrap matrix according the cutting boundaries
            bootstrapMatrix = testBootstrap(ctx, &quality, set_method, xmat, wmat, boundaries, i, bootiter, q_method,
                                            p_method, convergence, log_convergence, maxIter, maxSeconds, *inputParams);
        }
        // ---- Case where there's no cuts ---- //
        else
        {
            // int *boundaries = Calloc(2, int);
            bootstrapMatrix = testBootstrap(ctx, &quality, set_method, xmat, wmat, boundaries, -1, bootiter, q_method,
                                            p_method, convergence, log_convergence, maxIter, maxSeconds, *inputParams);
        }
        if (verbose && quality)
        {
            Rprintf("Standard deviation matrix:\n");
            printMatrix(&bootstrapMatrix);
            Rprintf("Statistic value:\t%.4f\n----------\n", quality);
        }
        // --- Case it converges
        if (quality <= set_threshold && quality != INFINITY)
        {
            for (int b = 0; (b < i) & (i != 1); b++)
            {
                results[b] = boundaries[b];
            }

            // Border cases
            *cuts = i;
            results[0] = i == 1 ? wmat->cols - 1 : results[0];

            if ((i != 1) & (boundaries != NULL))
            {
                Free(boundaries);
            }
            return bootstrapMatrix;
        }
        // --- Case it is a better candidate than before
        if (quality < bestValue && quality != INFINITY)
        {
            for (int b = 0; (b < i) & (i != 1); b++)
            {
                results[b] = boundaries[b];
            }
            if ((i != 1) & (boundaries != NULL))
                Free(boundaries);
            bestMatrix = bootstrapMatrix;
            results[0] = i == 1 ? wmat->cols - 1 : results[0];
            *cuts = i;
            bestValue = quality;
        }
        else
        {
            if (i != 1)
                Free(boundaries);
            freeMatrix(&bootstrapMatrix);
        }
    }
    freeMatrix(&lastReward);
    int totalMacrogroups = *cuts;
    totalMacrogroups += *cuts == -1 ? 2 : 0;
    if (verbose)
    {
        Rprintf("\nNo group aggregation yielded a standard deviation matrix statistic below the specified "
                "threshold. The aggregation with the lowest statistic was [");
        for (int k = 0; k < totalMacrogroups; k++)
        {
            Rprintf("%d", results[k] + 1);
            if (k != totalMacrogroups - 1)
                Rprintf(", ");
        }
        Rprintf("] with a value of %.4f — still above the threshold of %.4f. ", bestValue, set_threshold);
        if (!feasible)
        {
            Rprintf("If "
                    "you would like to retrieve this group aggregation despite its standard deviation matrix statistic "
                    "being above the threshold, set feasible = FALSE.\n");
        }
        else
        {
            Rprintf("Because "
                    "'feasibile' parameter is set to FALSE, the group aggregation will be returned anyway.\n");
        }
    }
    // ---...--- //
    *bestResult = true;
    return bestMatrix;
}

/**
 *
 * GREEDY APPROACH: We try every possible combination, it's of order O(2^{G-1})
 *
 */

/*
  Global variables to track the best parameters of the 'winning' EM.
*/
typedef struct
{
    double bestLogLikelihood;
    double *bestq;
    Matrix *bestMat;
    Matrix *bestBootstrap;
    double bestTime;
    int bestFinishReason;
    int bestIterTotal;
    int *bestBoundaries;
    double *bestExpected;
    int bestGroupCount;
} ExhaustiveResult;

typedef struct
{
    Matrix *xmat;
    Matrix *wmat;
    const char *set_method;
    int bootiter;
    const char *p_method;
    const char *q_method;
    double max_qual;
    double convergence;
    double log_convergence;
    bool verbose;
    int maxIter;
    double maxSeconds;
    QMethodInput inputParams;
} ExhaustiveOptions;

static void enumerateAllPartitions(int start, int G, int *currentBoundaries, int currentSize, ExhaustiveOptions *opts,
                                   ExhaustiveResult *res)
{
    if (start == G)
    {
        if (opts->verbose)
        {
            Rprintf("----------\nGroup aggregation:\t[");
            for (int k = 0; k < currentSize - 1; k++)
            {
                // Sum 1 to the index for using R's indexing
                Rprintf("%d, ", currentBoundaries[k] + 1);
            }
            Rprintf("%d]\n", currentBoundaries[currentSize - 1] + 1);
        }
        // ---- Build the matrix according the partition
        // ---- Comment this line IF we want to account the matrix of group size 1
        // if (currentSize == 1)
        //{ // The one aggregation case we'll do it with mult
        //    opts->q_method = "mult";
        //}
        else if (currentSize == 0)
            return;

        Matrix merged = mergeColumns(opts->wmat, currentBoundaries, currentSize);

        // ---- Run the EM Algorithm
        // setParameters(opts->xmat, &merged);

        // Matrix initP = getInitialP(opts->p_method);

        double timeUsed = 0.0;
        double logLLs = 0.0; // TODO: Change this when the array stops being required
        double *qvals = NULL;
        int finishingReason = 0, totalIter = 0;

        EMContext *ctx = EMAlgoritm(opts->xmat, &merged, opts->p_method, opts->q_method, opts->convergence,
                                    opts->log_convergence, opts->maxIter, opts->maxSeconds, false, &timeUsed,
                                    &totalIter, &logLLs, &finishingReason, &opts->inputParams);

        size_t n_q = (size_t)ctx->B * ctx->G * ctx->C;
        qvals = (double *)Calloc(n_q, double);
        memcpy(qvals, ctx->q, n_q * sizeof(double));
        // qvals = &ctx->q;

        double currentLL = (totalIter > 0) ? logLLs : -DBL_MAX;
        if (opts->verbose)
            Rprintf("Log-likelihood:\t%f\n", currentLL);

        // ---- Clean every allocated memory ---- //
        // cleanup(ctx);
        /*
        if (strcmp(opts->q_method, "exact") == 0)
        {
            cleanExact();
        }
        else if (strcmp(opts->q_method, "mcmc") == 0)
        {
            cleanHitAndRun();
        }
        */
        // else if (strcmp(opts->q_method, "mult") == 0)
        //{
        //    cleanMultinomial();
        //}
        // free the merged aggregator:
        // freeMatrix(&merged);

        // ---- Save the results if the value is better

        if (currentLL > res->bestLogLikelihood)
        {
            // --- Now it would be convenient to evaluate within the SD, later, goto notBestValue
            double qual;
            Matrix bootstrapedMat =
                testBootstrap(ctx, &qual, opts->set_method, opts->xmat, opts->wmat, currentBoundaries, currentSize,
                              opts->bootiter, opts->q_method, opts->p_method, opts->convergence, opts->log_convergence,
                              opts->maxIter, opts->maxSeconds, opts->inputParams);
            // freeMatrix(&bootstrapedMat);
            if (opts->verbose && qual != INFINITY)
                Rprintf("Standard deviation statistic:\t%f\n", qual);
            // compute bootstrap & qual …
            if (qual <= opts->max_qual && qual != INFINITY)
            {
                // free old best
                if (res->bestBootstrap)
                    freeMatrix(res->bestBootstrap), Free(res->bestBootstrap);
                if (res->bestMat)
                    freeMatrix(res->bestMat), Free(res->bestMat);
                if (res->bestq)
                    Free(res->bestq);
                if (res->bestExpected)
                    Free(res->bestExpected);
                if (res->bestBoundaries)
                    Free(res->bestBoundaries);

                // store new best
                res->bestLogLikelihood = currentLL;
                res->bestTime = timeUsed;
                res->bestFinishReason = finishingReason;
                res->bestIterTotal = totalIter;

                res->bestq = qvals;
                res->bestExpected = ctx->predicted_votes;
                ctx->predicted_votes = NULL; // prevent double‐free
                qvals = NULL;                // prevent double‐free

                // deep‐copy matrices
                res->bestMat = (Matrix *)Calloc(1, Matrix);
                *res->bestMat = copMatrix(&ctx->probabilities);
                res->bestBootstrap = (Matrix *)Calloc(1, Matrix);
                *res->bestBootstrap = copMatrix(&bootstrapedMat);
                freeMatrix(&bootstrapedMat);

                // copy boundaries
                res->bestGroupCount = currentSize;
                res->bestBoundaries = (int *)Calloc(currentSize, int);
                for (int i = 0; i < currentSize; i++)
                    res->bestBoundaries[i] = currentBoundaries[i];
            }
        }
        // else
        //  {
        //  cleanup non‐best branch
        //  freeMatrix(&finalP);
        cleanup(ctx);
        freeMatrix(&merged);
        // if (qvals)
        //    Free(qvals);
        //}
        return;
    }

    // recursion
    for (int end = start; end < G; end++)
    {
        currentBoundaries[currentSize] = end;
        enumerateAllPartitions(end + 1, G, currentBoundaries, currentSize + 1, opts, res);
    }
}

// 4. Rewrite your “driver” to use the struct instead of globals
Matrix aggregateGroupsExhaustive(Matrix *xmat, Matrix *wmat, int *results, int *cuts, const char *set_method,
                                 int bootiter, double max_qual, const char *p_method, const char *q_method,
                                 double convergence, double log_convergence, bool verbose, int maxIter,
                                 double maxSeconds, QMethodInput *inputParams, double *outBestLL, double **outBestQ,
                                 double **bestExpected, Matrix **bestBootstrap, double *outBestTime,
                                 int *outFinishReason, int *outIterTotal)
{
    int G = wmat->cols;

    // Initialize options and result
    ExhaustiveOptions opts = {.xmat = xmat,
                              .wmat = wmat,
                              .set_method = set_method,
                              .bootiter = bootiter,
                              .p_method = p_method,
                              .q_method = q_method,
                              .max_qual = max_qual,
                              .convergence = convergence,
                              .log_convergence = log_convergence,
                              .verbose = verbose,
                              .maxIter = maxIter,
                              .maxSeconds = maxSeconds,
                              .inputParams = *inputParams};

    ExhaustiveResult res = {.bestLogLikelihood = -DBL_MAX,
                            .bestq = NULL,
                            .bestMat = NULL,
                            .bestBootstrap = NULL,
                            .bestTime = 0.0,
                            .bestFinishReason = 0,
                            .bestIterTotal = 0,
                            .bestBoundaries = NULL,
                            .bestGroupCount = 0};

    // temp working buffer
    int *tempBoundaries = Calloc(G, int);

    // recurse
    enumerateAllPartitions(0, G, tempBoundaries, 0, &opts, &res);
    Free(tempBoundaries);

    // unpack result
    if (res.bestGroupCount == 0)
    {
        *cuts = 0;
        if (results)
            results[0] = -1;
        if (verbose)
            Rprintf("\nNo group aggregation yielded a standard deviation matrix statistic below the specified "
                    "threshold. "
                    "The aggregation with the lowest statistic was still above the threshold of %.4f. As a result, no "
                    "group aggregation is returned.",
                    max_qual);
        return createMatrix(1, 1);
    }

    *cuts = res.bestGroupCount;
    for (int i = 0; i < res.bestGroupCount; i++)
        results[i] = res.bestBoundaries[i];

    if (outBestLL)
        *outBestLL = res.bestLogLikelihood;
    if (outBestQ)
        *outBestQ = res.bestq, res.bestq = NULL;
    if (outBestTime)
        *outBestTime = res.bestTime;
    if (outFinishReason)
        *outFinishReason = res.bestFinishReason;
    if (outIterTotal)
        *outIterTotal = res.bestIterTotal;
    if (bestBootstrap)
        *bestBootstrap = res.bestBootstrap;
    if (bestExpected)
        *bestExpected = res.bestExpected, res.bestExpected = NULL;

    // copy & clean up
    Matrix bestCopy = copMatrix(res.bestMat);
    freeMatrix(res.bestMat);
    Free(res.bestMat);
    Free(res.bestBoundaries);

    return bestCopy;
}
