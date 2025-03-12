#include "dynamic_program.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <dirent.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#ifndef BLAS_INT
#define BLAS_INT int
#endif

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
Matrix testBootstrap(double *quality, const char *set_method, const Matrix *xmat, const Matrix *wmat,
                     const int *boundaries, int A, int bootiter, const char *q_method, const char *p_method,
                     const double convergence, const double log_convergence, const int maxIter, const double maxSeconds,
                     QMethodInput inputParams)
{

    // ---- Merge within macrogroups ---- //
    Matrix mergedMat = A == wmat->cols ? *wmat : mergeColumns(wmat, boundaries, A); // Boundaries is of length A
    // ---...--- //

    // ---- Obtain the bootstrapped results ---- //
    Matrix standardMat = bootstrapA(xmat, &mergedMat, bootiter, q_method, p_method, convergence, log_convergence,
                                    maxIter, maxSeconds, false, inputParams);
    // ---...--- //

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
    const Matrix *xmat, const Matrix *wmat,

    // ---- Results
    int *results, // Array with cutting indices
    int *cuts,    // Amount of cuts

    // ---- EM and Bootstrap parameters
    double set_threshold, const char *set_method, int bootiter, const char *p_method, const char *q_method,
    const double convergence, const double log_convergence, const int maxIter, double maxSeconds, const bool verbose,
    QMethodInput inputParams)
{

    // ---- Define initial parameters ---- //
    double bestValue = DBL_MAX;
    Matrix bestMatrix, bootstrapMatrix;
    Matrix lastReward = buildRewards(xmat, wmat, wmat->cols);
    int *boundaries;
    double quality;
    // ---...--- //

    // ---- Loop through all possible macrogroups, starting from |G| cuts to 2 ---- //
    for (int i = wmat->cols; i > 1; i--)
    { // --- For every macrogroup cut
        // --- Base case, try with |A| = |G|, basically, get the bootstrap of the whole matrix.
        if (verbose)
        {
            Rprintf("\nCalculating %d macro-groups\n", i);
        }
        double bestVal;
        boundaries = solveDP(wmat->cols, i, &lastReward, &bestVal);

        if (verbose)
        {
            Rprintf("Optimal actions:\t[");
            for (int k = 0; k < i - 1; k++)
            {
                // Sum 1 to the index for using R's indexing
                Rprintf("%d, ", boundaries[k] + 1);
            }
            Rprintf("%d]\n", boundaries[i - 1]);
            Rprintf("Objective function:\t%.4f\n", bestVal);
        }
        // ---- Calculate the bootstrap matrix according the cutting boundaries
        bootstrapMatrix = testBootstrap(&quality, set_method, xmat, wmat, boundaries, i, bootiter, q_method, p_method,
                                        convergence, log_convergence, maxIter, maxSeconds, inputParams);
        if (verbose)
        {
            Rprintf("Bootstrapped matrix:\n");
            printMatrix(&bootstrapMatrix);
            Rprintf("Threshold value:\t%.4f\n", quality);
        }
        // --- Case it converges
        if (quality <= set_threshold)
        {
            for (int b = 0; b < i; b++)
            {
                results[b] = boundaries[b];
            }
            *cuts = i;
            Free(boundaries);
            return bootstrapMatrix;
        }
        // --- Case it is a better candidate than before
        if (quality < bestValue)
        {
            for (int b = 0; b < i; b++)
            {
                results[b] = boundaries[b];
            }
            Free(boundaries);
            bestMatrix = bootstrapMatrix;
            *cuts = i;
            bestValue = quality;
        }
        else
        {
            Free(boundaries);
            freeMatrix(&bootstrapMatrix);
        }
    }
    freeMatrix(&lastReward);
    if (verbose)
    {
        Rprintf("\nThe maximum threshold value was not accomplished. Returning the results of having %d macro-groups, "
                "having an statistic of %.4f, corresponding to the lesser threshold.\n",
                *cuts, bestValue);
    }
    // ---...--- //
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
static double g_bestLogLikelihood = -DBL_MAX;
static double *g_bestq = NULL;
static Matrix *g_bestMat = NULL;
static double g_besttime = 0;
static int g_bestFinishReason = 0;
static int g_bestIterTotal = 0;
static int *g_bestBoundaries = NULL; // will hold best partition indices
static int g_bestGroupCount = 0;     // how many groups in best partition

/*
   Recursive function that enumerates and evaluates all partitions of [0..G-1] into contiguous
   subranges.

   - `start` tells us where the current block begins
   - `currentBoundaries` up to `currentSize-1` are the "cutting indices" so far,
     meaning each boundary is the *end-index* of one subrange.
       e.g. if the subrange is [s..b], we store b in currentBoundaries.
   - Once `start == G`, we have a complete partition to evaluate.
*/
static void enumerateAllPartitions(int start, const int G, int *currentBoundaries, int currentSize, Matrix *xmat,
                                   const Matrix *wmat, const char *p_method, const char *q_method, double convergence,
                                   double log_convergence, bool verbose, int maxIter, double maxSeconds,
                                   QMethodInput inputParams)
{
    // ---- BASE CASE: We have closed a combination ---- //
    if (start == G)
    {
        // ---- Build the matrix according the partition
        // ---- Comment this line IF we want to account the matrix of group size 1
        if (currentSize == 1)
            return;
        Matrix merged = mergeColumns(wmat, currentBoundaries, currentSize);

        // ---- Run the EM Algorithm
        setParameters(xmat, &merged);

        Matrix initP = getInitialP(p_method);

        double timeUsed = 0.0;
        double logLLs[maxIter]; // TODO: Change this when the array stops being required
        double *qvals = NULL;
        int finishingReason = 0, totalIter = 0;

        Matrix finalP = EMAlgoritm(&initP, q_method, convergence, log_convergence, maxIter, maxSeconds, false,
                                   &timeUsed, &totalIter, logLLs, &qvals, &finishingReason, inputParams);

        double currentLL = (totalIter > 0) ? logLLs[totalIter - 1] : -DBL_MAX;
        Rprintf("current LL:%f\n", currentLL);
        // ---- Save the results if the value is better
        if (currentLL > g_bestLogLikelihood)
        {
            // -- Free previous best data if it exists -- //
            if (g_bestMat != NULL)
            {
                freeMatrix(g_bestMat);
                Free(g_bestMat);
                g_bestMat = NULL;
            }
            if (g_bestq != NULL)
            {
                Free(g_bestq);
                g_bestq = NULL;
            }
            if (g_bestBoundaries != NULL)
            {
                Free(g_bestBoundaries);
                g_bestBoundaries = NULL;
            }

            g_bestLogLikelihood = currentLL;

            g_bestGroupCount = currentSize;
            g_bestq = qvals;
            g_bestFinishReason = finishingReason;
            g_besttime = timeUsed;
            g_bestIterTotal = totalIter;
            g_bestMat = (Matrix *)Calloc(1, Matrix);
            *g_bestMat = finalP; // shallow copy of the struct, since it's on the stack buffer
            // g_bestq = qvals;
            qvals = NULL; // so we don't free it below if it's our best

            g_bestBoundaries = (int *)Calloc(currentSize, int);
            if (verbose)
                Rprintf("\n----------\nFound a feasible candidate, with a log-likelihood of %.6f and the optimal cuts "
                        "being:\n[",
                        g_bestLogLikelihood);
            for (int i = 0; i < currentSize; i++)
            {
                g_bestBoundaries[i] = currentBoundaries[i];
                // Sum 1 to the index for using R's indexing
                if (verbose && i != currentSize - 1)
                    Rprintf("%d, ", g_bestBoundaries[i] + 1);
                if (verbose && i == currentSize - 1)
                    Rprintf("%d]\n", g_bestBoundaries[i] + 1);
            }
            if (verbose)
            {
                Rprintf("- Candidate group size: %d\n- Candidate time: %f\n- Candidate total iterations: %d\n- "
                        "Candidate finish reason: %d\n- Candidate probability matrix:\n",
                        g_bestGroupCount, g_besttime, g_bestIterTotal, g_bestFinishReason);
                printMatrix(g_bestMat);
                Rprintf("----------\n");
            }
        }
        else
        {
            freeMatrix(&finalP);
            if (qvals != NULL)
            {
                Free(qvals);
                qvals = NULL;
            }
        }

        // ---- Clean every allocated memory ---- //
        cleanup();
        if (strcmp(q_method, "exact") == 0)
        {
            cleanExact();
        }
        else if (strcmp(q_method, "hnr") == 0)
        {
            cleanHitAndRun();
        }
        freeMatrix(&initP);
        // free the merged aggregator:
        freeMatrix(&merged);

        return;
    }

    // ---- RECURSION: There are still combinations to be made ---- //
    for (int end = start; end < G; end++)
    {
        // ---- We'll define a subrange [start..end]. We store 'end' as the boundary:
        /*
        if (currentSize == 1)
        {
            Rprintf("Skipping partition where group size is 1: start=%d, end=%d\n", start, end);
            continue; // Skip this iteration
        }
        */
        currentBoundaries[currentSize] = end;

        // ---- Recurse with the next subrange starting at 'end+1':
        enumerateAllPartitions(end + 1, G, currentBoundaries, currentSize + 1, xmat, wmat, p_method, q_method,
                               convergence, log_convergence, verbose, maxIter, maxSeconds, inputParams);
    }
}

/*
  Main “brute force” function that tries *every* possible grouping of columns [0..G-1].
  It returns the aggregator (merged) matrix *for the best log-likelihood partition*,
  and it also writes out the boundaries of that partition into `results[]` and the
  number of macro-groups into `*cuts`.
*/
Matrix aggregateGroupsExhaustive(
    // ---- Matrices
    Matrix *xmat, const Matrix *wmat,

    // ---- Partition boundaries
    int *results, // Cutting array
    int *cuts,    // Amount of cuts

    // ---- EM parameters
    const char *p_method, const char *q_method, double convergence, double log_convergence, bool verbose, int maxIter,
    double maxSeconds, QMethodInput inputParams,

    // ---- Out parameters
    double *outBestLL, double **outBestQ, double *outBestTime, int *outFinishReason, int *outIterTotal)
{
    const int G = wmat->cols;

    // Reset global tracking of best solution:
    g_bestLogLikelihood = -DBL_MAX;
    if (g_bestBoundaries)
    {
        Free(g_bestBoundaries);
        g_bestBoundaries = NULL;
    }
    if (g_bestMat)
    {
        freeMatrix(g_bestMat);
    }
    if (g_bestq)
    {
        Free(g_bestq);
    }
    g_bestGroupCount = 0;
    g_bestIterTotal = 0;
    g_bestFinishReason = 0;
    g_besttime = 0;

    // --- Array for enumerating possible partitions
    int *tempBoundaries = (int *)Calloc(G, int);

    // ---- Calls the recursion and update the global variables according the best values
    enumerateAllPartitions(0, // start from column 0
                           G, tempBoundaries,
                           0, // current partition is empty initially
                           xmat, wmat, p_method, q_method, convergence, log_convergence, verbose, maxIter, maxSeconds,
                           inputParams);

    Free(tempBoundaries);

    // g_bestBoundaries[] has the best partition, and g_bestGroupCount is how many groups.
    // --- Case where any group is feasible. I'm not even sure if it's possible...
    if (g_bestGroupCount == 0)
    {
        // Means we never found anything feasible or G=0 case.
        *cuts = 0;
        if (results)
        {
            results[0] = -1; // or some sentinel
        }
        Matrix empty = createMatrix(0, 0);
        return empty;
    }

    // Copy out the best partition boundaries:
    *cuts = g_bestGroupCount;
    for (int i = 0; i < g_bestGroupCount; i++)
    {
        results[i] = g_bestBoundaries[i];
    }

    // ---- Return the best EM data via arguments ----
    if (outBestLL)
    {
        *outBestLL = g_bestLogLikelihood;
    }
    if (outBestQ)
    {
        *outBestQ = g_bestq;
        // Clear the global pointer so we don’t double-free
        g_bestq = NULL;
    }
    if (outBestTime)
    {
        *outBestTime = g_besttime;
    }
    if (outFinishReason)
    {
        *outFinishReason = g_bestFinishReason;
    }
    if (outIterTotal)
    {
        *outIterTotal = g_bestIterTotal;
    }
    // Clean globals and make a shallow copy
    // Free(g_bestMat);
    // g_bestMat = NULL;
    Free(g_bestBoundaries);
    g_bestBoundaries = NULL;
    g_bestGroupCount = 0;
    Matrix returnMat = copyMatrix(g_bestMat);
    freeMatrix(g_bestMat);

    return returnMat;
}
