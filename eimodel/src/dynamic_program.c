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

/*
 * Computes, via Dynamic Programming, the optimal slicing for group aggregations towards standard deviations.
 *
 */
/*
int *idealSet(const Matrix *xmat, const Matrix *wmat, int A)
{

    int G = wmat->cols;
    // Call the reward table:
    Matrix reward = buildRewards(xmat, wmat, G);
    Rprintf("The reward matrix is:\n");
    printMatrix(&reward);
    Matrix choice = createMatrix(G + 1, A + 1);
    // Let DP be a (t x a) matrix that defines the maximum reward if we partition the first 't' groups in 'a'
    // macrogroups.
    // As base conditions;
    // 1. DP[0,0] = 0. Partition 0 groups in 0 macrogroups is null.
    // 2. DP[0,a] = -inf. It's impossible to make macrogroups without partitioning
    // 3. DP[t,0] = -inf. It's impossible to make 0 macrogroups if we make a partition.
    Matrix DP = createMatrix(G + 1, A + 1);
    fillMatrix(&DP, -1.0e9);
    MATRIX_AT(DP, G, 0) = 0.0;
    // IF DP[G][a>0], -\infty

    for (int t = G - 1; t >= 0; t--)
    {
        for (int a = 1; a <= A; a++)
        {
            if (a == 1)
            { // Only one macro group left
                MATRIX_AT(DP, t, a) = MATRIX_AT(reward, t, G - 1);
                MATRIX_AT(choice, t, a) = (double)G;
            }
            else
            { // This will ONLY execute when a \in [2,A], i.e, I can make a decition to either close or open
                double max_value = -1.0e9;
                int best_choice = -1;
                for (int k = t + 1; k <= G; k++)
                { // ---- For each group CONNECTION that can be made
                    // i.e, this will check the reward for closing the current group to a given
                    // group (must be bigger)
                    double candidate = MATRIX_AT(reward, t, k - 1) + MATRIX_AT(DP, k, a - 1);
                    // DP[k][a-1] IS THE + REMAINING OF CLOSING THE REST
                    // This will set the maximum possible value, i.e, the best choice of closing
                    if (candidate > max_value)
                    {
                        max_value = candidate;
                        best_choice = k;
                    }
                }
                MATRIX_AT(DP, t, a) = max_value;
                MATRIX_AT(choice, t, a) = best_choice;
            }
        }
    }

    Rprintf("\nThe DP matrix is:\n");
    printMatrix(&DP);
    Rprintf("\nThe choice matrix is:\n");
    printMatrix(&choice);
    int *boundaries = Calloc(A, int);
    int t0 = 0;
    for (int a = A; a > 0; a--)
    {
        // boundary from choice[t0,a0]
        int k = (int)MATRIX_AT(choice, t0, a);
        // store in boundaries[a0-1]
        boundaries[a - 1] = k;
        t0 = k; // next segment starts at k
        if (t0 >= G)
        {
            break;
        }
    }
    freeMatrix(&DP);
    freeMatrix(&choice);
    freeMatrix(&reward);

    return boundaries;
}
*/
/*
Matrix bellman(Matrix state, Matrix decition, Matrix lastReward, int lt, int ut, int t, int G, int A) {


    if (t == G) {
        return MATRIX_AT(lastReward, G-lt, G-1);
    }
    else {

    }

}
*/
/**
 * DP made with memoization
 *   s: start index of the current macro-group (0-based)
 *   t: current single-age-range index we are deciding on. Basically, the current state where we want to see if we'll
 * close. u: how many macro-groups have been closed so far. Note that u <= A. G: total # of single-age-ranges A: total #
 * of macro-groups we want
 *
 * lastReward is a matrix of the values of the std deviations from closing group 'i' to 'j'.
 *
 * memo[] and used[] are 1D arrays sized for (G+1)*(G+1)*(A+1).
 * action[] is likewise sized and will store 0 or 1 (open/close).
 */
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
int *idealSet(const Matrix *xmat, const Matrix *wmat, int A)
{
// There are 1 to G-1 stages 't'.
// Each stage represented by the tuple s = (l, u), where:
// l \in \lbrace 1,\dots, \min\lbrace t, G-1\rbrace\rbrace
    // u \in \lbrace 1, \dots, A-1\rbrace
    // l is the number of CONSECUTIVE INDIVIDUAL RANGES from the last open macro group
    // u is the number of macro groups already formed
    // The action is a binary that takes 1 if group `t` closes the last open macro group

    // Amount of groups
int G = wmat->cols;
    // Defining the state matrix
    Matrix state = createMatrix(G, A); // Even though l goes to \min (t, G-1)
    MATRIX_AT(state, 0, 0) =
    //
    // The utility from stage t given st and action at is
    // Vt(st, at) = \sigma_{t-lt+1, t}\cdot at + V_t*(st+1(st, at))


}
*/
/*
 * Obtain the bootstrapping values of the group aggregations and the convergence value
 *
 */
Matrix testBootstrap(double *quality, const char *set_method, const Matrix *xmat, const Matrix *wmat,
                     const int *boundaries, int A, int bootiter, const char *q_method, const char *p_method,
                     const double convergence, const int maxIter, const double maxSeconds, QMethodInput inputParams)
{

    // ---- Merge within macrogroups ---- //
    Matrix mergedMat = A == wmat->cols ? *wmat : mergeColumns(wmat, boundaries, A); // Boundaries is of length A
    // ---...--- //

    // ---- Obtain the bootstrapped results ---- //
    Matrix standardMat = bootstrapA(xmat, &mergedMat, bootiter, q_method, p_method, convergence, maxIter, maxSeconds,
                                    false, inputParams);
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
 * @param[in, out] results An array with the slicing indices.
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
Matrix aggregateGroups(const Matrix *xmat, const Matrix *wmat, int *results, int *cuts, double set_threshold,
                       const char *set_method, int bootiter, const char *p_method, const char *q_method,
                       const double convergence, const int maxIter, double maxSeconds, const bool verbose,
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
                                        convergence, maxIter, maxSeconds, inputParams);
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
