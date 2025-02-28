#include "bootstrap.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Corregir, no llama bootstrap....
// NOTE: Do not use a pointer on wmat so we can remove columns
double getSigmaForRange(const char *set_method, const Matrix *xmat, Matrix wmat, int g1, int g2, int bootiter,
                        const char *q_method, const char *p_method, const double convergence, const int maxIter,
                        const double maxSeconds, QMethodInput inputParams)
{

    // ---- Remove columns that do not belong [g1, g2] ---- //
    for (int i = 0; i < wmat.cols; i++)
    {
        if (i < g1 || i > g2)
            removeColumn(&wmat, i);
    }
    // ---...--- //

    Matrix standardMat =
        bootstrapA(xmat, &wmat, bootiter, q_method, p_method, convergence, maxIter, maxSeconds, false, inputParams);

    // ---- Maximum method ---- //
    if (strcmp(p_method, "maximum") == 0)
    {
        double maxval = 0;
        for (int j = 0; j < standardMat.rows; j++)
        {
            for (int k = 0; k < standardMat.cols; k++)
            {
                if (MATRIX_AT(standardMat, j, k) > maxval)
                    maxval = MATRIX_AT(standardMat, j, k);
            }
        }
        return maxval;
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
        return mean / (double)(standardMat.rows * standardMat.cols);
    }
    // ---...--- //
}

Matrix buildRewards(const Matrix *xmat, const Matrix *wmat, int setSize, double set_threshold, const char *set_method,
                    int bootiter, const char *q_method, const char *p_method, const double convergence,
                    const int maxIter, const double maxSeconds, QMethodInput inputParams)
{
    // Table to store the rewards values
    Matrix table = createMatrix(setSize, setSize);

    // Only store the pairs where g1 <= g2
    for (int g1 = 0; g1 < setSize; g1++)
    {
        for (int g2 = g1; g2 < setSize; g2++)
        {
            MATRIX_AT(table, g1, g2) = getSigmaForRange(set_method, xmat, *wmat, g1, g2, bootiter, q_method, p_method,
                                                        convergence, maxIter, maxSeconds, inputParams);
        }
    }

    return table;
}

void idealSet(const Matrix *xmat, const Matrix *wmat, int G, int A, double set_threshold, const char *set_method,
              int bootiter, const char *q_method, const char *p_method, const double convergence, const int maxIter,
              const double maxSeconds, const bool verbose, QMethodInput inputParams)

{

    // Call the reward table:
    Matrix reward = buildRewards(xmat, wmat, G, set_threshold, set_method, bootiter, q_method, p_method, convergence,
                                 maxIter, maxSeconds, inputParams);

    Matrix decition = createMatrix(G, A);
    // Let DP be a (t x a) matrix that defines the maximum reward if we partition the first 't' groups in 'a'
    // macrogroups.
    // As base conditions;
    // 1. DP[0,0] = 0. Partition 0 groups in 0 macrogroups is null.
    // 2. DP[0,a] = -inf. It's impossible to make macrogroups without partitioning
    // 3. DP[t,0] = -inf. It's impossible to make 0 macrogroups if we make a partition.
    Matrix DP = createMatrix(G, A);
    fillMatrix(&DP, -1.0e9);
    MATRIX_AT(DP, G, 0) = 0.0;
    // IF DP[G][a>0], -\infty

    int best_choice = -1;
    for (int t = G - 1; t >= 0; t--)
    {
        for (int a = 1; a <= A; a++)
        {

            if (a == 1)
            { // Only one macro group left
                MATRIX_AT(DP, t, 1) = MATRIX_AT(reward, t, G - 1);
            }
            else
            { // This will ONLY execute when a \in [2,A], i.e, I can make a decition to either close or open
                double max_value = -1.0e9;
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
                MATRIX_AT(decition, t, a) = best_choice;
            }
        }
    }
    if (verbose)
        Rprintf("The best macro groups for a set size %d are:\t", A);
    int finalArr[A];
    for (int i = 0; i < A; i++)
    {
        finalArr[i] = (int)MATRIX_AT(decition, i, A);
        if (verbose)
            Rprintf("%d, ", finalArr[i]);
    }

    // The choices would be the maximum values from DP column wise.

    /*
    for (int t = 0; t < G; t++)
    {
        for (int g = 0; g < G; g++)
        {
            for (int a = 0; a < A; a++)
            {
                MATRIX_AT(DP, 0, 0) = 0;
                MATRIX_AT(DP, g, a) = -1.0e9;

                MATRIX_AT(decition, g, a) = -1;
            }
        }
        DParray[t] = DP;
        decitionArray[t] = decition;
    }
*/
    /*
    // Starting from the edge, in T, DP = reward
    for (int t = G - 1; t > -1; t--)
    {
        for (int g = 0; g < G; g++)
        {
            // int *bestPartitions;
            // int tie = 2;
            // R_max_col(reward.data, &G, &G, bestPartitions, &tie);
            for (int a = 0; a < A; a++)
            {
                if (a == A - 1) // Get the maximum value from std
            }
        }
        DParray[t] = DP;
        decitionArray[t] = decition;
    }

    // For 1 macrogroup and 1 partition, the reward would be the biggest
    // NOTE: USE R_max_col
    int *bestPartitions;
    int tie = 2;
    R_max_col(reward.data, &G, &G, bestPartitions, &tie);
    DP[1, 1] = MATRIX_AT(reward, 0, G);
*/
    //
}
