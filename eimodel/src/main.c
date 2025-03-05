/*
Copyright (c) 2025 Daniel Hermosilla

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

#undef I
// ---- Inititalize global variables ---- //
uint32_t TOTAL_VOTES = 0;
uint32_t TOTAL_BALLOTS = 0;
uint16_t TOTAL_CANDIDATES = 0;
uint16_t TOTAL_GROUPS = 0;
uint16_t *BALLOTS_VOTES = NULL;    // Total votes per ballot
uint32_t *CANDIDATES_VOTES = NULL; // Total votes per candidate
uint32_t *GROUP_VOTES = NULL;      // Total votes per group
double *inv_BALLOTS_VOTES = NULL;  // BALLOTS_VOTES^{-1}
Matrix *X = NULL;
Matrix *W = NULL;
// ---...--- //

/**
 * @brief Yields the global parameters of the process. Usually this should be done once for avoiding
 * computing a loop over ballots. It also changes the parameters in case it's called with other `x` and `w` matrix.
 *
 * Gets the total amount of votes in the process. This should only be donce once for avoiding computing loops
 * over ballots.
 *
 * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
 * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
 *
 * @return void. Will edit the static values in the file
 *
 * @note This should only be used once, later to be declared as a static value in the program.
 *
 * @warning
 * - Pointers shouldn't be NULL.
 * - `x` and `w` dimensions must be coherent.
 *
 */
void setParameters(Matrix *xCPP, Matrix *wCPP)
{
    // Must generate a copy because of R's gc() in c++

    Matrix *x = copyMatrixPtr(xCPP);
    Matrix *w = copyMatrixPtr(wCPP);

    // ---- Validation checks ---- //
    // ---- Check if there's a NULL pointer ----
    if (!x->data || !w->data)
    {
        error("Constructor: A NULL pointer was handed.\n");
    }

    // ---- Check for dimentional coherence ----
    if (x->cols != w->rows && x->cols > 0)
    {
        error("Constructor: The dimensions of the matrices handed are incorrect; `x` columns and `w` "
              "rows length "
              "must be the same, but they're %d and %d respectivately.\n",
              x->cols, w->rows);
    }

    // ---- Allocate memory for the global variables ---- //
    // ---- Since they're integers it will be better to operate with integers rather than a cBLAS operations (since it
    // receives doubles) ----
    TOTAL_CANDIDATES = x->rows;
    TOTAL_GROUPS = w->cols;
    TOTAL_BALLOTS = w->rows;
    CANDIDATES_VOTES = (uint32_t *)Calloc(TOTAL_CANDIDATES, uint32_t);
    GROUP_VOTES = (uint32_t *)Calloc(TOTAL_GROUPS, uint32_t);
    BALLOTS_VOTES = (uint16_t *)Calloc(TOTAL_BALLOTS, uint16_t);
    inv_BALLOTS_VOTES = (double *)Calloc(TOTAL_BALLOTS, double);

    // ---- Allocate memory for the matrices
    X = Calloc(1, Matrix);
    *X = createMatrix(x->rows, x->cols);
    memcpy(X->data, x->data, sizeof(double) * x->rows * x->cols);

    W = Calloc(1, Matrix);
    *W = createMatrix(w->rows, w->cols);
    memcpy(W->data, w->data, sizeof(double) * w->rows * w->cols);
    // ---...--- //

    // ---- Fill the variables ---- //
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate given a ballot box
            // ---- Add the candidate votes, ballot votes and total votes ----
            CANDIDATES_VOTES[c] += (uint32_t)MATRIX_AT_PTR(X, c, b);
            TOTAL_VOTES += (uint32_t)MATRIX_AT_PTR(
                X, c, b); // Usually, TOTAL_CANDIDATES < TOTAL_GROUPS, hence, it's better to make less sums.
            BALLOTS_VOTES[b] += (uint16_t)MATRIX_AT_PTR(X, c, b);
        } // --- End candidate loop

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            // ---- Add the group votes ----
            GROUP_VOTES[g] += (uint32_t)MATRIX_AT_PTR(W, b, g);
        } // --- End group loop

        // ---- Compute the inverse of the ballot votes, at this point the `b` votes are ready ----
        inv_BALLOTS_VOTES[b] = 1.0 / (double)BALLOTS_VOTES[b];
    } // --- End ballot box loop
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
Matrix getInitialP(const char *p_method)
{

    // ---- Validation: check the method input ----//
    if (strcmp(p_method, "uniform") != 0 && strcmp(p_method, "proportional") != 0 &&
        strcmp(p_method, "group_proportional") != 0)
    {
        error("Compute: The method `%s` to calculate the initial probability doesn't exist.\nThe supported methods "
              "are: `uniform`, `proportional` and `group_proportional`.\n",
              p_method);
    }
    // ---...--- //

    Matrix probabilities = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);

    // ---- Compute the uniform method ---- //
    // ---- It assumes a uniform distribution among candidates ----
    if (strcmp(p_method, "uniform") == 0)
    {
        fillMatrix(&probabilities, 1.0 / (double)TOTAL_CANDIDATES);
    }
    // ---...--- //

    // ---- Compute the proportional method ---- //
    // ---- It calculates the proportion of votes of each candidate, and assigns that probability to every demographic
    // group ----
    else if (strcmp(p_method, "proportional") == 0)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate
            double ratio =
                (double)CANDIDATES_VOTES[c] / (double)TOTAL_VOTES; // Proportion of candidates votes per total votes.
            for (int g = 0; g < TOTAL_GROUPS; g++)
            { // --- For each group, given a candidate
                MATRIX_AT(probabilities, g, c) = ratio;
            } // --- End group loop
        } // --- End candidate loop
    }
    // ---...--- //

    // ---- Compute the group_proportional method ---- //
    // ---- Considers the proportion of candidates votes and demographic groups aswell ----
    else
    {
        // ---- Create a temporary matrix to store the first results ----
        Matrix ballotProbability = createMatrix(TOTAL_BALLOTS, TOTAL_CANDIDATES);
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot vote
            double den = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate, given a ballot box
                MATRIX_AT(ballotProbability, b, c) = MATRIX_AT_PTR(X, c, b);
                den += MATRIX_AT(ballotProbability, b, c);
            }
            // ---- Handle border case ----
            if (den != 0)
            {
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // --- For each candidate, given a ballot box
                    MATRIX_AT(ballotProbability, b, c) /= den;
                }
            }
        } // --- End ballot box loop

        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot box
            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            { // --- For each group given a ballot box
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // --- For each candidate, given a ballot box and a group
                    MATRIX_AT(probabilities, g, c) += MATRIX_AT(ballotProbability, b, c) * MATRIX_AT_PTR(W, b, g);
                }
            }
        }

        // ---- Add the final values to the matrix
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group
                // ---- Handle border case ----
                if (GROUP_VOTES[g] == 0)
                    MATRIX_AT(probabilities, g, c) = 0;
                else
                    MATRIX_AT(probabilities, g, c) /= GROUP_VOTES[g];
            }
        }
        freeMatrix(&ballotProbability);
    }
    // ---...--- //
    return probabilities;
}

/*
 *
 * @brief Sets the configuration for the Q method, using a modularized approach.
 *
 * Given that different `Q` methods receive different parameters, a modularized approach is given towards each method
 *
 * @input[in] q_method A char with the q_method. Currently it supports "exact", "hnr", "mult", "mvn_cdf"
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
    else if (strcmp(q_method, "hnr") == 0)
    {
        config.computeQ = computeQHitAndRun;
    }
    else if (strcmp(q_method, "exact") == 0)
    {
        config.computeQ = computeQExact;
    }
    else if (strcmp(q_method, "mvn_cdf") == 0)
    {
        config.computeQ = computeQMultivariateCDF;
    }
    else if (strcmp(q_method, "mvn_pdf") == 0)
    {
        config.computeQ = computeQMultivariatePDF;
    }
    else
    {
        error("Compute: An invalid method was provided: `%s`\nThe supported methods are: `exact`, `hnr`"
              ", `mult`, `mvn_cdf` and `mvn_pdf`.\n",
              q_method);
    }

    // Directly store the input parameters
    config.params = inputParams;
    return config;
}

/*
 * @brief Computes the log-likelihood for a given probability and `q` array
 *
 * Given the conditional probability and the array, it computes the log-likelihood as defined in the paper.
 *
 * @param[in] q Array of matrices of dimension (bxgxc) that represents the probability that a voter of group "g" in
 * ballot box "b" voted for candidate "c" conditional on the observed result.
 * @param[in] prob A matrix of the calculated probabilities
 *
 * @return The value of the log-likelihood.
 *
 */

double logLikelihood(Matrix *prob, double *q)
{
    // ---- Define the summatory for the log-likelihood ---- //
    double logLL = 0;
    // ---- Outer summatory ----
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group, given a ballot box
            // ---- Inner summatory
            double cSummatory = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate, given a group and a ballot box
                double num = MATRIX_AT_PTR(prob, g, c);
                double den = Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
                if (den == 0)
                    continue;
                if (num == 0)
                    continue;
                double qval = den;
                cSummatory += qval * log(num / den);
            } // --- End c loop
            logLL += MATRIX_AT_PTR(W, b, g) * cSummatory;
        } // --- End g loop
    } // --- End b loop
    // ---...--- //
    return logLL;
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
Matrix getP(const double *q)
{
    // ---- Inititalize variables ---- //
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    // ---...--- //

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
                                 &W->data[g * TOTAL_BALLOTS], // Now correctly indexing W in column-major
                                 &newStride,                  // Column-major: stride is 1 for W
                                 baseY,                       // Column-major: index properly
                                 &stride                      // Stride: move down rows (1 step per row)
            );

            MATRIX_AT(toReturn, g, c) = val / GROUP_VOTES[g];
        }
    }
    // ---...--- //
    return toReturn;
}

/**
 * @brief Implements the whole EM algorithm.
 *
 * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
 * paper, it currently supports hnr, mult, mvn_cdf and mvn_pdf methods.
 *
 * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
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
Matrix EMAlgoritm(Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                  const double maxSeconds, const bool verbose, double *time, int *iterTotal, double *logLLarr,
                  double **qVal, int *finishing_reason, QMethodInput inputParams)
{
    // ---- Error handling is done on getQMethodConfig!
    if (verbose)
    {
        Rprintf("Starting the EM algorithm.\n");
        Rprintf("The method to calculate the conditional probability will be %s method with the following "
                "parameters:\nConvergence threshold:\t%.6f\nMaximum iterations:\t%d\n",
                q_method, convergence, maxIter);
    }

    // ---- Define the parameters for the main loop ---- //
    QMethodConfig config = getQMethodConfig(q_method, inputParams);
    double newLL;
    double oldLL = -999999;
    /*
    *qVal = config.computeQ(currentP, config.params); // Calcula el Q
    double oldLL = logLikelihood(currentP, *qVal); // Loglikelihood con P INICIAL
    double newLL;
    */
    // ---- Start timer
    struct timespec start, end, iter_start, iter_end; // Declare timers for overall and per-iteration
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double elapsed_total = 0;
    Matrix newProbability;
    // ---...--- //
    // ---- Execute the EM-iterations ---- //
    for (int i = 0; i < maxIter; i++)
    {
        // ---- Timer for the current iteration
        clock_gettime(CLOCK_MONOTONIC, &iter_start);
        *iterTotal = i;
        if (verbose)
        {
            Rprintf("\nThe current probability matrix at the %dth iteration is:\n", i);
            printMatrix(currentP);
        }

        // ---- Compute the `q`, `p` and `log-likelihood` ---- //
        Free(*qVal);
        *qVal = config.computeQ(currentP, config.params);
        newLL = logLikelihood(currentP, *qVal);
        // logLikelyArray[i] = newLL;
        // CONDICION DE PARADA ACÁ
        newProbability = copyMatrix(currentP);
        freeMatrix(currentP);
        *currentP = getP(*qVal); // MSTEP
        // ---...--- //

        // ---- Check convergence ---- //
        // Usually the CDF do really SMALL steps, so at least impose some iterations...
        int minIter = (strcmp(q_method, "mvn_cdf") == 0) ? 65 : 0;
        // if (fabs(newLL - oldLL) < convergence && i >= minIter)
        // if (false)
        if (convergeMatrix(&newProbability, currentP, convergence))
        {
            // ---- End timer ----
            clock_gettime(CLOCK_MONOTONIC_RAW, &end);

            if (verbose)
            {
                Rprintf("The convergence was found on iteration %d, with a log-likelihood of %.4f  and took %.5f "
                        "seconds!\n",
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
            Rprintf("Iteration %d took %.5f seconds.\n", i + 1, elapsed_iter);

        // ---- The maximum time was reached
        if (elapsed_total >= maxSeconds)
        {
            if (verbose)
                Rprintf("Maximum time limit reached.\n");
            *finishing_reason = 1;
            goto results;
        }
        oldLL = newLL;
    }
    Rprintf("Maximum iterations reached without convergence.\n"); // Print even if there's not verbose, might change
                                                                  // later.
    *finishing_reason = 2;
results:
    *logLLarr = newLL;
    *time = elapsed_total;
    // ---- Matrix must be returned without a pointer
    Matrix finalProbability = copyMatrix(currentP);
    freeMatrix(currentP);
    freeMatrix(&newProbability);
    return finalProbability;
    // return *currentP;
}

// ---- Clean all of the global variables ---- //
// __attribute__((destructor)) // Executes when the library is ready
void cleanup()
{
    TOTAL_VOTES = 0;
    TOTAL_BALLOTS = 0;
    TOTAL_CANDIDATES = 0;
    TOTAL_GROUPS = 0;

    if (CANDIDATES_VOTES != NULL)
    {
        Free(CANDIDATES_VOTES);
        CANDIDATES_VOTES = NULL;
    }
    if (GROUP_VOTES != NULL)
    {
        Free(GROUP_VOTES);
        GROUP_VOTES = NULL;
    }
    if (BALLOTS_VOTES != NULL)
    {
        Free(BALLOTS_VOTES);
        BALLOTS_VOTES = NULL;
    }
    if (inv_BALLOTS_VOTES != NULL)
    {
        Free(inv_BALLOTS_VOTES);
        inv_BALLOTS_VOTES = NULL;
    }

    if (X->data != NULL) // Note that the columns and rows are usually stack.
    {
        freeMatrix(X);
        Free(X);
        X = NULL;
    }
    if (W->data != NULL)
    {
        freeMatrix(W);
        Free(W);
        W = NULL;
    }
}
