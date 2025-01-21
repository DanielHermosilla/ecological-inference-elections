#include "exact.h"
#include "fileutils.h"
#include "globals.h"
#include "hitAndRun.h"
#include "instanceGenerator.h"
#include "matrixUtils.h"
#include "multinomial.h"
#include "multivariate-cdf.h"
#include "multivariate-pdf.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <cblas.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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
void setParameters(Matrix *x, Matrix *w)
{

    // ---- Validation checks ---- //
    // ---- Check if there's a NULL pointer ----
    if (!x->data || !w->data)
    {
        fprintf(stderr, "setParameters: A NULL pointer was handed.\n");
        exit(EXIT_FAILURE); // Frees used memory
    }

    // ---- Check for dimentional coherence ----
    if (x->cols != w->rows && x->cols > 0)
    {
        fprintf(stderr,
                "setParameters: The dimensions of the matrices handed to countVotes are incorrect; `x` columns and `w` "
                "rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x->cols, w->rows);
        exit(EXIT_FAILURE);
    }

    // ---- Allocate memory for the global variables ---- //
    // ---- Since they're integers it will be better to operate with integers rather than a cBLAS operations (since it
    // receives doubles) ----
    TOTAL_CANDIDATES = x->rows;
    TOTAL_GROUPS = w->cols;
    TOTAL_BALLOTS = w->rows;
    CANDIDATES_VOTES = (uint32_t *)calloc(TOTAL_CANDIDATES, sizeof(uint32_t));
    GROUP_VOTES = (uint32_t *)calloc(TOTAL_GROUPS, sizeof(uint32_t));
    BALLOTS_VOTES = (uint16_t *)calloc(TOTAL_BALLOTS, sizeof(uint16_t));
    inv_BALLOTS_VOTES = (double *)calloc(TOTAL_BALLOTS, sizeof(double));

    // ---- Allocate memory for the matrices
    X = malloc(sizeof(Matrix));
    *X = createMatrix(x->rows, x->cols);
    memcpy(X->data, x->data, sizeof(double) * x->rows * x->cols);

    W = malloc(sizeof(Matrix));
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
        strcmp(p_method, "group proportional") != 0)
    {
        fprintf(stderr, "getInitialP: The method `%s` passed to getInitialP doesn't exist.\n", p_method);
        exit(EXIT_FAILURE);
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

    // ---- Compute the group proportional method ---- //
    // ---- Considers the proportion of candidates votes and demographic groups aswell ----
    else
    {
        double numerator = 0.0; // Note that the denominator is already known
        uint16_t temp = 0;      // Will be used to cast multiplications on integers and add efficiency.

        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot vote
            if (BALLOTS_VOTES[b] == 0)
                continue; // Division by zero, even though it's very unlikely (a ballot doesn't have any vote).

            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            { // --- For each group given a ballot box
                if (MATRIX_AT_PTR(W, b, g) == 0.0)
                    continue; // Division by zero, case where a group doesn't vote in a ballot.
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                { // --- For each candidate given a group and a ballot box.
                    temp = (uint16_t)MATRIX_AT_PTR(X, c, b) * (uint16_t)MATRIX_AT_PTR(W, b, g); // w_bg * x_bg = a
                    numerator = temp * inv_BALLOTS_VOTES[b];                                    // (a/I_b)
                    MATRIX_AT(probabilities, g, c) += numerator;
                }
            }
        }

        // ---- Now do the division once per (g,c), this reduces roughly 2 millions divisions to 50 ----
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group
            // ---- Error handling ----
            if (GROUP_VOTES[g] == 0)
            {
                fprintf(stderr, "getInitialP: The %dth group does not have any vote assigned.", g);
                exit(EXIT_FAILURE);
            }

            double inv_gvotes = 1.0 / GROUP_VOTES[g];
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group
                MATRIX_AT(probabilities, g, c) *= inv_gvotes;
            }
        }
    }
    // ---...--- //
    return probabilities;
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
    double *ptrReturn = toReturn.data; // For OpenMP because they don't accept structs for reduction clauses
                                       // even though it's dereferenced

    // ---- Compute the dot products ---- //
    // ---- Due to the dot product, the parallelization is worth it ----
#pragma omp parallel for collapse(2) schedule(static)
    for (int g = 0; g < TOTAL_GROUPS; g++)
    { // --- For each group
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate given a group
            // Dot product over b=0..B-1 of W_{b,g} * Q_{b,g,c}
            double val = cblas_ddot(TOTAL_BALLOTS,
                                    &W->data[g],                    // points to W_{0,g}
                                    TOTAL_GROUPS,                   // stride: each next W_{b+1,g} is +G in memory
                                    &q[g * TOTAL_CANDIDATES + c],   // points to Q_{0,g,c}
                                    TOTAL_GROUPS * TOTAL_CANDIDATES // stride: each next Q_{b+1,g,c} is +(G*C) in memory
            );
            ptrReturn[g * TOTAL_CANDIDATES + c] = val; // Equivalent to MATRIX_AT(probabilities, g, c)
        }
    }

    // ---...--- //

    // ---- Perform the divisions apart ---- //
    // ---- Now divide by GROUP_VOTES[g] just once per (g,c) instead of doing it `b` times (note that the division is
    // usually expensive). Approximately reduces 400.000 double divisions to 50. ----
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // --- For each group
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate given a group
            ptrReturn[g * TOTAL_CANDIDATES + c] /= GROUP_VOTES[g];
        }
    }
    // ---...--- //

    return toReturn;
}

/**
 * @brief Implements the whole EM algorithm.
 *
 * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
 * paper, it currently supports Hit and Run, Multinomial, MVN CDF and MVN PDF methods.
 *
 * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "Multinomial", "MVN CDF", "MVN PDF" and "Exact" methods.
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
 * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
 * @param[in] verbose Wether to verbose useful outputs.
 *
 * @return Matrix: A matrix with the final probabilities. In case it doesn't converges, it returns the last
 * probability that was computed
 *
 * @note This is the main function that calls every other function for "q"
 *
 * @see getInitialP() for getting initial probabilities. Group proportional method is recommended.
 *
 * @warning
 * - Pointers shouldn't be NULL.
 * - `x` and `w` dimensions must be coherent.
 *
 */
Matrix EMAlgoritm(Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                  const bool verbose)
{

    if (strcmp(q_method, "Hit and Run") != 0 && strcmp(q_method, "Multinomial") != 0 &&
        strcmp(q_method, "MVN CDF") != 0 && strcmp(q_method, "MVN PDF") != 0 && strcmp(q_method, "Exact") != 0)
    {
        fprintf(stderr, "The method `%s` passed to EMAlgorithm doesn't exist.\n", q_method);
        exit(EXIT_FAILURE);
    }

    if (verbose)
    {
        printf("Starting the EM algorithm.\n The candidates matrix is:\n");
        printMatrix(X);
        printf("\nThe matrix with the demographic groups votation is:\n");
        printMatrix(W);
        printf("\nThe method to calculate the conditional probability will be %s method with the following "
               "parameters:\nConvergence threshold:\t%.6f\nMaximum iterations:\t%d\n",
               q_method, convergence, maxIter);
    }

    double *q;

    struct timespec start, end; // Start time

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < maxIter; i++)
    {
        if (verbose)
        {
            printf("\n%d of iterations have been done.\nThe current probability matrix is:\n", i);
            printMatrix(currentP);
        }

        // Exact method
        if (strcmp(q_method, "Exact") == 0)
        {
            q = computeQExact(currentP);
        }

        else if (strcmp(q_method, "Hit and Run") == 0)
        {
            q = computeQHitAndRun(currentP, 10000, 1000);
        }
        // Multinomial
        else if (strcmp(q_method, "Multinomial") == 0)
        {
            q = computeQMultinomial(currentP);
        }
        // MVN CDF
        else if (strcmp(q_method, "MVN CDF") == 0)
        {
            q = computeQMultivariateCDF(currentP, 10000, 0.00001, "Genz");
        }
        // MVN PDF
        else
        {
            q = computeQMultivariatePDF(currentP);
        }

        Matrix newProbability = getP(q);
        free(q);
        if (convergeMatrix(&newProbability, currentP, convergence))
        {

            // End timer
            clock_gettime(CLOCK_MONOTONIC, &end);
            double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

            if (verbose)
            {
                printf("The convergence was found on iteration %d and took %.5f seconds!\n", i, elapsed);
            }
            freeMatrix(currentP);
            return newProbability;
        }

        freeMatrix(currentP);
        *currentP = createMatrix(newProbability.rows, newProbability.cols);
        memcpy(currentP->data, newProbability.data, sizeof(double) * newProbability.rows * newProbability.cols);
        freeMatrix(&newProbability);
    }
    printf("Maximum iterations reached without convergence.\n"); // Print even if there's not verbose, might change
    // later.
    return *currentP;
}

void testProb()
{

    Matrix prob = getInitialP("uniform");
    printf("\nThe probability matrix for `uniform` method is:\n");
    printMatrix(&prob);
    freeMatrix(&prob);

    Matrix prob2 = getInitialP("proportional");
    printf("\nThe probability matrix for `proportional` method is:\n");
    printMatrix(&prob2);
    freeMatrix(&prob2);

    Matrix prob3 = getInitialP("group proportional");
    printf("\nThe probability matrix for `group proportional` method is:\n");
    printMatrix(&prob3);
    freeMatrix(&prob3);
}

__attribute__((destructor)) // Executes when the library is ready
void cleanup()
{
    if (CANDIDATES_VOTES != NULL)
    {
        free(CANDIDATES_VOTES);
        CANDIDATES_VOTES = NULL;
    }
    if (GROUP_VOTES != NULL)
    {
        free(GROUP_VOTES);
        GROUP_VOTES = NULL;
    }
    if (BALLOTS_VOTES != NULL)
    {
        free(BALLOTS_VOTES);
        BALLOTS_VOTES = NULL;
    }
    if (inv_BALLOTS_VOTES != NULL)
    {
        free(inv_BALLOTS_VOTES);
        inv_BALLOTS_VOTES = NULL;
    }

    if (X->data != NULL) // Note that the columns and rows are usually stack.
    {
        freeMatrix(X);
        free(X);
        X = NULL;
    }
    if (W->data != NULL)
    {
        freeMatrix(W);
        free(W);
        W = NULL;
    }
}

int main()
{
    printf("The program is running\n");
    /*
    Matrix XX = {.data = NULL, .rows = 0, .cols = 0};
    Matrix G = {.data = NULL, .rows = 0, .cols = 0};
    char *method = "multinomial";
    // createInstance(&XX, &G, 42, *method); // TODO: Arreglar esto para poder crear una instancia...
    //Matrix matrices[2] = {XX, G};
    //writeMatrices("matricesTest5.bin", matrices, 2);
*/
    /*
     */
    Matrix matrixArray[2];
    readMatrices("matricesTest5.bin", matrixArray, 2);
    Matrix XX = matrixArray[0];
    Matrix G = matrixArray[1];

    // Start timer
    struct timespec start, end; // Start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    setParameters(&XX, &G);
    Matrix P = getInitialP("group proportional");

    Matrix Pnew = EMAlgoritm(&P, "Hit and Run", 0.000001, 10000, true);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("The whole algorithm calculation took %.16f seconds!\n", elapsed);

    printMatrix(&Pnew);
    freeMatrix(&Pnew);
    freeMatrix(&XX);
    freeMatrix(&G);
    free(XX.data);
    free(G.data);

    return 1;
}
