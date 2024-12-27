#include "matrixUtils.h"
#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <cblas.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// note: Refer to https://www.seehuhn.de/pages/linear.html#blas
// double m[] = {
// 3, 1, 3,
// 1, 5, 9,
// 2, 6, 5
// };

// note: THIS IS THE STRUCT DEFINED IN THE matrixUtils.h file, this comment is just for clarity
// typedef struct
//{
//  double *data; // Pointer to matrix data array (row-major order)
//  int rows;     // Number of rows
//  int cols;     // Number of columns
//} Matrix;

// Macro for making an easier indexation.
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])
// Macro for evaluating maximum numbers, be aware of possible issues, try to use it on ints only.
#define MAX(a, b) ((a) > (b) ? a : b)
// Macro for evaluating minimum numbers, be aware of possible issues, try to use it on ints only.
#define MIN(a, b) ((a) < (b) ? a : b)

int32_t countVotes(Matrix x, Matrix w)
{
    /**
     * @brief Yields the total amount of votes in the process. Usually this should be done once for avoiding
     * computing a loop over ballots.
     *
     * Gets the total amount of votes in the process. This should only be donce once for avoiding computing loops
     * over ballots.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     *
     * @return An integer with the total amount of votes of the process.
     *
     * @note This should only be used once, later to be declared as a constant value in the program.
     *
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */

    // Validation, checks NULL pointer
    if (!x.data || !w.data)
    {
        fprintf(stderr, "A NULL pointer was handed to getInitialP.\n");
        exit(EXIT_FAILURE); // Frees used memory
    }

    // Validation, checks dimentional coherence
    if (x.cols != w.rows && x.cols > 0)
    {
        fprintf(stderr,
                "The dimensions of the matrices handed to getInitialP are incorrect; `x` columns and `w` rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x.cols, w.rows);
        exit(EXIT_FAILURE);
    }

    // Since they're integers it will be better to operate with integers rather than a cBLAS operations (since it
    // receives doubles).
    const int size1 = x.rows * x.cols;
    const int size2 = w.rows * w.cols;

    int smallerDimension = MIN(size1, size2);
    int32_t total = 0;

    if (smallerDimension == size1)
    { // Candidates have a lesser dimension
#pragma omp parallel for reduction(+ : total)
        for (int i = 0; i < size1; i++)
        {
            total += (int32_t)(x.data[i]);
        }
    }

    else
    { // Groups have a lesser dimension
#pragma omp parallel for reduction(+ : total)
        for (int i = 0; i < size2; i++)
        {
            total += (int32_t)(w.data[i]);
        }
    }
    return total;
}

Matrix getInitialP(Matrix x, Matrix w, const char *p_method, const int32_t totalVotes)
{

    /**
     * @brief Computes the initial probability of the EM algoritm.
     *
     * Given the observables results, it computes a convenient initial "p" value for initiating the
     * algorithm. Currently it supports the "uniform", "group_proportional" and "proportional" methods.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     * @param[in] p_method The method for calculating the initial parameter. Currently it supports "uniform",
     * "group_proportional" and "proportional" methods.
     * @param[in] totalVotes An 32-bit integer (assuming the votes doesn't exceed 2.1billions...) with the total amount
     * of votes. This parameter is used for computing only once the total amount of votes and being more efficient.
     * @return Matrix of dimension (gxc) with the initial probability for each demographic group "g" voting for a given
     * candidate "c".
     * @note This should be used only that the first iteration of the EM-algorithm.
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */

    // Validation, checks NULL pointer
    if (!x.data || !w.data)
    {
        fprintf(stderr, "A NULL pointer was handed to getInitialP.\n");
        exit(EXIT_FAILURE); // Frees used memory
    }

    // Validation, checks dimentional coherence
    if (x.cols != w.rows && x.cols > 0)
    {
        fprintf(stderr,
                "The dimensions of the matrices handed to getInitialP are incorrect; `x` columns and `w` rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x.cols, w.rows);
        exit(EXIT_FAILURE);
    }

    // Validation, checks the method string
    if (strcmp(p_method, "uniform") != 0 && strcmp(p_method, "proportional") != 0 &&
        strcmp(p_method, "group proportional") != 0)
    {
        fprintf(stderr, "The method `%s` passed to getInitialP doesn't exist.\n", p_method);
        exit(EXIT_FAILURE);
    }

    const int candidates = x.rows;
    const int groups = w.cols;
    const int ballots = x.cols;

    Matrix probabilities;
    probabilities = createMatrix(groups, candidates);

    // Asumes a Uniform distribution among candidates
    if (strcmp(p_method, "uniform") == 0)
    {
        fillMatrix(&probabilities, 1.0 / candidates);
    }
    // The proportional method calculates the proportion of votes of each candidate, and assigns that probability to
    // every demographic group;

    else if (strcmp(p_method, "proportional") == 0)
    {
        double canVotes[candidates];
        // Step 1: Calculate the total amount of votes per candidates.
        rowSum(&x, canVotes);

        // Step 2: Fill the matrix
        for (int c = 0; c < candidates; c++)
        {
            double ratio = (double)canVotes[c] / (double)totalVotes;
            for (int g = 0; g < groups; g++)
            {
                MATRIX_AT(probabilities, g, c) = ratio;
            }
        }
    }

    else // group proportional
         // Now it considers the proportion of candidates votes AND demographic group, it's an extension of the
         // past method
    {
        double votePerCandidate[candidates];
        int totalCandidate = 0;
        rowSum(&x, votePerCandidate);

        double votePerGroup[groups];
        int totalGroup = 0;
        colSum(&w, votePerGroup);

        // Step 1: Calculate the total amount of votes for candidate
        for (int c = 0; c < candidates; c++)
        {
            totalCandidate += votePerCandidate[c];
        }

        // Step 2: Calculate the total amount of votes per group
        for (int g = 0; g < groups; g++)
        {
            totalGroup += votePerGroup[g];
        }

        // Step 3: Fill the matrix
        for (int c = 0; c < candidates; c++)
        {
            for (int g = 0; g < groups; g++)
            {
                MATRIX_AT(probabilities, g, c) =
                    (votePerCandidate[c] * votePerGroup[g]) / (double)(totalCandidate * totalGroup);
            }
        }
    }
    return probabilities;
}

Matrix getP(Matrix w, double q, int32_t totalVotes)

/*
 * @brief Computes the optimal solution for the `M` step
 *
 * Given the conditional probability and the votations per demographic group, it calculates the new probability for
 * the next iteration.
 *
 * @param[in] q Array of matrices of dimension (bxgxc) that represents the probability that a voter of group "g" in
 * ballot box "b" voted for candidate "c" conditional on the observed result.
 * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
 *
 * @return A matrix with the optimal probabilities according the Log-likelihood.
 *
 * @see getInitialP() for getting initial probabilities. This method is recommended to be used exclusively for the EM
 * Algorithm, unless there's a starting "q" to start with.
 *
 */

{
    // Matrix returnMatrix = createMatrix(w.cols, q[0].cols);
    return w;
}

Matrix EMAlgoritm(Matrix *x, Matrix *w, Matrix *currentP, const char *q_method, const int32_t totalVotes,
                  const double convergence, const int maxIter, const bool verbose)
// NOTE: Use pointer on matrices to avoid copying huge matrices on stack memory.
{

    /**
     * @brief Implements the whole EM algorithm.
     *
     * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
     * paper, it currently supports Hit and Run, Multinomial, MVN CDF and MVN PDF methods.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "Multinomial", "MVN CDF" and "MVN PDF" methods.
     * @param[in] totalVotes An 32 bit sized integer that has the total amount of votes. It's there to simplify the
     * M-step.
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

    if (strcmp(q_method, "Hit and Run") != 0 && strcmp(q_method, "Multinomial") != 0 &&
        strcmp(q_method, "MVN CDF") != 0 && strcmp(q_method, "MVN PDF") != 0)
    {
        fprintf(stderr, "The method `%s` passed to EMAlgorithm doesn't exist.\n", q_method);
        exit(EXIT_FAILURE);
    }

    if (verbose)
    {
        printf("Starting the EM algorithm.\n The candidates matrix is:\n");
        printMatrix(x);
        printf("\nThe matrix with the demographic groups votation is:\n");
        printMatrix(w);
        printf("\nThe method to calculate the conditional probability will be %s method with the following "
               "parameters:\nConvergence threshold:\t%.6f\nMaximum iterations:\t%d\n",
               q_method, convergence, maxIter);
    }

    // There are some things that can be calculated ONCE and reused on the E-step:
    // Total amount of votes (i.e \sum_{b\in B, g\in G}w_{bg}). It can either be done with
    // the candidates vector or group vector. TODO: Reuse this value also for initial probability

    // Matrix *probabilityPtr = initialP;
    double q = 2.0;

    for (int i = 0; i < maxIter; i++)
    {
        // Hit and Run
        if (strcmp(q_method, "Hit and Run") == 0)
        {
            printf("Executing 'Hit and Run' method.\n");
            break;
        }
        // Multinomial
        else if (strcmp(q_method, "Multinomial") == 0)
        {
            printf("Executing 'Multinomial' method.\n");
            break;
        }
        // MVN CDF
        else if (strcmp(q_method, "MVN CDF") == 0)
        {
            printf("Executing 'MVN CDF' method.\n");
            break;
        }
        // MVN PDF
        else
        {
            printf("Executing 'MVN PDF' method.\n");
            break;
        }

        Matrix newProbability = getP(*currentP, q, totalVotes);

        if (convergeMatrix(&newProbability, currentP, convergence))
        {
            if (verbose)
            {
                printf("The convergence was found on iteration %d\n", i);
            }
            freeMatrix(currentP);
            freeMatrix(x);
            freeMatrix(w);
            return newProbability;

            freeMatrix(currentP);
            *currentP = newProbability;
        }
        // TODO: Compute de Log-likelihood and all the other stuff
    }
    printf("Maximum iterations reached without convergence.\n"); // Print even if there's not verbose, might change
    // later.
    return *currentP;
}
// SEXP hello_gsl()
//{
//     printf("Hello, World from C!\n");

//   return R_NilValue;
//
void testProb(Matrix X, Matrix G, int32_t votes)
{

    printf("Running test for the initial probability matrix\n The `X` matrix with the candidates votes is:\n");
    printMatrix(&X);
    printf("\nThe `w` matrix with the groups votes is:\n");
    printMatrix(&G);

    Matrix prob = getInitialP(X, G, "uniform", votes);
    printf("\nThe probability matrix for `uniform` method is:\n");
    printMatrix(&prob);
    freeMatrix(&prob);

    Matrix prob2 = getInitialP(X, G, "proportional", votes);
    printf("\nThe probability matrix for `proportional` method is:\n");
    printMatrix(&prob2);
    freeMatrix(&prob2);

    Matrix prob3 = getInitialP(X, G, "group proportional", votes);
    printf("\nThe probability matrix for `group proportional` method is:\n");
    printMatrix(&prob3);
    freeMatrix(&prob3);
}

int main()
{
    printf("The program is running\n");
    time_t start, end;
    start = clock();
    Matrix X = createMatrix(3, 3);
    Matrix G = createMatrix(3, 2);
    double xVal[9] = {0, 0, 2, 4, 5, 8, 1, 2, 3};
    double gVal[6] = {2, 5, 3, 1, 0, 9};
    const int32_t TOTAL_VOTES = countVotes(X, G);

    memcpy(X.data, xVal, sizeof(xVal));
    memcpy(G.data, gVal, sizeof(gVal));

    testProb(X, G, TOTAL_VOTES);
    freeMatrix(&X);
    freeMatrix(&G);
    end = clock();
    int t = (end - start) / CLOCKS_PER_SEC;
    printf("The program took %d seconds!", t);
    return 0;
}
