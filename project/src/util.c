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

static uint32_t TOTAL_VOTES = 0;
static uint32_t TOTAL_BALLOTS = 0;
static uint16_t TOTAL_CANDIDATES = 0;
static uint16_t TOTAL_GROUPS = 0;
static uint32_t *CANDIDATES_VOTES = NULL;
static uint32_t *GROUP_VOTES = NULL;
static Matrix *X = NULL;
static Matrix *W = NULL;

void setParameters(Matrix x, Matrix w)
{
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
                "The dimensions of the matrices handed to countVotes are incorrect; `x` columns and `w` rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x.cols, w.rows);
        exit(EXIT_FAILURE);
    }

    // Since they're integers it will be better to operate with integers rather than a cBLAS operations (since it
    // receives doubles).
    TOTAL_CANDIDATES = x.rows;
    TOTAL_GROUPS = w.cols;
    TOTAL_BALLOTS = w.rows;
    uint32_t *CANDIDATES_VOTES = (uint32_t *)malloc(TOTAL_CANDIDATES * sizeof(uint32_t));
    uint32_t *GROUP_VOTES = (uint32_t *)malloc(TOTAL_GROUPS * sizeof(uint32_t));
    *X = x;
    *W = w;

#pragma omp parallel for reduction(+ : CANDIDATES_VOTES[ : TOTAL_CANDIDATES])                                          \
    reduction(+ : GROUP_VOTES[ : TOTAL_GROUPS]) reduction(+ : TOTAL_VOTES)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            CANDIDATES_VOTES[c] += (uint32_t)MATRIX_AT(x, c, b);
            TOTAL_VOTES += (uint32_t)MATRIX_AT(
                x, c, b); // Usually it's always the candidate with lesser dimension, preferible to not
                          // compromise legibility over a really small and unprobable gain in efficiency
        }
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            GROUP_VOTES[g] += (uint32_t)MATRIX_AT(w, b, g);
        }
    }
}

Matrix getInitialP(const char *p_method)
{

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

    // Validation, checks the method string
    if (strcmp(p_method, "uniform") != 0 && strcmp(p_method, "proportional") != 0 &&
        strcmp(p_method, "group proportional") != 0)
    {
        fprintf(stderr, "The method `%s` passed to getInitialP doesn't exist.\n", p_method);
        exit(EXIT_FAILURE);
    }

    Matrix probabilities;
    probabilities = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);

    // Asumes a Uniform distribution among candidates
    if (strcmp(p_method, "uniform") == 0)
    {
        fillMatrix(&probabilities, 1.0 / (double)TOTAL_CANDIDATES);
    }
    // The proportional method calculates the proportion of votes of each candidate, and assigns that probability to
    // every demographic group;

    else if (strcmp(p_method, "proportional") == 0)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double ratio =
                (double)CANDIDATES_VOTES[c] / (double)TOTAL_VOTES; // Proportion of candidates votes per total votes.
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                MATRIX_AT(probabilities, g, c) = ratio;
            }
        }
    }

    else // group proportional
         // Now it considers the proportion of candidates votes AND demographic group, it's an extension of the
         // past method
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                MATRIX_AT(probabilities, g,
                          c) = // This considers a joint proportion; TODO: correct considering the paper equation.
                    (double)(CANDIDATES_VOTES[c] * GROUP_VOTES[g]) / (double)(TOTAL_CANDIDATES * TOTAL_GROUPS);
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
        CANDIDATES_VOTES = NULL;
    }
}

int main()
{
    printf("The program is running\n");
    time_t start, end;
    start = clock();
    Matrix X = createMatrix(3, 3);
    Matrix G = createMatrix(3, 2);
    double xVal[9] = {0, 0, 2, 4, 5, 8, 1, 2, 3};
    double gVal[6] = {2, 5, 3, 1, 0, 14};
    memcpy(X.data, xVal, sizeof(xVal));
    memcpy(G.data, gVal, sizeof(gVal));

    const int32_t TOTAL_VOTES = countVotes(X, G);
    printf("\nTotal votes are %d\n\n", TOTAL_VOTES);

    testProb(X, G, TOTAL_VOTES);
    freeMatrix(&X);
    freeMatrix(&G);
    end = clock();
    int t = (end - start) / CLOCKS_PER_SEC;
    printf("The program took %d seconds!", t);
    return 0;
}
