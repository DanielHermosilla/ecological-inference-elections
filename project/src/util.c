#include "fileutils.h"
#include "instanceGenerator.h"
#include "matrixUtils.h"
#include "multinomial.h"
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

// note: THIS IS THE STRUCT DEFINED IN THE matrixUtils.h file, this comment is just for clarity
// typedef struct
//{
//  double *data; // Pointer to matrix data array (row-major order)
//  int rows;     // Number of rows
//  int cols;     // Number of columns
//} Matrix;
//

// Macro for accessing a 3D flattened array (b x g x c)
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])

// Macro for making an easier indexation.
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])
#define MATRIX_AT_PTR(matrix, i, j) (matrix->data[(i) * (matrix->cols) + (j)])

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

void setParameters(Matrix *x, Matrix *w)
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
    if (!x->data || !w->data)
    {
        fprintf(stderr, "A NULL pointer was handed to getInitialP.\n");
        exit(EXIT_FAILURE); // Frees used memory
    }

    // Validation, checks dimentional coherence
    if (x->cols != w->rows && x->cols > 0)
    {
        fprintf(stderr,
                "The dimensions of the matrices handed to countVotes are incorrect; `x` columns and `w` rows length "
                "must be the same, but they're %d and %d respectivately.\n",
                x->cols, w->rows);
        exit(EXIT_FAILURE);
    }

    // Since they're integers it will be better to operate with integers rather than a cBLAS operations (since it
    // receives doubles).
    TOTAL_CANDIDATES = x->rows;
    TOTAL_GROUPS = w->cols;
    TOTAL_BALLOTS = w->rows;
    CANDIDATES_VOTES = (uint32_t *)calloc(TOTAL_CANDIDATES, sizeof(uint32_t));
    GROUP_VOTES = (uint32_t *)calloc(TOTAL_GROUPS, sizeof(uint32_t));
    BALLOTS_VOTES = (uint16_t *)calloc(TOTAL_BALLOTS, sizeof(uint16_t));
    inv_BALLOTS_VOTES = (double *)calloc(TOTAL_BALLOTS, sizeof(double));

    X = x;
    W = w;

    //#pragma omp parallel for reduction(+ : CANDIDATES_VOTES[ : TOTAL_CANDIDATES])                                          \
    reduction(+ : GROUP_VOTES[ : TOTAL_GROUPS]) reduction(+ : TOTAL_VOTES)                                             \
    reduction(+ : BALLOTS_VOTES[ : TOTAL_BALLOTS])
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            CANDIDATES_VOTES[c] += (uint32_t)MATRIX_AT_PTR(X, c, b);
            TOTAL_VOTES += (uint32_t)MATRIX_AT_PTR(
                X, c, b); // Usually it's always the candidate with lesser dimension, preferible to not
                          // compromise legibility over a really small and unprobable gain in efficiency
            BALLOTS_VOTES[b] += (uint16_t)MATRIX_AT_PTR(X, c, b);
        }
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            GROUP_VOTES[g] += (uint32_t)MATRIX_AT_PTR(W, b, g);
        }
    }

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        inv_BALLOTS_VOTES[b] = 1.0 / (double)BALLOTS_VOTES[b];
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
         // Now it considers the proportion of candidates votes AND demographic groups, it's an extension of the
         // past method
    {
        double numerator = 0.0;                 // Note that the denominator is already known
        uint32_t temp = 0;                      // Will be used to cast multiplications on integers and add efficiency.
        double *prob_data = probabilities.data; // For OpenMP because they don't accept structs for reduction clauses
                                                // even though it's dereferenced
#pragma omp parallel for reduction(+ : prob_data[ : TOTAL_GROUPS * TOTAL_CANDIDATES])
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            if (BALLOTS_VOTES[b] == 0)
                continue; // Division by zero, even though it's very unlikely (a ballot doesn't have any vote).

            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            {
                if (MATRIX_AT_PTR(W, b, g) == 0.0)
                    continue; // Division by zero, case where a group doesn't vote in a ballot.
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    temp = (uint32_t)MATRIX_AT_PTR(X, c, b) * (uint32_t)MATRIX_AT_PTR(W, b, g); // w_bg * x_bg = a
                    numerator = (double)temp * inv_BALLOTS_VOTES[b];                            // (a/I_b)
                    prob_data[g * TOTAL_CANDIDATES + c] += numerator; // Equivalent to MATRIX_AT(probabilities, g, c)
                }
            }
        }

        // Now do the division once per (g,c), this reduces roughly 2 millions divisions to 50
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            if (GROUP_VOTES[g] == 0)
            {
                fprintf(stderr, "The %dth group does not have any vote assigned.", g);
                exit(EXIT_FAILURE);
            }

            double inv_gvotes = 1.0 / (double)GROUP_VOTES[g];
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                prob_data[g * TOTAL_CANDIDATES + c] *= inv_gvotes;
            }
        }
    }
    return probabilities;
}

Matrix getP(const double *q, const bool continuous)

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

{
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *ptrReturn = toReturn.data; // For OpenMP because they don't accept structs for reduction clauses
                                       // even though it's dereferenced

#pragma omp parallel for collapse(2)
    for (int g = 0; g < TOTAL_GROUPS; g++)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {

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

    // Now divide by GROUP_VOTES[g] just once per (g,c) instead of doing it `b` times (note that the division is
    // usually expensive). Approximately reduces 400.000 double divisions to 50.
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            ptrReturn[g * TOTAL_CANDIDATES + c] /= GROUP_VOTES[g];
        }
    }
    return toReturn;
}

Matrix EMAlgoritm(Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                  const bool verbose)
// NOTE: Use pointer on matrices to avoid copying huge matrices on stack memory.
{

    /**
     * @brief Implements the whole EM algorithm.
     *
     * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
     * paper, it currently supports Hit and Run, Multinomial, MVN CDF and MVN PDF methods.
     *
     * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "Multinomial", "MVN CDF" and "MVN PDF" methods.
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
        printMatrix(X);
        printf("\nThe matrix with the demographic groups votation is:\n");
        printMatrix(W);
        printf("\nThe method to calculate the conditional probability will be %s method with the following "
               "parameters:\nConvergence threshold:\t%.6f\nMaximum iterations:\t%d\n",
               q_method, convergence, maxIter);
    }

    // There are some things that can be calculated ONCE and reused on the E-step:
    double *q;

    struct timespec start, end; // Start time

    // Start timer
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (int i = 0; i < maxIter; i++)
    {
        if (i % 10 == 0 && verbose)
        {
            printf("%.0f%% of iterations have been done.\n", (i / (double)maxIter) * 100);
        }
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
            printf("The probability matrix handed is:\n");
            printMatrix(currentP);
            q = computeQMultinomial(currentP);
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

        Matrix newProbability = getP(q, true);
        printf("\nThe new probability calculated with getP is:\n\n");
        printMatrix(&newProbability);
        printf("\nThe current probability is:");
        printMatrix(currentP);

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
            freeMatrix(X);
            freeMatrix(W);
            free(q);
            return newProbability;
        }

        freeMatrix(currentP);
        *currentP = createMatrix(newProbability.rows, newProbability.cols);
        memcpy(currentP->data, newProbability.data, sizeof(double) * newProbability.rows * newProbability.cols);
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
void testProb()
{

    // printf("Running test for the initial probability matrix\n The `X` matrix with the candidates votes is:\n");
    // printMatrix(X);
    // printf("\nThe `w` matrix with the groups votes is:\n");
    // printMatrix(W);

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
    if (X != NULL)
    {
        freeMatrix(X);
        X = NULL;
    }
    if (W != NULL)
    {
        freeMatrix(W);
        W = NULL;
    }
}

int main()
{
    printf("The program is running\n");

    Matrix XX = {.data = NULL, .rows = 0, .cols = 0};
    Matrix G = {.data = NULL, .rows = 0, .cols = 0};
    char *method = "multinomial";
    createInstance(&XX, &G, 42, *method); // TODO: Arreglar esto para poder crear una instancia...

    Matrix matrices[2] = {XX, G};

    writeMatrices("matricesTest.bin", matrices, 2);

    setParameters(&XX, &G);
    testProb();
    Matrix P = getInitialP("uniform");

    Matrix Pnew = EMAlgoritm(&P, "Multinomial", 0.0001, 10000, true);
    printMatrix(&Pnew);
    freeMatrix(&Pnew);
    freeMatrix(&XX);
    freeMatrix(&G);

    return 1;
}
