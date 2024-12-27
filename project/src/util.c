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

// Macro for making an easier indexation
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])

Matrix getInitialP(Matrix x, Matrix w, const char *p_method)
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
        int total = 0;

        // Step 1: Calculate the total amount of votes per candidates.
        rowSum(&x, canVotes);

        // Step 2: Calculate the total amount of votes
        for (int c = 0; c < candidates; c++)
        {
            total += canVotes[c];
        }

        // Step 3: Fill the matrix
        for (int c = 0; c < candidates; c++)
        {
            for (int g = 0; g < groups; g++)
            {
                MATRIX_AT(probabilities, g, c) = (double)canVotes[c] / (double)total;
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

Matrix getP(Matrix a, double q)
{
    return a;
}

void EMAlgoritm(Matrix x, Matrix w, Matrix initialP, const char *q_method, double convergence, int maxIter,
                bool verbose)
{

    /**
     * @brief Implements the whole EM algorithm.
     *
     * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
     * paper, it currently supports Hit and Run, Multinomial, MVN CDF and MVN PDF methods.
     *
     * @param[in] x Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
     * @param[in] w Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
     * @param[in] initialP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "Multinomial", "MVN CDF" and "MVN PDF" methods.
     * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
     * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
     * @param[in] verbose Wether to verbose useful outputs.
     *
     * @return TODO: On Python it's a dict, maybe a file with results written.
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

    if (verbose == true)
    {
        printf("Starting the EM algorithm.\n The candidates matrix is:\n");
        printMatrix(&x);
        printf("\nThe matrix with the demographic groups votation is:\n");
        printMatrix(&w);
        printf("\nThe method to calculate the conditional probability will be %s method with the following "
               "parameters:\nConvergence threshold:\t%.6f\nMaximum iterations:\t%d\n",
               q_method, convergence, maxIter);
    }

    Matrix *probabilityPtr = &initialP;
    double q = 2.0;

    for (int i; i < maxIter; i++)
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

        newProbability = getP(probabilityPtr, q);
    }
}
// SEXP hello_gsl()
//{
//     printf("Hello, World from C!\n");

//   return R_NilValue;
//
void testProb(Matrix X, Matrix G)
{

    printf("Running test for the initial probability matrix\n The `X` matrix with the candidates votes is:\n");
    printMatrix(&X);
    printf("\nThe `w` matrix with the groups votes is:\n");
    printMatrix(&G);

    Matrix prob = getInitialP(X, G, "uniform");
    printf("\nThe probability matrix for `uniform` method is:\n");
    printMatrix(&prob);
    freeMatrix(&prob);

    Matrix prob2 = getInitialP(X, G, "proportional");
    printf("\nThe probability matrix for `proportional` method is:\n");
    printMatrix(&prob2);
    freeMatrix(&prob2);

    Matrix prob3 = getInitialP(X, G, "group proportional");
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

    memcpy(X.data, xVal, sizeof(xVal));
    memcpy(G.data, gVal, sizeof(gVal));

    testProb(X, G);
    freeMatrix(&X);
    freeMatrix(&G);
    end = clock();
    int t = (end - start) / CLOCKS_PER_SEC;
    printf("The program took %d seconds!", t);
    return 0;
}
