#include "multivariateUtils.h"
#include "globals.h"
#include "matrixUtils.h"
#include <cblas.h>
#include <lapacke.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Computes the parameters of the unconditional probability
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeRows()" function from matrixUtils.
 * @param[in, out] *mu An array of size c-1 to store the results of the average.
 * @param[in, out] *sigma A matrix of size (c-1, c-1) to store the sigma matrix.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getParams(int b, int g, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma)
{

    // mu MUST be of size TOTAL_CANDIDATES-1
    // sigma MUST be of size TOTAL_CANDIDATES-1xTOTAL_CANDIDATES-1
    if (probabilitiesReduced->cols != TOTAL_CANDIDATES - 1)
    {
        fprintf(
            stderr,
            "The probability matrix handed should consider C-1 candidates, but it has %d columns. Consider using the "
            "`removeLastCols()` function from matrixUtils.\n",
            probabilitiesReduced->cols);
        exit(EXIT_FAILURE);
    }
    if (W == NULL && X == NULL)
    {
        fprintf(stderr, "The `w` and `x` matrices aren't defined.\n");
        exit(EXIT_FAILURE);
    }

    double *groupVotesPerBallot = (double *)malloc(TOTAL_GROUPS * sizeof(double));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        groupVotesPerBallot[g] = MATRIX_AT_PTR(W, b, g);
    }

    // ---- Computation of mu ---- //
    cblas_dgemv(CblasRowMajor,              // Row-major order
                CblasNoTrans,               // Transpose the matrix (p^T)
                TOTAL_GROUPS,               // Rows of original matrix (G)
                TOTAL_CANDIDATES - 1,       // Columns of original matrix (C)
                1.0,                        // alpha
                probabilitiesReduced->data, // Matrix p
                TOTAL_GROUPS,               // Leading dimension (rows in original p)
                groupVotesPerBallot,        // Vector w_b
                1,                          // Increment for x
                0.0,                        // beta
                mu,                         // Output vector μ
                1                           // Increment for y
    );
    // free(groupVotesPerBallot);

    // ---- Computation of sigma ---- // Here there is an error
    Matrix diagonalVotesPerBallot = createDiagonalMatrix(groupVotesPerBallot, TOTAL_GROUPS);
    // printf("The diagonal votes per ballot matrix is\n");
    // printMatrix(&diagonalVotesPerBallot);
    // printf("The probabilities reduced matrix is\n");
    // printMatrix(probabilitiesReduced);

    Matrix temp = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_GROUPS);
    printf("Created temp, an empty matrix, that is:");
    printMatrix(&temp);
    printf("Going to multiplicate p^T*diag(w_b) as\n");
    printMatrix(probabilitiesReduced);
    printMatrix(&diagonalVotesPerBallot);
    // Compute p^T * diag(w_b)
    cblas_dgemm(CblasRowMajor,               // Row-major order
                CblasNoTrans,                // No transpose (Matrix A => p^T)
                CblasNoTrans,                // No transpose
                TOTAL_CANDIDATES - 1,        // Rows
                TOTAL_GROUPS,                // Columns
                TOTAL_GROUPS,                // Inner dimension
                1.0,                         // alpha
                probabilitiesReduced->data,  // Matrix p
                TOTAL_GROUPS,                // Leading dimension
                diagonalVotesPerBallot.data, // Matrix diag(w_b)
                TOTAL_GROUPS,                // Leading dimension
                0.0,                         // beta
                temp.data,                   // Output matrix
                TOTAL_GROUPS);               // Leading dimension

    // Result is (c-1, g)
    printf("The first product is\n");
    printMatrix(&temp);
    //  Compute (p^T * diag(w_b)) * p
    printf("Going to multiplicate [p^T*diag(w_b)]*p as\n");
    printMatrix(&temp);
    printMatrix(probabilitiesReduced);
    cblas_dgemm(CblasRowMajor,              // Row-major order
                CblasNoTrans,               // No transpose
                CblasTrans,                 // No transpose
                TOTAL_CANDIDATES - 1,       // Rows
                TOTAL_CANDIDATES - 1,       // Columns
                TOTAL_GROUPS,               // Inner dimension
                1.0,                        // alpha
                temp.data,                  // Matrix (p^T * diag(w_b))
                TOTAL_GROUPS,               // Leading dimension
                probabilitiesReduced->data, // Matrix p
                TOTAL_GROUPS,               // Leading dimension
                0.0,                        // beta
                sigma->data,                // Output Σ_b
                TOTAL_CANDIDATES - 1);      // Leading dimension

    printf("The second product is\n");
    printMatrix(sigma);
    //  Subtract diag(μ_b)
    printf("\nSubstracting mu of the diagonal, whom elements are:\n");
    for (int i = 0; i < TOTAL_CANDIDATES - 1; i++)
    {
        printf("%.3f, ", mu[i]);
        MATRIX_AT_PTR(sigma, i, i) = mu[i] - MATRIX_AT_PTR(sigma, i, i);
    }
    printf("\nFinal matrix is");
    printMatrix(sigma);
    freeMatrix(&temp);
    freeMatrix(&diagonalVotesPerBallot);
    free(groupVotesPerBallot);
}

/**
 * @brief Computes the parameters of the conditional probability
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the results
 * of a group.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeRows()" function from matrixUtils.
 * @param[in, out] *newMu An array of size c-1 to store the results of the average.
 * @param[in, out] *newSigma A matrix of size (c-1, c-1) to store the sigma matrix.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getAverageConditional(int b, int g, const Matrix *probabilitiesReduced, Matrix *conditionalMu,
                           Matrix *conditionalSigma)
{
    // ---- Get the parameters of the unconditional probability ---- //
    double *newMu = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
    Matrix newSigma = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
    getParams(b, g, probabilitiesReduced, newMu, &newSigma);
    // ---- ... ----

    // ---- Computation for mu ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
        {
            MATRIX_AT_PTR(conditionalMu, g, c) = newMu[c] - MATRIX_AT_PTR(probabilitiesReduced, g, c);
        }
    }
    free(newMu);
    printf("\nThe new mu for the conditional probability is mu - p_g:\n");
    printMatrix(conditionalMu);
    // ---- ... ---- //

    // ---- Get the parameters for the conditional sigma ---- //

    // 1. Get the diagonal probabilities.
    double **probabilitiesForG = (double **)malloc(TOTAL_GROUPS * sizeof(double *));
    Matrix *diagonalProbabilities = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        probabilitiesForG[g] = (double *)malloc(((TOTAL_CANDIDATES - 1) * sizeof(double)));
        for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
        {
            probabilitiesForG[g][c] = MATRIX_AT_PTR(probabilitiesReduced, g, c);
        }
        diagonalProbabilities[g] = createDiagonalMatrix(probabilitiesForG[g], TOTAL_CANDIDATES - 1);
    }

    // 2. Get the matrix multiplications, they're esentially outer products
    Matrix *matrixMultiplications = (Matrix *)malloc((TOTAL_GROUPS) * sizeof(Matrix));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        Matrix mult = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
        cblas_dger(CblasRowMajor, TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1, 1.0, probabilitiesForG[g], 1,
                   probabilitiesForG[g], 1, mult.data, TOTAL_CANDIDATES - 1);
        matrixMultiplications[g] = mult;
    }

    // 3. Add the results to the final array of matrices.
    // Esentially computes:
    // $$\sigma_b = diag(p_{g}^{t})-p^{t}_{g}p_{g}$$
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        Matrix sigmaG = copyMatrix(&newSigma);
        for (uint16_t i = 0; i < TOTAL_CANDIDATES - 1; i++)
        {
            for (uint16_t j = 0; j < TOTAL_CANDIDATES - 1; j++)
            {
                MATRIX_AT(conditionalSigma[g], i, j) =
                    MATRIX_AT(sigmaG, i, j) - MATRIX_AT(matrixMultiplications[g], i, j);
                if (i == j)
                {
                    MATRIX_AT(conditionalSigma[g], i, j) = MATRIX_AT(diagonalProbabilities[g], i, j);
                }
            }
            freeMatrix(&sigmaG);
        }
        freeMatrix(&matrixMultiplications[g]);
        freeMatrix(&diagonalProbabilities[g]);
    }
    free(matrixMultiplications);
    free(probabilitiesForG);
    freeMatrix(&newSigma);
    // ---- ... ----
}
