#include "multivariateUtils.h"
#include <cblas.h>
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
    if (probabilitiesReduced->rows != TOTAL_CANDIDATES - 1)
    {
        fprintf(stderr, "The probability matrix handed should consider C-1 candidates. Consider using the "
                        "`removeRows()` function from matrixUtils.\n");
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
    cblas_dgemv(CblasRowMajor,        // Row-major order
                CblasTrans,           // Transpose the matrix (p^T)
                TOTAL_GROUPS,         // Rows of original matrix (G)
                TOTAL_CANDIDATES - 1, // Columns of original matrix (C)
                1.0,                  // alpha
                probabilitiesReduced, // Matrix p
                TOTAL_GROUPS,         // Leading dimension (rows in original p)
                groupVotesPerBallot,  // Vector w_b
                1,                    // Increment for x
                0.0,                  // beta
                mu,                   // Output vector μ
                1                     // Increment for y
    );
    // free(groupVotesPerBallot);

    // ---- Computation of sigma ---- //
    Matrix diagonalVotesPerBallot = createDiagonalMatrix(groupVotesPerBallot, TOTAL_GROUPS);

    Matrix temp = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_GROUPS);

    // Compute p^T * diag(w_b)
    cblas_dgemm(CblasRowMajor,               // Row-major order
                CblasNoTrans,                // No transpose
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

    // Compute (p^T * diag(w_b)) * p
    cblas_dgemm(CblasRowMajor,              // Row-major order
                CblasNoTrans,               // No transpose
                CblasNoTrans,               // No transpose
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

    // Subtract diag(μ_b)
    for (int i = 0; i < TOTAL_CANDIDATES - 1; i++)
    {
        MATRIX_AT_PTR(sigma, i, i) = mu[i] - MATRIX_AT_PTR(sigma, i, i);
    }

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

void getAverageConditional(int b, int g, const Matrix *probabilitiesReduced, double *newMu, Matrix *newSigma)
{
    getParams(b, g, probabilitiesReduced, newMu, newSigma);

    // Extract p_g^* (Group vector from probabilitiesReduced)

    // ---- Computation for mu ----
    double *probabilitiesForAGroup = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        probabilitiesForAGroup[c] = MATRIX_AT_PTR(probabilitiesReduced, c, g);
    }

    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        newMu[c] -= probabilitiesForAGroup[c];
    }

    // ---- Computation for sigma ----
    // Subtract diag(p_g^T)
    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        MATRIX_AT_PTR(newSigma, c, c) -= probabilitiesForAGroup[c];
    }

    // Subtract p_g^T * p_g (Outer Product)
    Matrix outerProduct = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);

    cblas_dger(CblasRowMajor,          // Row-major order
               TOTAL_CANDIDATES - 1,   // Number of rows
               TOTAL_CANDIDATES - 1,   // Number of columns
               1.0,                    // Alpha (scalar multiplier)
               probabilitiesForAGroup, // Vector x (row vector)
               1,                      // Increment for x
               probabilitiesForAGroup, // Vector y (column vector)
               1,                      // Increment for y
               outerProduct.data,      // Output matrix (result of p_g^T * p_g)
               TOTAL_CANDIDATES - 1    // Leading dimension
    );

#pragma omp for collapse(2)
    for (uint16_t i = 0; i < TOTAL_CANDIDATES - 1; i++)
    {
        for (uint16_t j = 0; j < TOTAL_CANDIDATES - 1; j++)
        {
            MATRIX_AT_PTR(newSigma, i, j) -= MATRIX_AT(outerProduct, i, j);
        }
    }

    free(probabilitiesForAGroup);
    freeMatrix(&outerProduct);
}
