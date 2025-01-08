#include "globals.h"
#include "matrixUtils.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>

void *getParams(int b, int g, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma)
{
    // mu MUST be of size TOTAL_CANDIDATES-1
    // sigma MUST be of size TOTAL_CANDIDATES-1
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

double *getAverageConditional(int b, int g, const Matrix *probabilitiesReduced)
{
    double *unconditionalMu = getAverageConditional(b, g, probabilitiesReduced);

    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        unconditionalMu[c] -= MATRIX_AT_PTR(probabilitiesReduced, g, c);
    }
    return unconditionalMu;
}

double *getSigma(int b, int g, const double *mu, const Matrix)
