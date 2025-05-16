/*
Copyright (c) 2025 fastei team

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

#include "utils_multivariate.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <Rinternals.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

#ifndef BLAS_INT
#define BLAS_INT int
#endif
/**
 * @brief Computes the parameters of the unconditional probability without the last candidate
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

void getParams(int b, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma)
{

    // ---- Check parameters ---- //
    // ---- Note: Mu must be of size TOTAL_CANDIDATES-1 and sigma of size (TOTAL_CANDIDATES-1xTOTAL_CANDIDATES-1) ----
    if (probabilitiesReduced->cols != TOTAL_CANDIDATES - 1)
    {
        error("Multivariate utils: The probability matrix handed should consider C-1 candidates, but it has %d "
              "columns. Consider using the "
              "`removeLastCols()` function from matrixUtils.\n",
              probabilitiesReduced->cols);
    }
    if (W == NULL && X == NULL)
    {
        error("Multivariate utils: The `w` and `x` matrices aren't defined.\n");
    }
    // --- ... --- //

    // --- Calculations --- //
    // ---- The votes that a group has made on a given ballot ----
    double *groupVotesPerBallot = getRow(W, b);

    // ---- Computation of mu ----
    // ---- Performing the matrix multiplication of p^T * w_b

    char trans = 'T';
    int m = TOTAL_GROUPS;         // G
    int n = TOTAL_CANDIDATES - 1; // C - 1
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1, incy = 1;

    // Remember: Fortran BLAS is always column-major.
    int lda = m; // = G

    F77_CALL(dgemv)
    (&trans, // 'T'
     &m,     // G
     &n,     // C-1
     &alpha, probabilitiesReduced->data,
     &lda, // G
     groupVotesPerBallot, &incx, &beta, mu, &incy FCONE);

    // ---- Computation of sigma ----
    // ---- Get a diagonal matrix with the group votes on a given ballot ----
    Matrix diagonalVotesPerBallot = createDiagonalMatrix(groupVotesPerBallot, TOTAL_GROUPS);
    // ---- Temporary matrix to store results ----
    Matrix temp = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_GROUPS);
    // ---- Calculates the matrix multiplication of p^T * diag(w_b); result must be (C-1 x G) ----

    char transA = 'T';        // p^T
    char transB = 'N';        // diag(w_b) as-is
    m = TOTAL_CANDIDATES - 1; // (C-1)
    n = TOTAL_GROUPS;         // G
    int k = TOTAL_GROUPS;     // G

    // Leading dimensions in column-major:
    lda = TOTAL_GROUPS;     // p is G x (C-1) => LDA = G
    int ldb = TOTAL_GROUPS; // diag(w_b) is G x G => LDB = G
    int ldc = m;            // = C-1, for the result which is (C-1) x G in column-major

    F77_CALL(dgemm)
    (&transA, &transB, &m, &n, &k, &alpha, probabilitiesReduced->data, &lda, diagonalVotesPerBallot.data, &ldb, &beta,
     temp.data, &ldc FCONE FCONE);

    transA = 'N';             // no transpose
    transB = 'N';             // no transpose
    m = TOTAL_CANDIDATES - 1; // C-1
    n = TOTAL_CANDIDATES - 1; // C-1
    k = TOTAL_GROUPS;         // G
    alpha = 1.0;
    beta = 0.0;

    lda = m;            // = C-1, A is (C-1) x G
    ldb = TOTAL_GROUPS; // B is G x (C-1)
    ldc = m;            // = C-1, C is (C-1) x (C-1)

    F77_CALL(dgemm)
    (&transA, &transB, &m, &n, &k, &alpha, temp.data, &lda, probabilitiesReduced->data, &ldb, &beta, sigma->data,
     &ldc FCONE FCONE);

    // ---- Substract the diagonal with the average ----
    // ---- Note: This could be optimized with a cBLAS call too ----

    for (int j = 0; j < TOTAL_CANDIDATES - 1; j++)
    { // ---- For each candidate
        for (int i = 0; i < TOTAL_CANDIDATES - 1; i++)
        { // ---- For each candidate given another candidate
            if (i == j)
            { // ---- If it corresponds to a diagonal, substract diagonal
                MATRIX_AT_PTR(sigma, i, j) = mu[i] - MATRIX_AT_PTR(sigma, i, j);
                continue;
            }
            MATRIX_AT_PTR(sigma, i, j) = -MATRIX_AT_PTR(sigma, i, j);
        }
    }
    //  ---- Free alocated memory ----
    Free(groupVotesPerBallot);
    freeMatrix(&temp);
    freeMatrix(&diagonalVotesPerBallot);
    // --- ... --- //
}

/**
 * @brief Computes the parameters of the conditional probability WITHOUT the last candidate
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

void getAverageConditional(int b, const Matrix *probabilitiesReduced, Matrix *conditionalMu, Matrix **conditionalSigma)
{
    // ---- Get the parameters of the unconditional probability ---- //
    double *newMu = (double *)Calloc((TOTAL_CANDIDATES - 1), double);
    Matrix newSigma = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
    getParams(b, probabilitiesReduced, newMu, &newSigma);
    // ---- ... ----

    // ---- Computation for mu ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
        { // ---- For each candidate given a group
            MATRIX_AT_PTR(conditionalMu, g, c) = newMu[c] - MATRIX_AT_PTR(probabilitiesReduced, g, c);
        }
    }
    // ---- The original mu isn't needed anymore ---- //
    Free(newMu);
    //  ---- ... ---- //

    // ---- Get the parameters for the conditional sigma ---- //

    // ---- Get the diagonal probabilities ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store the probabilities for a given group ----
    double **probabilitiesForG = (double **)Calloc(TOTAL_GROUPS, double *);
    // ---- Create an array of size `TOTAL_GROUPS` that will store diagonal matrices with the probabilities ----
    Matrix *diagonalProbabilities = (Matrix *)Calloc((TOTAL_GROUPS), Matrix);

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        probabilitiesForG[g] = getRow(probabilitiesReduced, g);
        diagonalProbabilities[g] = createDiagonalMatrix(probabilitiesForG[g], TOTAL_CANDIDATES - 1);
    }
    // --- ... --- //

    // ---- Get the matrix multiplications ---- //
    // ---- This multiplications are esentially outer products ----
    // ---- Create an array of size `TOTAL_GROUPS` that will store each outer product ----
    Matrix *matrixMultiplications = (Matrix *)Calloc((TOTAL_GROUPS), Matrix);

    char trans = 'N';             // (There's no separate transpose flag in dger,
                                  //  but we keep a placeholder for clarity)
    int m = TOTAL_CANDIDATES - 1; // M
    int n = TOTAL_CANDIDATES - 1; // N
    double alpha = 1.0;
    int incx = 1, incy = 1;

    // For column-major, LDA = #rows = m
    int lda = m;

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        matrixMultiplications[g] = createMatrix(m, n);

        F77_CALL(dger)
        (&m, &n, &alpha, probabilitiesForG[g], &incx, probabilitiesForG[g], &incy, matrixMultiplications[g].data, &lda);

        Free(probabilitiesForG[g]);
    }

    // --- ... --- //

    // ---- Add the results to the final array of matrices ----
    // ---- Esentially computes: $$\sigma_b = diag(p_{g}^{t})-p^{t}_{g}p_{g}$$ ----
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        for (uint16_t i = 0; i < TOTAL_CANDIDATES - 1; i++)
        { // ---- For each candidate given a group
            for (uint16_t j = 0; j < TOTAL_CANDIDATES - 1; j++)
            { // ---- For each candidate given a group and a candidate
                // ---- Add the multiplication of probabilities ----
                MATRIX_AT_PTR(conditionalSigma[g], i, j) =
                    MATRIX_AT(newSigma, i, j) + MATRIX_AT(matrixMultiplications[g], i, j);
                if (i == j)
                { // ---- If it's a diagonal
                    // ---- Substract the diagonal probabilities ----
                    MATRIX_AT_PTR(conditionalSigma[g], i, j) -= MATRIX_AT(diagonalProbabilities[g], i, j);
                }
            }
        }
        // ---- Free unnecesary space ----
        freeMatrix(&matrixMultiplications[g]);
        freeMatrix(&diagonalProbabilities[g]);
    }
    // ---- Free space ----
    Free(matrixMultiplications);
    Free(diagonalProbabilities);
    Free(probabilitiesForG);
    freeMatrix(&newSigma);

    // --- ... --- //
}

/**
 * @brief Computes the Mahalanobis distance with last candidate adjustment.
 *
 * @param[in] x Pointer to the input feature vector (size C-1).
 * @param[in] mu Pointer to the mean vector (size C-1).
 * @param[in] inverseSigma Pointer to the inverse covariance matrix (size (C-1) x (C-1)).
 * @param[out] maha Pointer to the resulting Mahalanobis distance array (could be C-1 or C)
 * @param[in] size Size of the truncated candidate space (C-1).
 * @param[in] truncate Boolean value to see if the *maha parameter fixes for a lineally dependent sistem and returns a
 * size `C` array.
 */
void getMahanalobisDist(double *x, double *mu, Matrix *inverseSigma, double *maha, int size, bool reduced)
{
    // ---- Initialize temporary arrays ----
    double diff[size];
    double temp[size];

    // ---- Compute the difference vector ---- //
    for (int i = 0; i < size; i++)
    { // ---- For each truncated element
        diff[i] = x[i] - mu[i];
    }
    // --- ... --- //

    // ---- Compute the multiplication  ---- //
    // ---- inverseSigma * diff

    // ---- Note: The upper triangle is filled on the inverse sigma aswell for the Cholensky method.

    char uplo = 'L'; // lower triangle
    int n = size;
    double alpha = 1.0;
    double beta = 0.0;
    int incx = 1, incy = 1;

    // In column-major, LDA = #rows = n
    int lda = n;

    F77_CALL(dsymv)(&uplo, &n, &alpha, inverseSigma->data, &lda, diff, &incx, &beta, temp, &incy FCONE);

    // --- ... --- //

    // ---- Compute Mahalanobis distance (truncated) ---- //
    double mahanobisTruncated = 0.0;
    for (int i = 0; i < size; i++)
    { // ---- For each truncated element
        // ---- The first parenthesis ----
        maha[i] = diff[i] * temp[i]; // Store intermediate results
        mahanobisTruncated += maha[i];
    }
    // --- ... ---//
    if (!reduced)
    {
        // ---- Compute the Mahanalobis distance with the last candidate ---- //
        /*
         * The updated mahanalobis distance for "C" candidates can be written as:
         *
         * $$D_{i}^{2}=D^{2}_{baseline}-\sigma^{-1}(x_i-\mu_i)+diag(\sigma^{-1})$$
         *
         * The baseline would be the mahanobis distance for the "C-1" candidate (mahanobisTruncated)
         * */
        maha[size] = mahanobisTruncated; // Last element is used as a reference
        for (int c = 0; c < size; c++)
        { // ---- For each candidate (doesn't consider the last one)
            maha[c] = mahanobisTruncated - 2 * temp[c] + MATRIX_AT_PTR(inverseSigma, c, c);
        }
    }
    // ---...--- //
}
