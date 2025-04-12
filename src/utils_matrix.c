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

#include "utils_matrix.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <math.h>
#ifdef _OPENMP
#include <omp.h> // Parallelization
#endif
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef BLAS_INT
#define BLAS_INT int
#endif
/**
 * @brief Transpose from row-major -> column-major (or vice versa) into a separate buffer.
 *
 * @param[in]  src   Pointer to the source matrix data
 * @param[in]  rows  Number of rows in the source
 * @param[in]  cols  Number of cols in the source
 * @param[out] dst   Pointer to the destination buffer (size rows*cols)
 *
 * After this, dst will hold the transpose of src.
 */
static void transposeMatrix(const double *src, int rows, int cols, double *dst)
{
    for (int r = 0; r < rows; r++)
    {
        for (int c = 0; c < cols; c++)
        {
            dst[c * rows + r] = src[r * cols + c];
        }
    }
}

// ----------------------------------------------------------------------------
// Utility functions
// ----------------------------------------------------------------------------

/**
 * @brief Make an array of a constant value.
 *
 * Given a value, it fills a whole array with a constant value.
 *
 * @param[in, out] array Pointer to the array to be filled.
 * @param[in] N The size of the array.
 * @param[in] value The constant value to fill
 *
 * @return void Written on the input array
 *
 * @note
 * - It uses cBLAS for optimization
 *
 * @example
 * Example usage:
 * @code
 * double array[5];
 * makeArray(array, 5, 3.14);
 * // array -> now constains [3.14, 3.14, ..., 3.14]
 * @endcode
 */

void makeArray(double *array, int N, double value)
{
    if (!array)
    {
        error("Matrix handling: A NULL pointer was handed as an array.\n");
    }

    if (N < 0)
    {
        error("Matrix handling: A incoherent dimension was handen for making the array.\n");
    }

    // Fill the array with the specified constant value
    for (int i = 0; i < N; i++)
    {
        array[i] = value;
    }
}

/**
 * @brief Checks if the matrix is well defined
 *
 * Given a pointer to a matrix, it verifies if the matrix is well alocated and defined and throws an error if there's
 * something wrong.
 *
 * @param[in] m A pointer to the matrix
 *
 * @return void
 *
 * @note
 * - This will just throw errors, note that EXIT_FAILURE will dealocate memory
 *
 * @warning
 * - The pointer may be NULL.
 * - The dimensions may be negative.
 */

void checkMatrix(const Matrix *m)
{

    // Validation, checks NULL pointer
    if (!m || !m->data)
    {
        error("Matrix handling: A NULL pointer was handed as a matrix argument.\n");
    }

    // Checks dimensions
    if (m->rows <= 0 || m->cols <= 0)
    {
        error("Matrix handling: Invalid matrix dimensions: rows=%d, cols=%d\n", m->rows, m->cols);
    }
}

/**
 * @brief Creates an empty dynamically allocated memory matrix of given dimensions.
 *
 * Given certain dimensions of rows and colums, creates an empty Matrix with allocated memory towards the data.
 *
 * @param[in] rows The number of rows of the new matrix.
 * @param[in] cols The number of columns of the new matrix.
 *
 * @return Matrix Empty matrix of dimensions (rows x cols) with allocated memory for its data.
 *
 * @note
 * - Remember to free the memory! It can be made with freeMatrix() call
 *
 * @warning
 * - The memory may be full.
 * - If dimensions are negative.
 */

Matrix createMatrix(int rows, int cols)
{
    if (rows <= 0 || cols <= 0)
    {
        error("Matrix handling: Invalid matrix dimensions: rows=%d, cols=%d\n", rows, cols);
    }

    Matrix m;
    m.rows = rows;
    m.cols = cols;

    m.data = Calloc(rows * cols, double);

    if (!m.data)
    {
        error("Matrix handling: Failed to allocate matrix data\n");
    }

    return m;
}

/**
 * @brief Liberates the allocated matrix data.
 *
 * @param[in] m The matrix to free the data.
 *
 * @return void Changes to be made on the input matrix and memory.
 *
 */

void freeMatrix(Matrix *m)
{
    // TODO: Implement a validation warning.
    if (m != NULL && m->data != NULL)
    {
        Free(m->data);
        m->data = NULL;
    }
    m->rows = 0;
    m->cols = 0;
}

/**
 * @brief Prints the matrix data.
 *
 * @param[in] m The matrix to print the data.
 *
 * @return void No return, prints a message on the console.
 *
 * @note
 * - Use the function mainly for debugging.
 */

void printMatrix(const Matrix *matrix)
{
    checkMatrix(matrix); // Assertion

    Rprintf("Matrix (%dx%d):\n", matrix->rows, matrix->cols);

    for (int i = 0; i < matrix->rows; i++)
    {
        Rprintf("| ");
        for (int j = 0; j < matrix->cols - 1; j++)
        {
            Rprintf("%.3f\t", MATRIX_AT_PTR(matrix, i, j));
        }
        Rprintf("%.3f", MATRIX_AT_PTR(matrix, i, matrix->cols - 1));
        Rprintf(" |\n");
    }
}

/**
 * @brief Computes a row-wise sum.
 *
 * Given a matrix, it computes the sum over all the rows and stores them in an array.
 * @param[in] matrix Pointer to the input matrix.
 * @param[out] result Pointer of the resulting array of length `rows`.
 *
 * @return void Written on *result
 *
 * @note
 * - Matrix should be in col-major order
 * - This function uses cBLAS library, where the operation can be written as a matrix product
 *   of X * 1.
 * - Just support double type
 *
 * @example
 * Example usage:
 * @code
 *
 * double data[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 *
 * Matrix matrix = {
 * .data = values,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * double result[matrix->rows]
 *
 * rowSum(matrix, result);
 * // result now contains [6.0, 15.0]
 * @endcode
 */

void rowSum(const Matrix *matrix, double *result)
{
    checkMatrix(matrix); // Assertion

    // We will use malloc to avoid stack overflows. Usually, we don't know
    // what's the size that will be given to the matrix, but if it exceeds
    // 1.000.000 columns (8MB for a `double` type matrix) it will overflow

    double *ones = (double *)Calloc(matrix->cols, double);
    if (!ones)
    {
        // In theory, R's Calloc would take care of that and throw its own error
        error("Matrix handling: Failed to allocate memory to the rowSum function. \n");
    }

    makeArray(ones, matrix->cols, 1.0);
    // BLAS parameters
    char trans = 'T'; // No transpose needed (matrix already column-major)
    double alpha = 1.0, beta = 0.0;
    int incX = 1, incY = 1;

    // Perform Matrix-Vector Multiplication (Matrix * Ones = Row Sums)
    F77_CALL(dgemv)
    (&trans, &(matrix->rows), &(matrix->cols), // Matrix dimensions (M, N)
     &alpha, matrix->data, &(matrix->rows),    // Matrix and leading dimension
     ones, &incX,                              // Vector of ones
     &beta, result, &incY FCONE);              // Output row sum

    Free(ones);
}

/**
 * @brief Computes a column-wise sum.
 *
 * Given a matrix, it computes the sum over all the columns and stores them in an array.
 *
 * @param[in] matrix Pointer to a Matrix structure.
 * @param[out] result Array for the resulting Matrix. It must have the dimensions of the matrix columns.
 *
 * @return void Written on *result
 *
 * @note
 * - Matrix should be in col-major order
 * - It will use cBLAS for operations, where it will do a matrix product of X * 1^T. It'll use the
 *   already algorithm implemented in rowSum
 * - Just support double as a result due to cBLAS
 *
 * @warning
 * - The matrix or array pointer may be NULL.
 *
 * @example
 * Example usage:
 * @code
 *
 * double data[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 *
 * Matrix matrix = {
 * .data = data,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * double result[matrix->cols]
 *
 * colSum(matrix, result);
 * // result now contains [5.0, 7.0, 9.0]
 * @endcode
 */

void colSum(const Matrix *matrix, double *result)
{
    checkMatrix(matrix); // Assertion

    double *ones = (double *)Calloc(matrix->rows, double);
    if (!ones)
    {
        error("Matrix handling: Failed to allocate memory in colSum function.\n");
    }

    makeArray(ones, matrix->rows, 1.0);

    // BLAS parameters
    char trans = 'N'; // Transpose to sum columns
    double alpha = 1.0, beta = 0.0;
    int incX = 1, incY = 1;

    // Perform Matrix-Vector Multiplication: (Transpose(Matrix) * Ones = Column Sums)
    F77_CALL(dgemv)
    (&trans, &(matrix->rows), &(matrix->cols), // (M, N) dimensions
     &alpha, matrix->data, &(matrix->rows),    // Matrix and leading dimension
     ones, &incX,                              // Vector of ones
     &beta, result, &incY FCONE);              // Output column sum

    Free(ones);
}

/**
 * @brief Fills matrix with a constant value.
 *
 * Given a matrix, it fills a whole matrix with a constant value.
 *
 * @param[in, out] matrix Pointer to matrix to be filled.
 * @param[in] value The constant value to fill
 *
 * @return void Written on the input matrix
 *
 * @note
 * - Matrix should be in col-major order.
 *
 * @example
 * Example usage:
 * @code
 * double values[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 * Matrix matrix = {
 * .data = values,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * fillMatrix(matrix, 9);
 * // matrix->data now contains [9.0, 9.0, 9.0, ..., 9.0]
 * @endcode
 */

void fillMatrix(Matrix *matrix, const double value)
{
    checkMatrix(matrix); // Assertion
    int size = matrix->rows * matrix->cols;

    makeArray(matrix->data, size, value);
}

/**
 * @brief Checks if the difference of two matrices converge to a value
 *
 * Given two matrices, it performs de absolute difference and evaluate the convergence towards a given
 * arbitrary values: |x1 - x2| < epsilon. If there's a value whom convergence is greater than epsilon, the convergence
 * is not achieved.
 *
 * @param[in] matrix Matrix to perform the substraction.
 * @param[in] matrix Matrix to perform the substraction.
 * @param[in] double Arbitrary value to evaluate the convergence
 *
 * @return bool Boolean value to see if it converges.
 *
 * @warning:
 * - Both matrices should be from the same dimention.
 * @note
 * - Matrix should be in col-major order.
 *
 * @example
 * Example usage:
 * @code
 * double values[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 * double values2[6] = {
 * 		1.1, 2.1, 2.9,
 * 		3.9, 5.1, 6.1
 * }
 * Matrix matrix = {
 * .data = values,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * Matrix matrix2 = {
 * .data = values2,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * bool converges = convergeMatrix(matrix, matrix2,  0.02);
 * // bool->true
 * @endcode
 */
bool convergeMatrix(const Matrix *matrixA, const Matrix *matrixB, const double convergence)
{

    checkMatrix(matrixA);
    checkMatrix(matrixB);

    if (convergence <= 0)
    {
        return false; // We will acept -Inf values
    }

    if (matrixA->cols != matrixB->cols || matrixA->rows != matrixB->rows)
    {
        error("Matrix handling: The dimensions of both matrices doesn't match.\n");
    }

    int size = matrixA->rows * matrixB->cols;
    int incX = 1;
    int incY = 1;
    double alpha = -1.0;

    double *diff = (double *)Calloc(size, double);

    F77_CALL(dcopy)(&(size), matrixA->data, &incX, diff, &incY);
    F77_CALL(daxpy)(&(size), &alpha, matrixB->data, &incX, diff, &incY);
    for (int i = 0; i < size; i++)
    {
        // If there's a value whom convergence is greater than epsilon, the convergence
        // isn't achieved.
        if (fabs(diff[i]) >= convergence)
        {
            Free(diff);
            return false;
        }
    }

    Free(diff);
    return true;
}

/**
 * @brief Retrieves the maximum element of the matrix.
 *
 * @param[in] matrix Matrix to find the maximum element.
 *
 * @return double The maximum element
 *
 * @note
 * - Matrix should be in col-major order.
 *
 * @example
 * Example usage:
 * @code
 * double values[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };

 * Matrix matrix = {
 * .data = values,
 * .rows = 2,
 * .cols = 3
 * }
 *
 * double maximum = maxElement(&matrix);
 *
 * // maximum=6.0
 * @endcode
 */

double maxElement(const Matrix *m)
{

    checkMatrix(m);
    int size = m->cols * m->rows;

    double max = m->data[0];
    for (int i = 0; i < size; i++)
    {
        if (max < m->data[i])
        {
            max = m->data[i];
        }
    }
    return max;
}

/**
 * @brief Removes the last row of a matrix.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @return Matrix A new matrix with one less row.
 *
 * @note
 * - The original matrix remains unchanged.
 * - The memory for the new matrix is dynamically allocated; remember to free it.
 */
Matrix removeLastRow(const Matrix *matrix)
{
    checkMatrix(matrix); // Validate the matrix

    if (matrix->rows <= 1)
    {
        error("Matrix handling: Matrix must have at least two rows to remove one.\n");
    }

    // Create a new matrix with one less row
    Matrix newMatrix = createMatrix(matrix->rows - 1, matrix->cols);

    // Copy all rows except the last one
    for (int i = 0; i < matrix->rows - 1; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            MATRIX_AT(newMatrix, i, j) = MATRIX_AT_PTR(matrix, i, j);
        }
    }

    return newMatrix;
}

/**
 * @brief Removes the last column of a matrix.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @return Matrix A new matrix with one less column.
 *
 * @note
 * - The original matrix remains unchanged.
 * - The memory for the new matrix is dynamically allocated; remember to free it.
 */
Matrix removeLastColumn(const Matrix *matrix)
{
    checkMatrix(matrix); // Validate the matrix

    if (matrix->cols <= 1)
    {
        error("Matrix handling: Matrix must have at least two columns to remove one.\n");
    }

    // Create a new matrix with one less column
    Matrix newMatrix = createMatrix(matrix->rows, matrix->cols - 1);

    // Copy all columns except the last one
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols - 1; j++)
        {
            MATRIX_AT(newMatrix, i, j) = MATRIX_AT_PTR(matrix, i, j);
        }
    }

    return newMatrix;
}

/**
 * @brief Creates a diagonal matrix given a 1D array
 *
 * @param[in] vector Pointer to the array.
 * @param[in] size The size of the array.
 * @return Matrix A new matrix (size x size) with each element of the array as a diagonal.
 *
 * @note
 * - The original array remains unchanged.
 */
Matrix createDiagonalMatrix(const double *vector, int size)
{
    Matrix diag = createMatrix(size, size);

    if (!diag.data)
    {
        error("Matrix handling: Failed to allocate memory for diagonal matrix.\n");
    }

    for (int i = 0; i < size; i++)
    {
        MATRIX_AT(diag, i, i) = vector[i];
    }

    return diag;
}

/**
 * @brief Computes the inverse of a symmetric positive-definite matrix using Cholesky decomposition.
 *
 * @param[in, out] matrix Pointer to the input matrix (overwritten with the inverse).
 *
 * @note The input matrix must be square and symmetric.
 */
void inverseSymmetricPositiveMatrix(Matrix *matrix)
{
    checkMatrix(matrix);
    char lCh = 'L';
    if (matrix->rows != matrix->cols)
    {
        error("Matrix handling: Matrix must be square.\n");
    }

    int n = matrix->rows;
    if (n == 1)
    {
        double val = matrix->data[0];
        if (val != 0.0)
            matrix->data[0] = 1.0 / val;
        return;
    }

    // copy for "emergency" fallback
    Matrix emergencyMat = copyMatrix(matrix);

    int info;
    F77_CALL(dpotrf)(&lCh, &n, matrix->data, &n, &info FCONE);
    // Cholesky => ('L', n, matrix->data)

    if (info < 0)
    {
        error("Matrix handling: dpotrf illegal value argument. info=%d\n", info);
    }
    if (info > 0)
    {
        Rprintf("Cholesky decomposition failed. Leading minor not positive definite.\n"
                "Retrying with +1 on diagonal.\n");

        // restore from emergencyMat, add small diagonal
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(emergencyMat, i, j);
                if (i == j)
                    MATRIX_AT_PTR(matrix, i, j) += 1.0;
            }
        }
        freeMatrix(&emergencyMat);
        inverseSymmetricPositiveMatrix(matrix);
        return; // be sure to stop here so the next steps don't run again
    }

    // Now invert the Cholesky => rowMajor_dpotri('L', n, matrix->data)
    int info2;
    F77_CALL(dpotri)(&lCh, &n, matrix->data, &n, &info2 FCONE);
    if (info2 != 0)
    {
        error("Matrix handling: Matrix inversion failed after Cholesky. info=%d\n", info2);
    }

    // Fill upper triangle = mirror of lower triangle
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT_PTR(matrix, j, i);
        }
    }
    freeMatrix(&emergencyMat);
}

/**
 * @brief Inverts a real symmetric NxN matrix (overwrites the input).
 *
 * Uses an eigen-decomposition (dsyev) to invert A = Q * diag(vals) * Q^T.
 * The input matrix must be square and invertible (no zero eigenvalues).
 *
 * @param[in,out] matrix Pointer to the NxN symmetric matrix in col-major layout.
 */

void inverseMatrixEigen(Matrix *matrix)
{
    checkMatrix(matrix);
    if (matrix->rows != matrix->cols)
    {
        error("Matrix handling: The matrix must be square for getting the inverse and Eigen values.\n");
    }
    int n = matrix->rows;

    double *eigenvals = (double *)Calloc(n, double);
    if (!eigenvals)
    {
        error("Matrix handling: Cannot allocate eigenvals.\n");
    }

    // First, we query the working array to get the optimal size...
    int info;
    double query_work;
    int lwork = -1; // Query mode
    char vCh = 'V';
    char uCh = 'U';

    F77_CALL(dsyev)(&vCh, &uCh, &n, matrix->data, &n, eigenvals, &query_work, &lwork, &info FCONE FCONE);

    if (info != 0)
    {
        error("Matrix handling: dsyev workspace query failed with info = %d\n", info);
    }

    // Allocate the working array
    lwork = (int)query_work; // LAPACK returns optimal size in query_work
    double *work = (double *)Calloc(lwork, double);

    // Step 3: Compute Eigen decomposition
    F77_CALL(dsyev)(&vCh, &uCh, &n, matrix->data, &n, eigenvals, work, &lwork, &info FCONE FCONE);

    if (info != 0)
    {
        error("Matrix handling: dsyev failed with info = %d\n", info);
    }

    // Invert eigenvalues
    for (int i = 0; i < n; i++)
    {
        if (fabs(eigenvals[i]) < 1e-15)
        {
            error("Matrix handling: Zero or near-zero eigenvalue => not invertible.\n");
            Free(eigenvals);
        }
        eigenvals[i] = 1.0 / eigenvals[i];
    }

    // Build diagonal matrix from eigenvals
    Matrix Dinv = createDiagonalMatrix(eigenvals, n);
    Free(eigenvals);

    // temp = Q * Dinv
    Matrix temp = createMatrix(n, n);
    double alpha = 1.0;
    double beta = 1.0;
    char noTranspose = 'N';
    char yTranspose = 'T';

    F77_CALL(dgemm)
    (&noTranspose, &noTranspose,           // Column-major perspective => No-trans x No-trans
     &n, &n, &n, &alpha, matrix->data, &n, // Q (eigenvectors) in column-major
     Dinv.data, &n, &beta, temp.data, &n FCONE FCONE);

    // A_inv = temp * Q^T
    Matrix temp2 = createMatrix(n, n);
    F77_CALL(dgemm)
    (&noTranspose, &yTranspose,         // Column-major perspective => No-trans x Transposed
     &n, &n, &n, &alpha, temp.data, &n, // `temp` in column-major
     matrix->data, &n,                  // `Q^T` is achieved by using 'T' (BLAS transposes internally)
     &beta, temp2.data, &n FCONE FCONE);

    // copy result back
    memcpy(matrix->data, temp2.data, n * n * sizeof(double));

    freeMatrix(&temp2);
    freeMatrix(&temp);
    freeMatrix(&Dinv);
}

/**
 * @brief Computes the inverse of a general square matrix using LU decomposition.
 *
 * @param[in, out] matrix Pointer to the input matrix (overwritten with the inverse).
 *
 * @note The input matrix must be square and invertible.
 */
void inverseMatrixLU(Matrix *matrix)
{
    checkMatrix(matrix);
    if (matrix->rows != matrix->cols)
    {
        error("Matrix handling: Matrix must be square for inversion.\n");
    }
    int n = matrix->rows;

    int *ipiv = (int *)Calloc(n, int);
    if (!ipiv)
    {
        error("Matrix handling: Failed to allocate pivot array.\n");
    }

    // LU decomposition
    int info;
    F77_CALL(dgetrf)(&n, &n, matrix->data, &n, ipiv, &info);
    if (info != 0)
    {
        error("Matrix handling: LU decomposition failed. info=%d\n", info);
        Free(ipiv);
    }

    // Invert from LU
    // ---- Needs to assign a workspace ----
    // Query optimal workspace size
    double query_work;
    int lwork = -1;
    F77_CALL(dgetri)(&n, matrix->data, &n, ipiv, &query_work, &lwork, &info);

    // Allocate optimal work array
    lwork = (int)query_work;
    double *work = (double *)Calloc(lwork, double);

    if (!work)
    {
        error("Matrix handling: Failed to allocate workspace for dgetri.\n");
    }

    // Compute matrix inverse
    F77_CALL(dgetri)(&n, matrix->data, &n, ipiv, work, &lwork, &info);

    // Free workspace
    Free(work);
    if (info != 0)
    {
        error("Matrix handling: Matrix inversion failed. info=%d\n", info);
        Free(ipiv);
    }

    Free(ipiv);
}

Matrix copyMatrix(const Matrix *original)
{
    checkMatrix(original); // Ensure the original matrix is valid

    // Create a new matrix with the same dimensions
    Matrix copy = createMatrix(original->rows, original->cols);

    // Copy the data from the original matrix
    for (int i = 0; i < original->rows; i++)
    {
        for (int j = 0; j < original->cols; j++)
        {
            MATRIX_AT(copy, i, j) = MATRIX_AT_PTR(original, i, j);
        }
    }

    return copy;
}

/**
 * @brief Extracts the n-th row of a matrix as a dynamically allocated array.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @param[in] rowIndex The index of the row to extract (0-based).
 * @return double* A dynamically allocated array containing the row elements.
 *
 * @note The caller is responsible for freeing the returned array.
 */
double *getRow(const Matrix *matrix, int rowIndex)
{
    checkMatrix(matrix); // Ensure the matrix is valid

    if (rowIndex < 0 || rowIndex >= matrix->rows)
    {
        error("Matrix handling: Row index out of bounds: %d\n", rowIndex);
    }

    // Allocate memory for the row
    double *row = (double *)Calloc(matrix->cols, double);
    if (!row)
    {
        error("Matrix handling: Failed to allocate memory for row.\n");
    }

    // Copy the elements of the row
    for (int j = 0; j < matrix->cols; j++)
    {
        row[j] = MATRIX_AT_PTR(matrix, rowIndex, j);
    }

    return row;
}

/**
 * @brief Extracts the n-th column of a matrix as a dynamically allocated array.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @param[in] colIndex The index of the column to extract (0-based).
 * @return double* A dynamically allocated array containing the column elements.
 *
 * @note The caller is responsible for freeing the returned array.
 */
double *getColumn(const Matrix *matrix, int colIndex)
{
    checkMatrix(matrix); // Ensure the matrix is valid

    if (colIndex < 0 || colIndex >= matrix->cols)
    {
        error("Matrix handling: Column index out of bounds: %d\n", colIndex);
    }

    // Allocate memory for the column
    double *column = (double *)Calloc(matrix->rows, double);
    if (!column)
    {
        error("Matrix handling: Failed to allocate memory for column.\n");
    }

    // Copy the elements of the column
    for (int i = 0; i < matrix->rows; i++)
    {
        column[i] = MATRIX_AT_PTR(matrix, i, colIndex);
    }

    return column;
}

/**
 * @brief Adds a new row to a given matrix by reallocating its memory.
 *
 * The new row is appended to the end of the matrix. The function modifies the
 * input matrix struct in-place.
 *
 * @param[in, out] matrix Pointer to the matrix to which the row will be added.
 * @param[in] newRow Pointer to the array containing the elements of the new row.
 *
 * @return void
 *
 * @note
 * - The newRow must have the same number of elements as the matrix's columns.
 * - The matrix must be valid and well-defined.
 *
 * @example
 * Example usage:
 * @code
 * Matrix m = createMatrix(2, 3); // 2x3 matrix
 * double newRow[3] = {4.0, 5.0, 6.0};
 * addRowToMatrix(&m, newRow);
 * // m is now a 3x3 matrix
 * @endcode
 */
void addRowToMatrix(Matrix *matrix, const double *newRow)
{
    checkMatrix(matrix); // Ensure the matrix is valid

    if (!newRow)
    {
        error("Matrix handling: The new row pointer is NULL on the function to add a new row to the matrix.\n");
    }

    // Reallocate memory for the new row
    Matrix temp = copyMatrix(matrix);
    size_t newSize = (matrix->rows + 1) * matrix->cols;
    double *newData = Realloc(matrix->data, newSize, double);

    if (!newData)
    {
        error("Matrix handling: Failed to reallocate memory for the matrix.\n");
    }
    matrix->data = newData;

    // Shift existing columns down to make space for the new row
    for (int j = 0; j < matrix->cols; j++)
    {
        matrix->data[j * (matrix->rows + 1) + matrix->rows] = newRow[j];
    }

    for (int i = 0; i < matrix->rows + 1; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            if (i != matrix->rows)
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(temp, i, j);
            else
                MATRIX_AT_PTR(matrix, i, j) = newRow[j];
        }
    }

    // Update the matrix dimensions
    matrix->rows++;
    freeMatrix(&temp);
}

/**
 * @brief Removes a specific row from a matrix in place.
 *
 * This function modifies the input matrix to remove the specified row.
 *
 * @param[in,out] matrix Pointer to the matrix to modify.
 * @param[in] rowIndex The index of the row to remove (0-based).
 */
void removeRow(Matrix *matrix, int rowIndex)
{
    checkMatrix(matrix); // Validate the input matrix

    if (rowIndex < 0 || rowIndex >= matrix->rows)
    {
        error("Matrix handling: Row index out of bounds: %d\n", rowIndex);
    }

    // Shift rows up to overwrite the specified row
    Matrix temp = copyMatrix(matrix);
    matrix->rows -= 1;
    matrix->data = Realloc(matrix->data, matrix->rows * matrix->cols, double);

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            if (i < rowIndex)
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(temp, i, j);
            else if (i > rowIndex)
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(temp, i + 1, j);
            else
                continue;
        }
    }

    freeMatrix(&temp);
}

/**
 * @brief Adds a row of zeros at a specific index in a matrix in place.
 *
 * This function modifies the input matrix to add a row of zeros at the specified index.
 *
 * @param[in,out] matrix Pointer to the matrix to modify.
 * @param[in] rowIndex The index where the new row should be added (0-based).
 */
void addRowOfZeros(Matrix *matrix, int rowIndex)
{
    checkMatrix(matrix); // Validate the input matrix

    if (rowIndex < 0 || rowIndex > matrix->rows)
    {
        error("Matrix handling: Row index out of bounds: %d\n", rowIndex);
    }

    // Resize the matrix to have one additional row
    Matrix temp = copyMatrix(matrix);
    matrix->rows += 1;
    matrix->data = Realloc(matrix->data, matrix->rows * matrix->cols, double);
    if (!matrix->data)
    {
        error("Matrix handling: Memory reallocation failed while resizing the matrix.\n");
    }

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            if (i < rowIndex)
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(temp, i, j);
            else if (i > rowIndex)
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(temp, i - 1, j);
            else
                MATRIX_AT_PTR(matrix, i, j) = 0.0;
        }
    }
    freeMatrix(&temp);
}

/**
 * @brief Removes a specific column from a matrix in place.
 *
 * This function modifies the input matrix to remove the specified column.
 *
 * @param[in,out] matrix Pointer to the matrix to modify.
 * @param[in] colIndex The index of the column to remove (0-based).
 */
void removeColumn(Matrix *matrix, int colIndex)
{
    checkMatrix(matrix); // Validate the input matrix

    if (colIndex < 0 || colIndex >= matrix->cols)
    {
        error("Matrix handling: Column index out of bounds: %d\n", colIndex);
    }

    // Shift columns left to overwrite the specified column
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = colIndex; j < matrix->cols - 1; j++)
        {
            MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT_PTR(matrix, i, j + 1);
        }
    }

    // Resize the matrix to have one less column
    matrix->cols -= 1;
    matrix->data = Realloc(matrix->data, matrix->rows * matrix->cols, double);
    if (!matrix->data)
    {
        error("Matrix handling: Memory reallocation failed while resizing the matrix.\n");
    }
}

/**
 * @brief Adds a column of zeros at a specific index in a matrix in place.
 *
 * This function modifies the input matrix to add a column of zeros at the specified index.
 *
 * @param[in,out] matrix Pointer to the matrix to modify.
 * @param[in] colIndex The index where the new column should be added (0-based).
 */
void addColumnOfZeros(Matrix *matrix, int colIndex)
{
    checkMatrix(matrix); // Validate the input matrix

    if (colIndex < 0 || colIndex > matrix->cols)
    {
        error("Matrix handling: Column index out of bounds: %d\n", colIndex);
    }

    // Resize the matrix to have one additional column
    matrix->cols += 1;
    matrix->data = Realloc(matrix->data, matrix->rows * matrix->cols, double);
    if (!matrix->data)
    {
        error("Matrix handling: Memory reallocation failed while resizing the matrix.\n");
    }

    // Shift existing columns right to make space for the new column
    for (int j = matrix->cols - 1; j > colIndex; j--)
    {
        for (int i = 0; i < matrix->rows; i++)
        {
            MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT_PTR(matrix, i, j - 1);
        }
    }

    // Fill the new column with zeros
    for (int i = 0; i < matrix->rows; i++)
    {
        MATRIX_AT_PTR(matrix, i, colIndex) = 0.0;
    }
}

/**
 * @brief Creates a copy of the given Matrix.
 *
 * @param orig Pointer to the original Matrix.
 * @return Pointer to a new Matrix that is a copy of orig.
 *
 * This function uses malloc to allocate memory for both the Matrix struct and its data array.
 * The caller is responsible for freeing the memory (using free) when it is no longer needed.
 */
Matrix *copyMatrixPtr(const Matrix *orig)
{
    // Allocate memory for the new Matrix structure.
    Matrix *copy = (Matrix *)Calloc(1, Matrix);
    if (copy == NULL)
    {
        error("Memory allocation error in copyMatrix: could not allocate Matrix struct");
    }

    // Copy the dimensions.
    copy->rows = orig->rows;
    copy->cols = orig->cols;

    // Compute the total number of elements.
    int totalElements = orig->rows * orig->cols;

    // Allocate memory for the data array.
    copy->data = (double *)Calloc(totalElements, double);
    if (copy->data == NULL)
    {
        free(copy);
        error("Memory allocation error in copyMatrix: could not allocate data array");
    }

    // Copy the data from the original matrix.
    memcpy(copy->data, orig->data, totalElements * sizeof(double));

    return copy;
}

/*
 * Given an array of actions, it merges columns by summing the row values, generating a new matrix.
 * For example, if boundaries = {2, 4, 8} and wmat has 8 columns,
 * the function will merge columns as follows:
 * - New column 0: sum of columns 0 to 2
 * - New column 1: sum of columns 3 to 4
 * - New column 2: sum of columns 5 to 7
 *
 * @param[in] wmat A pointer to the original matrix. Won't be changed
 * @param[in] boundaries An array with the indices of boundaries.
 * @param[in] numBoundaries The size of the 'boundaries' array.
 *
 * @return A new matrix merged by columns
 */
Matrix mergeColumns(const Matrix *wmat, const int *boundaries, int numBoundaries)
{
    int newCols = numBoundaries; // Total merged column groups
    int rows = wmat->rows;
    Matrix newMat = createMatrix(rows, newCols);

    int start = 0; // First group starts from column 0

    for (int g = 0; g < newCols; g++)
    {
        int end = (g < numBoundaries) ? boundaries[g] : wmat->cols; // Adjust end index

        for (int r = 0; r < rows; r++)
        {
            double sum = 0.0;
            for (int c = start; c <= end; c++) // Ensure correct summation range
            {
                sum += MATRIX_AT_PTR(wmat, r, c);
            }
            MATRIX_AT(newMat, r, g) = sum;
        }

        start = end + 1; // Move start to the next segment
    }

    return newMat;
}

/*
 * @brief Checks if two matrices are equal
 *
 * @param[in] The first matrix to check
 * @param[in] The second matrix to check
 */
bool matricesAreEqual(Matrix *a, Matrix *b)
{
    checkMatrix(a);
    checkMatrix(b);
    for (int g = 0; g < TOTAL_GROUPS; g++)
    {
        for (int c = 0; c < TOTAL_CANDIDATES; c++)
        {
            if (MATRIX_AT_PTR(a, g, c) != MATRIX_AT_PTR(b, g, c))
                return false;
        }
    }
    return true;
}
