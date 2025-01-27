#include "matrixUtils.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <omp.h> // Parallelization
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
// Macro for easier matrix indexation
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])

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
        fprintf(stderr, "A NULL pointer was handed as an array.\n");
        exit(EXIT_FAILURE);
    }

    if (N < 0)
    {
        fprintf(stderr, "A incoherent dimension was handen for making the array.\n");
        exit(EXIT_FAILURE);
    }

// Fill the array with the specified constant value
#pragma omp parallel for
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
        fprintf(stderr, "A NULL pointer was handed as a matrix argument.\n");
        exit(EXIT_FAILURE);
    }

    // Checks dimensions
    if (m->rows <= 0 || m->cols <= 0)
    {
        fprintf(stderr, "Invalid matrix dimensions: rows=%d, cols=%d\n", m->rows, m->cols);
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Invalid matrix dimensions: rows=%d, cols=%d\n", rows, cols);
        exit(EXIT_FAILURE);
    }

    Matrix m;
    m.rows = rows;
    m.cols = cols;

    m.data = calloc(rows * cols, sizeof(double));

    if (!m.data)
    {
        fprintf(stderr, "Failed to allocate matrix data\n");
        exit(EXIT_FAILURE);
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
        free(m->data);
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

    printf("Matrix (%dx%d) of type double\n", matrix->rows, matrix->cols);

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            printf("%.5f ", matrix->data[(i) * (matrix->cols) + (j)]);
        }
        printf("\n");
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
 * - Matrix should be in row-major order
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

    double *ones = (double *)malloc(matrix->cols * sizeof(double));
    if (!ones)
    {
        fprintf(stderr, "Failed to allocate memory to the rowSum function. \n");
        exit(EXIT_FAILURE);
    }

    makeArray(ones, matrix->cols, 1.0);
    // Perform Matrix-Vector Multiplication (Matrix * Ones = Row Sums)
    cblas_dgemv(CblasRowMajor, // Row-major storage
                CblasNoTrans,  // Don't transpose the matrix
                matrix->rows,  // Number of rows
                matrix->cols,  // Number of columns
                1.0,           // Scalar multiplier => y = 1.0 *  (A * x) + beta * y
                matrix->data,  // Matrix pointer
                matrix->cols,  // Leading dimension (number of columns in case of row-major)
                ones,          // Vector of ones
                1,             // Increment for vector (1 = contiguous)
                0.0,           // Beta multiplier => y = alpha * (A*x) + 0.0 * y
                result,        // Output row sums
                1              // Step size for writing the results
    );

    free(ones);
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
 * - Matrix should be in row-major order
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

    double *ones = (double *)malloc(matrix->rows * sizeof(double));
    makeArray(ones, matrix->rows, 1.0);

    // Perform Matrix-Vector Multiplication (Matrix * Ones = Row Sums)
    cblas_dgemv(CblasRowMajor, // Row-major storage
                CblasTrans,    // Don't transpose the matrix
                matrix->rows,  // Number of rows
                matrix->cols,  // Number of columns
                1.0,           // Scalar multiplier => y = 1.0 *  (A * x) + beta * y
                matrix->data,  // Matrix pointer
                matrix->cols,  // Leading dimension (number of columns in case of row-major)
                ones,          // Vector of ones
                1,             // Increment for vector (1 = contiguous)
                0.0,           // Beta multiplier => y = alpha * (A*x) + 0.0 * y
                result,        // Output row sums
                1              // Step size for writing the results
    );
    free(ones);
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
 * - Matrix should be in row-major order.
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
 * - Matrix should be in row-major order.
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
        fprintf(stderr, "An invalid value was handed to convergence, it must be greater than cero.\n");
        exit(EXIT_FAILURE);
    }

    if (matrixA->cols != matrixB->cols || matrixA->rows != matrixB->rows)
    {
        fprintf(stderr, "The dimensions of both matrices doesn't match.\n");
        exit(EXIT_FAILURE);
    }

    int size = matrixA->rows * matrixB->cols;

    double *diff = (double *)malloc(size * sizeof(double));

    cblas_dcopy(size, matrixA->data, 1, diff, 1);
    cblas_daxpy(size, -1.0, matrixB->data, 1, diff, 1);

    for (int i = 0; i < size; i++)
    {
        // If there's a value whom convergence is greater than epsilon, the convergence
        // isn't achieved.
        if (fabs(diff[i]) >= convergence)
        {
            free(diff);
            return false;
        }
    }

    free(diff);
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
 * - Matrix should be in row-major order.
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
        fprintf(stderr, "Matrix must have at least two rows to remove one.\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Matrix must have at least two columns to remove one.\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Failed to allocate memory for diagonal matrix.\n");
        exit(EXIT_FAILURE);
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
    checkMatrix(matrix); // Ensure the matrix is valid

    if (matrix->rows != matrix->cols)
    {
        fprintf(stderr, "Matrix must be square for inversion.\n");
        exit(EXIT_FAILURE);
    }

    if (matrix->rows == 1 && matrix->cols == 1)
    {
        double currentVal = MATRIX_AT_PTR(matrix, 0, 0);
        if (currentVal != 0)
            MATRIX_AT_PTR(matrix, 0, 0) = 1 / currentVal;
        return;
    }

    int n = matrix->rows;
    int lda = n; // Leading dimension (number of columns in row-major storage)
    int info;
    Matrix emergencyMat = copyMatrix(matrix);
    // Cholesky Decomposition (L * L^T = A)
    info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, matrix->data, lda);
    if (info < 0)
    {
        fprintf(stderr,
                "Cholesky decomposition failed. The %d argument towards LAPACKE_dpotrf is had an illegal value.\n",
                info);
        exit(EXIT_FAILURE);
    }
    if (info > 0)
    {
        fprintf(stderr,
                "Cholesky decomposition failed. The leading minor of order %d is not positive definite.\nNote: if "
                "dealing with variance matrices of a Multivariate Normal, remember that there are selected cases where "
                "the Cholesky decomposition could fail in case of an eigenvalue being zero.\nRetrying by adding a "
                "small perturbation to the diagonals\n",
                info);

        for (int i = 0; i < matrix->rows; i++)
        {
            for (int j = 0; j < matrix->cols; j++)
            {
                MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT(emergencyMat, i, j);
                if (i == j)
                    MATRIX_AT_PTR(matrix, i, j) += 1;
            }
        }
        freeMatrix(&emergencyMat);
        inverseSymmetricPositiveMatrix(matrix);
    }

    // Invert the Cholesky Factorization

    info = LAPACKE_dpotri(LAPACK_ROW_MAJOR, 'L', n, matrix->data, lda);
    if (info != 0)
    {
        fprintf(stderr, "Matrix inversion failed after Cholesky decomposition. Error code: %d\n", info);
        printMatrix(matrix);
        exit(EXIT_FAILURE);
    }

    // Fill the upper triangle of the inverse matrix, this is not really necessary, but would prevent future problems.
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
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
 * @param[in,out] matrix Pointer to the NxN symmetric matrix in row-major layout.
 */
void inverseMatrixEigen(Matrix *matrix)
{
    checkMatrix(matrix);

    if (matrix->rows != matrix->cols)
    {
        fprintf(stderr, "inverseMatrixEigen: Matrix must be square.\n");
        exit(EXIT_FAILURE);
    }
    int n = matrix->rows;

    // Allocate space for eigenvalues
    double *eigenvals = (double *)malloc(n * sizeof(double));
    if (!eigenvals)
    {
        fprintf(stderr, "inverseMatrixEigen: cannot allocate eigenvals.\n");
        exit(EXIT_FAILURE);
    }

    // 2) dsyev: compute all eigenvalues and eigenvectors of a real symmetric matrix
    //    - 'V' means we want both eigenvalues and eigenvectors
    //    - 'U' means the matrix is stored in the upper part (row-major).
    //    On exit, matrix->data holds the eigenvectors in columns, eigenvals[] holds the eigenvalues.
    int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, // row-major storage
                             'V',              // compute Eigenvalues & Eigenvectors
                             'U',              // 'U' => input matrix is in the upper triangle
                             n,                // dimension
                             matrix->data,     // in/out: on exit, columns = eigenvectors
                             n,                // leading dimension (row-major)
                             eigenvals         // out: eigenvalues
    );

    if (info != 0)
    {
        fprintf(stderr, "inverseMatrixEigen: dsyev failed (info = %d)\n", info);
        free(eigenvals);
        exit(EXIT_FAILURE);
    }

    // Invert the eigenvalues => 1 / lambda_i (check none are zero too)
    for (int i = 0; i < n; i++)
    {
        if (fabs(eigenvals[i]) < 1e-15)
        {
            fprintf(stderr, "inverseMatrixEigen: Zero or near-zero eigenvalue => not invertible.\n");
            free(eigenvals);
            exit(EXIT_FAILURE);
        }
        eigenvals[i] = 1.0 / eigenvals[i];
    }

    // ---- Calculation of A^{-1} as Q *1/\lambda * Q^T ---- //
    // Build the diagonal matrix from eigenvals
    Matrix Dinv = createDiagonalMatrix(eigenvals, n);

    // temp = Q * Dinv
    Matrix temp = createMatrix(n, n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, // Q is NxN, Dinv is NxN
                n, n, n, 1.0, matrix->data, n,             // Q (eigenvectors in columns)
                Dinv.data, n, 0.0, temp.data, n);

    // A_inv = temp * Q^T
    // I will use a temporary matrix since it generate errors when recicling a variable
    Matrix temp2 = createMatrix(n, n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, temp.data, n, // first operand
                matrix->data, n,                                                     // second operand
                0.0, temp2.data, n);

    memcpy(matrix->data, temp2.data, n * n * sizeof(double));
    freeMatrix(&temp2);
    freeMatrix(&temp);
    freeMatrix(&Dinv);
    free(eigenvals);
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
    checkMatrix(matrix); // Validate the matrix

    if (matrix->rows != matrix->cols)
    {
        fprintf(stderr, "Matrix must be square for inversion.\n");
        exit(EXIT_FAILURE);
    }

    int n = matrix->rows;
    int *ipiv = malloc(n * sizeof(int)); // Pivot indices for LU decomposition
    if (!ipiv)
    {
        fprintf(stderr, "Failed to allocate memory for pivot indices.\n");
        exit(EXIT_FAILURE);
    }

    int info;

    // Perform LU decomposition (matrix is overwritten with LU factors)
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matrix->data, n, ipiv);
    if (info != 0)
    {
        fprintf(stderr, "LU decomposition failed. Error code: %d\n", info);
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    // Compute the inverse using the LU decomposition (matrix is overwritten with its inverse)
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, matrix->data, n, ipiv);
    if (info != 0)
    {
        fprintf(stderr, "Matrix inversion failed. Error code: %d\n", info);
        free(ipiv);
        exit(EXIT_FAILURE);
    }

    free(ipiv);
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
        fprintf(stderr, "Row index out of bounds: %d\n", rowIndex);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the row
    double *row = (double *)malloc(matrix->cols * sizeof(double));
    if (!row)
    {
        fprintf(stderr, "Failed to allocate memory for row.\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Column index out of bounds: %d\n", colIndex);
        exit(EXIT_FAILURE);
    }

    // Allocate memory for the column
    double *column = (double *)malloc(matrix->rows * sizeof(double));
    if (!column)
    {
        fprintf(stderr, "Failed to allocate memory for column.\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "addRowToMatrix: The new row pointer is NULL.\n");
        exit(EXIT_FAILURE);
    }

    // Reallocate memory for the new row
    size_t newSize = (matrix->rows + 1) * matrix->cols * sizeof(double);
    double *newData = realloc(matrix->data, newSize);

    if (!newData)
    {
        fprintf(stderr, "addRowToMatrix: Failed to reallocate memory for the matrix.\n");
        exit(EXIT_FAILURE);
    }

    matrix->data = newData;

    // Append the new row to the matrix
    memcpy(&matrix->data[matrix->rows * matrix->cols], newRow, matrix->cols * sizeof(double));

    // Update the matrix dimensions
    matrix->rows++;
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
        fprintf(stderr, "Row index out of bounds: %d\n", rowIndex);
        exit(EXIT_FAILURE);
    }

    // Shift rows up to overwrite the specified row
    for (int i = rowIndex; i < matrix->rows - 1; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT_PTR(matrix, i + 1, j);
        }
    }

    // Resize the matrix to have one less row
    matrix->rows -= 1;
    matrix->data = realloc(matrix->data, matrix->rows * matrix->cols * sizeof(double));
    if (!matrix->data)
    {
        fprintf(stderr, "Memory reallocation failed while resizing the matrix.\n");
        exit(EXIT_FAILURE);
    }
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
        fprintf(stderr, "Row index out of bounds: %d\n", rowIndex);
        exit(EXIT_FAILURE);
    }

    // Resize the matrix to have one additional row
    matrix->rows += 1;
    matrix->data = realloc(matrix->data, matrix->rows * matrix->cols * sizeof(double));
    if (!matrix->data)
    {
        fprintf(stderr, "Memory reallocation failed while resizing the matrix.\n");
        exit(EXIT_FAILURE);
    }

    // Shift rows down to make space for the new row
    for (int i = matrix->rows - 1; i > rowIndex; i--)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            MATRIX_AT_PTR(matrix, i, j) = MATRIX_AT_PTR(matrix, i - 1, j);
        }
    }

    // Fill the new row with zeros
    for (int j = 0; j < matrix->cols; j++)
    {
        MATRIX_AT_PTR(matrix, rowIndex, j) = 0.0;
    }
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
        fprintf(stderr, "Column index out of bounds: %d\n", colIndex);
        exit(EXIT_FAILURE);
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
    matrix->data = realloc(matrix->data, matrix->rows * matrix->cols * sizeof(double));
    if (!matrix->data)
    {
        fprintf(stderr, "Memory reallocation failed while resizing the matrix.\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Column index out of bounds: %d\n", colIndex);
        exit(EXIT_FAILURE);
    }

    // Resize the matrix to have one additional column
    matrix->cols += 1;
    matrix->data = realloc(matrix->data, matrix->rows * matrix->cols * sizeof(double));
    if (!matrix->data)
    {
        fprintf(stderr, "Memory reallocation failed while resizing the matrix.\n");
        exit(EXIT_FAILURE);
    }

    // Shift columns right to make space for the new column
    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = matrix->cols - 1; j > colIndex; j--)
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
