#include "matrixUtils.h"
#include <cblas.h>
#include <math.h>
#include <omp.h> // Parallelization
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

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

void checkMatrix(Matrix *m)
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
 * @brief Creates a 4-dimensional matrix.
 *
 * @param[in] dim1 The size of the first dimension.
 * @param[in] dim2 The size of the second dimension.
 * @param[in] dim3 The size of the third dimension.
 * @param[in] dim4 Th size of the fourth dimension.
 *
 * @return int Set of pointers to each dimensional array. This matrix is not from the Matrix struct.
 *
 */

double ****create4DMatrix(int dim1, int dim2, int dim3, int dim4)
{
    double ****matrix = calloc(dim1, sizeof(double ***));
    for (int i = 0; i < dim1; i++)
    {
        matrix[i] = calloc(dim2, sizeof(double **));
        for (int j = 0; j < dim2; j++)
        {
            matrix[i][j] = calloc(dim3, sizeof(double *));
            for (int k = 0; k < dim3; k++)
            {
                matrix[i][j][k] = calloc(dim4, sizeof(double));
            }
        }
    }
    return matrix;
}

/**
 * @brief Liberates the allocated 4D-matrix data.
 *
 * @param[in] ****m The matrix pointer to free the data.
 * @param[in] dim1 The size of the first dimension.
 * @param[in] dim2 The size of the second dimension.
 * @param[in] dim3 The size of the third dimension.
 *
 * @return void Changes to be made on the input matrix and memory.
 *
 */

void free4DMatrix(int ****matrix, int dim1, int dim2, int dim3)
{
    for (int i = 0; i < dim1; i++)
    {
        for (int j = 0; j < dim2; j++)
        {
            for (int k = 0; k < dim3; k++)
            {
                free(matrix[i][j][k]);
            }
            free(matrix[i][j]);
        }
        free(matrix[i]);
    }
    free(matrix);
}

/**
 * @brief Liberates the allocated 4D-matrix data.
 *
 * @param[in] ****m The matrix pointer to free the data.
 * @param[in] dim1 The size of the first dimension.
 * @param[in] dim2 The size of the second dimension.
 * @param[in] dim3 The size of the third dimension.
 *
 * @return void Changes to be made on the input matrix and memory.
 *
 */

void free4DMatrix(int ****matrix, int dim1, int dim2, int dim3);

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

void printMatrix(Matrix *matrix)
{
    checkMatrix(matrix); // Assertion

    printf("Matrix (%dx%d) of type double\n", matrix->rows, matrix->cols);

    for (int i = 0; i < matrix->rows; i++)
    {
        for (int j = 0; j < matrix->cols; j++)
        {
            printf("%.3f ", matrix->data[(i) * (matrix->cols) + (j)]);
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

void rowSum(Matrix *matrix, double *result)
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

void colSum(Matrix *matrix, double *result)
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

void fillMatrix(Matrix *matrix, double value)
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
bool convergeMatrix(Matrix *matrixA, Matrix *matrixB, double convergence)
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
