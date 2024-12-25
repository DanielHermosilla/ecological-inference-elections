#include "matrixUtils.h"
#include <omp.h> // Parallelization
#include <stdio.h>

/**
 * @brief Creates an empty matrix of given dimensions.
 *
 * Given certain dimensions of rows and colums, creates an empty Matrix with allocated memory towards the data.
 *
 * @param[in] rows The number of rows of the new matrix.
 * @param[in] cols The number of columns of the new matrix.
 *
 * @return Matrix Empty matrix of dimensions (rows x cols) with allocated memory for its data.
 *
 * @warning
 * - The memory may be full.
 */

Matrix createMatrix(int rows, int cols)
{
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double *)malloc(rows * cols * sizeof(double));
    if (!m.data)
    {
        perror("Failed to allocate matrix data");
        exit(EXIT_FAILURE);
    }
    return m;
}
/**
 * @brief Computes a row-wise sum.
 *
 * Given a matrix, it computes the sum over all the rows and stores them in an array.
 *
 * @param[in] matrix Pointer to an array that represents a (rows x col) matrix.
 * @param[out] result Pointer of the resulting array of length `rows`.
 * @param[in] rows The number of rows of the matrix.
 * @param[in] cols The number of columns of the matrix.
 *
 * @return void Written on *result
 *
 * @note
 * - Matrix should be in row-major order
 * - The size of the array **cannot** be checked. Be careful with this.
 *
 * @example
 * Example usage:
 * @code
 * double matrix[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 * double row_sums[2];
 * rowSum(matrix, row_sums, 2, 3);
 * // row_sums now contains [6.0, 15.0]
 * @endcode
 */
void rowSum(double *matrix, double *result, int rows, int cols)
{
// For parallelization, note that the array has its own identifier on the loop, hence,
// there shouldn't be a critical/coliding problem. The thread must be only made on "j"
// though, otherwise, there would be colissions with "i".
#pragma omp parallel for
    for (int i = 0; i < rows; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++)
        {
            result[i] += matrix[i * cols + j];
        }
    }
}

/**
 * @brief Computes a column-wise sum.
 *
 * Given a matrix, it computes the sum over all the columns and stores them in an array.
 *
 * @param[in] Matrix Pointer to a Matrix structure.
 * @param[out] Matrix Pointer of the resulting Matrix. It must have only 1 row.
 *
 * @return void Written on *result
 *
 * @note
 * - Matrix should be in row-major order.
 * - The resulting matrix must be of dimension 1 x `cols`
 *
 * @example
 * Example usage:
 * @code
 * double matrix[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 * double col_sum[3];
 * colSum(matrix, col_sum, 2, 3);
 * // col_sum now contains [5.0, 7.0, 9.0]
 * @endcode
 */

void colSum(Matrix *matrix, Matrix *result)
{
    // Validation, checks NULL pointer
    if (!matrix || !matrix->data || !result || !result->data)
    {
        fprintf(stderr, "A NULL pointer was handed to colSum.\n");
        return;
    }

    // Checks dimensions; result matrix must be of one row.
    if (result->rows != 1 || result->cols != matrix->cols)
    {
        fprintf(stderr, "The resulting matrix must have dimensions (1 x %d) for computing colSum.\n", matrix->cols);
        return;
    }

    // For parallelization, note that the array has its own identifier on the loop, hence,
    // there shouldn't be a critical/coliding problem. The thread must be only made on "j"
    // though, otherwise, there would be colissions with "i".

#pragma omp parallel for
    for (int j = 0; j < matrix->cols; j++)
    {
        result[j] = 0.0;
        for (int i = 0; i < matrix->rows; i++)
        {
            result->data[j] += matrix->data[i * matrix->cols + j];
        }
    }
}

/**
 * @brief Fills matrix with a constant value.
 *
 * Given a matrix, it fills a whole matrix with a constant value.
 * *
 * @param[in, out] matrix Pointer to an array that represents a (rows x col) matrix.
 * @param[in] rows The number of rows of the matrix.
 * @param[in] cols The number of columns of the matrix.
 * @param[in] value The constant value to fill
 *
 * @return void Written on the input matrix
 *
 * @note
 * - Matrix should be in row-major order.
 * - The size of the array **cannot** be checked. Be careful with this.
 *
 * @example
 * Example usage:
 * @code
 * double matrix[6] = {
 *     1.0, 2.0, 3.0,
 *     4.0, 5.0, 6.0
 * };
 *
 * fillMatrix(matrix, 2, 3, 9);
 * // matrix now contains [9.0, 9.0, 9.0, ..., 9.0]
 * @endcode
 */
void fillMatrix(double *matrix, int rows, int cols, double value)
{
#pragma omp parallel for
    for (int i = 0; i < rows * cols; i++)
    {
        matrix[i] = value;
    }
}
