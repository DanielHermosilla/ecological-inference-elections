#include "matrixUtils.h"
#include <omp.h> // Parallelization
#include <stdio.h>

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
 * @param[in] matrix Pointer to an array that represents a (rows x col) matrix.
 * @param[out] result Pointer of the resulting array of length `col`.
 * @param[in] rows The number of rows of the matrix.
 * @param[in] cols The number of columns of the matrix.
 *
 * @return void Written on *result
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
 * double col_sum[3];
 * colSum(matrix, col_sum, 2, 3);
 * // col_sum now contains [5.0, 7.0, 9.0]
 * @endcode
 */

void colSum(double *matrix, double *result, int rows, int cols)
{
// For parallelization, note that the array has its own identifier on the loop, hence,
// there shouldn't be a critical/coliding problem. The thread must be only made on "j"
// though, otherwise, there would be colissions with "i".
#pragma omp parallel for
    for (int j = 0; j < cols; j++)
    {
        result[j] = 0.0;
        for (int i = 0; i < rows; i++)
        {
            result[j] += matrix[i * cols + j];
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
