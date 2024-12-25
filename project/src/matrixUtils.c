#include "matrixUtils.h"
#include <cstdlib>
#include <omp.h> // Parallelization
#include <stdio.h>

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
    if (m->data)
    {
        free(m->data);
        m->data = NULL;
    }
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

void printMatrix(Matrix *m)
{

    // Don't use parallelization, otherwise, the prints would be messed up.
    for (int i = 0; i < m->rows; i++)
    {
        for (int j = 0; j < m->cols; j++)
        {
            printf("%.2f ", m->data[i * m->cols + j]);
        }
        printf("\n");
    }
}

/**
 * @brief Computes a row-wise sum.
 *
 * Given a matrix, it computes the sum over all the rows and stores them in an array.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @param[out] result Pointer of the resulting array of length `rows`.
 *
 * @return void Written on *result
 *
 * @note
 * - Matrix should be in row-major order
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
    // Validation, checks NULL pointer
    if (!matrix || !matrix->data || !result)
    {
        fprintf(stderr, "A NULL pointer was handed to rowSum.\n");
        exit(EXIT_FAILURE);
    }
    // For parallelization, note that the array has its own identifier on the loop, hence,
    // there shouldn't be a critical/coliding problem. The thread must be only made on "j"
    // though, otherwise, there would be colissions with "i".
    //
#pragma omp parallel for
    for (int i = 0; i < matrix->rows; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < matrix->cols; j++)
        {
            result[i] += matrix->data[i * matrix->cols + j];
        }
    }
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
 * .data = values,
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
    // Validation, checks NULL pointer
    if (!matrix || !matrix->data || !result)
    {
        fprintf(stderr, "A NULL pointer was handed to colSum.\n");
        exit(EXIT_FAILURE);
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
            result[j] += matrix->data[i * matrix->cols + j];
        }
    }
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
 * double matrix[6] = {
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
#pragma omp parallel for
    for (int i = 0; i < matrix->rows * matrix->cols; i++)
    {
        matrix->data[i] = value;
    }
}
