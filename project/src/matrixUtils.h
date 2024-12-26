#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    double *data; // Pointer to matrix data array (row-major order)
    int rows;     // Number of rows
    int cols;     // Number of columns
} Matrix;

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

Matrix createMatrix(int rows, int cols);

/**
 * @brief Liberates the allocated matrix data.
 *
 * @param[in] m The matrix to free the data.
 *
 * @return void Changes to be made on the input matrix and memory.
 *
 */

void freeMatrix(Matrix *m);

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

void printMatrix(Matrix *m);

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

void rowSum(Matrix *matrix, double *result);

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

void colSum(Matrix *matrix, double *result);

/**
 * @brief Make an array of a constant value.
 *
 * Given a value, it fills a whole array with a constant value.
 *
 * @param[in, out] array Pointer to matrix to be filled.
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

void makeArray(double *array, int N, double value);

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
void fillMatrix(Matrix *matrix, double value);

#endif // MATRIX_UTILS_H
