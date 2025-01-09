#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <stdbool.h>
#include <stdlib.h>

// Macro for making an easier indexation.
#define MATRIX_AT(matrix, i, j) (matrix.data[(i) * (matrix.cols) + (j)])
#define MATRIX_AT_PTR(matrix, i, j) (matrix->data[(i) * (matrix->cols) + (j)])

// All of the helper functions are made towards double type matrices
typedef struct
{
    double *data; // Pointer to matrix data array (row-major order)
    int rows;     // Number of rows
    int cols;     // Number of columns
} Matrix;

// The helper functions won't work towards this matrix
typedef struct
{
    size_t *data; // Pointer to matrix data array (row-major order)
    int rows;     // Number of rows
    int cols;     // Number of columns
} SizeTMatrix;

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

void checkMatrix(const Matrix *m);

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

void printMatrix(const Matrix *m);

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

void rowSum(const Matrix *matrix, double *result);

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

void colSum(const Matrix *matrix, double *result);

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

void fillMatrix(Matrix *matrix, const double value);

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
bool convergeMatrix(const Matrix *matrixA, const Matrix *matrixB, const double convergence);

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

double maxElement(const Matrix *m);

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
Matrix removeLastColumn(const Matrix *matrix);

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
Matrix removeLastRow(const Matrix *matrix);

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
Matrix createDiagonalMatrix(const double *vector, int size);

/**
 * @brief Computes the inverse of a symmetric positive-definite matrix using Cholesky decomposition.
 *
 * @param[in, out] matrix Pointer to the input matrix (overwritten with the inverse).
 *
 * @note The input matrix must be square and symmetric.
 */
void inverseSymmetricPositiveMatrix(Matrix *matrix);

/**
 * @brief Computes the inverse of a general square matrix using LU decomposition.
 *
 * @param[in, out] matrix Pointer to the input matrix (overwritten with the inverse).
 *
 * @note The input matrix must be square and invertible.
 */
void inverseMatrixLU(Matrix *matrix);

Matrix copyMatrix(const Matrix *original);

/**
 * @brief Extracts the n-th row of a matrix as a dynamically allocated array.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @param[in] rowIndex The index of the row to extract (0-based).
 * @return double* A dynamically allocated array containing the row elements.
 *
 * @note The caller is responsible for freeing the returned array.
 */
double *getRow(const Matrix *matrix, int rowIndex);

/**
 * @brief Extracts the n-th column of a matrix as a dynamically allocated array.
 *
 * @param[in] matrix Pointer to the input matrix.
 * @param[in] colIndex The index of the column to extract (0-based).
 * @return double* A dynamically allocated array containing the column elements.
 *
 * @note The caller is responsible for freeing the returned array.
 */
double *getColumn(const Matrix *matrix, int colIndex);

#endif // MATRIX_UTILS_H
