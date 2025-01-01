#ifndef FILEUTILS_H
#define FILEUTILS_H

#include "matrixUtils.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Writes an array of matrices to a binary file, including their metadata.
 *
 * Given a pointer to an array of matrices, it writes them to a binary file. Usually it's useful for saving the `X` and
 * `W` matrices with their probabilities matrix. It's also possible to save a single matrix.
 * @param filename Name of the binary file.
 * @param matrices Pointer to the array of matrices.
 * @param count Number of matrices in the array.
 */
void writeMatrices(const char *filename, Matrix *matrices, int count);

/**
 * @brief Reads an array of matrices from a binary file.
 *
 * @param filename Name of the binary file.
 * @param matrices Pointer to an array of matrices.
 * @param count Pointer to store the number of matrices read.
 */
void readMatrices(const char *filename, Matrix *matrices, int count);

#endif
