#ifndef FILEUTILS_H
#define FILEUTILS_H

#include "matrixUtils.h"

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

void readJSONAndStoreMatrices(const char *filename, Matrix *w, Matrix *x, Matrix *p);

/**
 * @brief Write the results to a .txt file.
 *
 * This function writes the name of the input JSON file and the contents of the
 * matrix `p` to a specified output text file.
 *
 * @param outputFileName Name of the output .txt file.
 * @param inputFileName Name of the input .json file.
 * @param *methodUsed The method that was used for the EM algorithm.
 * @param convergence The convergence value
 * @param iterations The amount of iterations occupied until the convergence
 * @param *pReal The ground truth probability matrix
 * @param *pCalculated The estimated probability matrix
 */
void writeResults(const char *outputFileName, const char *inputFileName, const char *methodUsed,
                  const double convergence, const double iterations, const double time, const int iterationsMade,
                  const Matrix *pReal, const Matrix *pCalculated, int S, int M, bool hit);

/**
 * @brief Write the results to a .json file.
 *
 * This function writes the name of the input JSON file, the method used, parameters, and the matrices
 * `pReal` and `pCalculated` to a specified output JSON file.
 *
 * @param outputFileName Name of the output .json file.
 * @param inputFileName Name of the input .json file.
 * @param methodUsed Name of the method used for calculations.
 * @param convergence Convergence threshold used.
 * @param iterations Maximum number of iterations allowed.
 * @param time Time taken to compute the results.
 * @param iterationsMade Actual number of iterations performed.
 * @param pReal Ground truth probability matrix.
 * @param pCalculated Estimated probability matrix.
 * @param S Samples (used in Hit and Run).
 * @param M Step size (used in Hit and Run).
 * @param hit Whether the method used was Hit and Run.
 */
void writeResultsJSON(const char *outputFileName, const char *inputFileName, const char *methodUsed,
                      const double convergence, const double iterations, const double time, const int iterationsMade,
                      const Matrix *pReal, const Matrix *pCalculated, const double *logLL, int S, int M, bool hit);

#endif
