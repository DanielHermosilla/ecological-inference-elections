#ifndef COMPUTE_EXACT_H
#define COMPUTE_EXACT_H

#include "globals.h"
#include "matrixUtils.h"
#include "memoizationUtil.h"

/**
 * @brief Calculate the value of `q_{bgc}`.
 *
 * It calculates all of the values for q_{bgc} by the definition on the paper. It returns the array of type `double`.
 *
 * @param[in] *probabilities A pointer to the matrix with the probabilities.
 *
 * @return *double: A pointer toward the array.
 *
 * @note: A single pointer is used to store the array continously. This is for using cBLAS operations later.
 *
 */

double *computeQExact(const Matrix *probabilities);

#endif
