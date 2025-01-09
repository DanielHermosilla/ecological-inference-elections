#ifndef COMPUTE_MULTIVARIATE_PDF_H
#define COMPUTE_MULTIVARIATE_PDF_H

#include "globals.h"
#include "matrixUtils.h"
#include "multivariateUtils.h"

/**
 * @brief Computes an approximate conditional probability using a Multivariate PDF approach.
 *
 * @param[in] probabilities    Matrix (g x c) - probabilities of each group and candidate.
 *
 * @return A pointer to a flattened 3D array (b x g x c) representing the probabilities.
 */
double *computeQMultivariatePDF(Matrix const *probabilities);

#endif
