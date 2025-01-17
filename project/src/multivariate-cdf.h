#ifndef COMPUTE_MULTIVARIATE_CDF_H
#define COMPUTE_MULTIVARIATE_CDF_H

#include "globals.h"
#include "matrixUtils.h"
#include "multivariateUtils.h"

/**
 * @brief Computes an approximate conditional probability using a Multivariate CDF approach.
 *
 * @param[in] probabilities Matrix (g x c) - probabilities of each group and candidate.
 * @param[in] monteCarloSamples Amount of samples to use in the Montecarlo simulation
 * @param[in] epsilon The error threshold used for the Genz Montecarlo method.
 * @param[in] *method The method for calculating the Montecarlo simulation. Currently available methods are `Plain`,
 * `Miser` and `Vegas`.

 * @return A pointer to a flattened 3D array (b x g x c) representing the probabilities.
 */
double *computeQMultivariateCDF(Matrix const *probabilities, int monteCarloSamples, double epsilon, const char *method);

#endif
