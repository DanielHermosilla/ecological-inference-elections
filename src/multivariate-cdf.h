#ifndef COMPUTE_MULTIVARIATE_CDF_H_EIM
#define COMPUTE_MULTIVARIATE_CDF_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_matrix.h"
#include "utils_multivariate.h"

    /**
     * @brief Computes an approximate conditional probability using a Multivariate CDF approach.
     *
     * @param[in] probabilities Matrix (g x c) - probabilities of each group and candidate.
     * @param[in] params A QMethodInput struct that has the Monte Carlo samples, epsilon (error threshold) and the
     method to run.
     * @param[in] monteCarloSamples Amount of samples to use in the Montecarlo simulation
     * @param[in] epsilon The error threshold used for the Genz Montecarlo method.
     * @param[in] *method The method for calculating the Montecarlo simulation. Currently available methods are `Plain`,
     * `Miser` and `Vegas`.

     * @return A pointer to a flattened 3D array (b x g x c) representing the probabilities.
     */
    void computeQMultivariateCDF(EMContext *ctx, QMethodInput params, double *ll);

#ifdef __cplusplus
}
#endif
#endif
