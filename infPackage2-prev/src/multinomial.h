#ifndef COMPUTE_MULTINOMIAL_H_EIM
#define COMPUTE_MULTINOMIAL_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_matrix.h"

    /**
     * @brief Computes an approximate of the conditional probability by using a Multinomial approach.
     *
     * Given the observables parameters and the probability matrix, it computes an approximation of `q` with the
     * Multinomial approach.
     *
     * @param[in] *probabilities Matrix of dimension (gxc) with the probabilities of each group and candidate.
     * @param[in] params The optional parameters to the function. For this specific method, there's no supported
     * optional parameters yet.
     *
     * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices,
     * but it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
     *
     */
    double *computeQMultinomial(Matrix const *probabilities, QMethodInput params);
#ifdef __cplusplus
}
#endif
#endif
