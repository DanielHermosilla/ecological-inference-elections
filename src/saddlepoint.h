#ifndef COMPUTE_SADDLEPOINT_H_EIM
#define COMPUTE_SADDLEPOINT_H_EIM

#ifdef __cplusplus
extern "C"
{
#endif

#include "globals.h"

    /**
     * @brief Computes an approximate conditional probability using a saddlepoint PMF approximation.
     *
     * It approximates the PMF of X_b = sum_g Z_bg, with Z_bg ~ Multinomial(n_bg, p_g),
     * through a multivariate saddlepoint approximation, and derives q_{bgc} from the implied tilted probabilities.
     */
    void computeQSaddlepoint(EMContext *ctx, QMethodInput params, double *ll);

#ifdef __cplusplus
}
#endif

#endif
