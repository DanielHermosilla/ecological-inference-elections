#ifndef COMPUTE_MULTIVARIATE_PDF_H_EIM
#define COMPUTE_MULTIVARIATE_PDF_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_matrix.h"
#include "utils_multivariate.h"

    /**
     * @brief Computes an approximate conditional probability using a Multivariate PDF approach.
     *
     * @param[in] probabilities Matrix (g x c) - probabilities of each group and candidate.
     *
     * @return A pointer to a flattened 3D array (b x g x c) representing the probabilities.
     */
    double *computeQMultivariatePDF(Matrix const *probabilities, QMethodInput params, double *ll);

#ifdef __cplusplus
}
#endif
#endif
