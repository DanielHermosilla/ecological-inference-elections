#ifndef DYNAMIC_H
#define DYNAMIC_H

#ifdef __cplusplus

extern "C"
{
#endif
#include "bootstrap.h"
#include "globals.h"
#include "utils_matrix.h"
    /*
     * Main function to obtain the heuristical best group aggregation, using dynamic programming. Tries every
     * combination using a standard deviation approximate. Given the approximate, computes the bootstrapped standard
     * deviation and checks if it accomplishes the proposed threshold.
     *
     * @param[in] xmat The candidate (c x b) matrix.
     * @param[in] wmat The group (b x g) matrix.
     * @param[in, out] results An array with the slicing indices.
     * @param[in] set_threshold The threshold of the proposed method
     * @param[in] set_method The method for evaluating the bootstrapping threshold.
     * @param[in] bootiter The amount of bootstrap iterations
     * @param[in] p_method The method for calculating the initial probability.
     * @param[in] q_method The method for calculating the EM algorithm of the boot samples.
     * @param[in] convergence The convegence threshold for the EM algorithm.
     * @param[in] maxIter The maximum amount of iterations to perform on the EM algorithm.
     * @param[in] maxSeconds The maximum amount of seconds to run the algorithm.
     * @param[in] verbose Boolean to whether verbose useful outputs.
     * @param[in] inputParams The parameters for specific methods.
     *
     * @return The heuristic optimal matrix with bootstrapped standard deviations.
     */
    Matrix aggregateGroups(const Matrix *xmat, const Matrix *wmat, int *results, int *cuts, double set_threshold,
                           const char *set_method, int bootiter, const char *p_method, const char *q_method,
                           const double convergence, const int maxIter, double maxSeconds, const bool verbose,
                           QMethodInput inputParams);

#ifdef __cplusplus
}
#endif
#endif
