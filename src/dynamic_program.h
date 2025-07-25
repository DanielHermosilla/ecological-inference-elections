#ifndef DYNAMIC_H
#define DYNAMIC_H

#ifdef __cplusplus

extern "C"
{
#endif
#include "bootstrap.h"
#include "globals.h"
#include "main.h"
#include "utils_matrix.h"
    /*
     * Main function to obtain the heuristical best group aggregation, using dynamic programming. Tries every
     * combination using a standard deviation approximate. Given the approximate, computes the bootstrapped standard
     * deviation and checks if it accomplishes the proposed threshold.
     *
     * @param[in] xmat The candidate (c x b) matrix.
     * @param[in] wmat The group (b x g) matrix.
     * @param[in, out] results An array with the slicing indices.
     * @param[in, out] cuts The size of results
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
    Matrix aggregateGroups(
        // ---- Matrices
        EMContext *ctx,

        // ---- Results
        int *results, // Array with cutting indices
        int *cuts,    // Amount of cuts
        bool *bestResult,

        // ---- EM and Bootstrap parameters
        double set_threshold, const char *set_method, bool feasible, int bootiter, const char *p_method,
        const char *q_method, const double convergence, const double log_convergence, const int maxIter,
        double maxSeconds, const bool verbose, QMethodInput *inputParams);

    /*
     * Function for testing all of the 2^{G-1} combinations, returning the best aggregation according the
     * log-likelihood.
     *
     * @param[in] xmat The candidate (c x b) matrix.
     * @param[in] wmat The group (b x g) matrix.
     * @param[in, out results An array with the slicing indices.
     * @param[in, out] cuts The size of results
     * @param[in] p_method The method for calculating the initial probability.
     * @param[in] q_method The method for calculating the EM algorithm of the boot samples.
     * @param[in] convergence The convegence threshold for the EM algorithm.
     * @param[in] log_convergence The log-convergence of the EM.
     * @param[in] verbose Boolean to whether verbose useful outputs.
     * @param[in] maxIter The maximum amount of iterations to perform on the EM algorithm.
     * @param[in] maxSeconds The maximum amount of seconds to run the algorithm.
     * @param[in] inputParams The parameters for specific methods.
     * @param[in] outBestLL The best log-likelihood.
     * @param[in] outbestQ The best Q-value.
     * @param[in] outBestTime The best time
     * @param[in] outFinishReason The finishing reason
     * @param[in] outIterTotal The total amount of iterations of the EM
     *
     */
    Matrix aggregateGroupsExhaustive(Matrix *xmat, Matrix *wmat, int *results, int *cuts, const char *set_method,
                                     int bootiter, double max_qual, const char *p_method, const char *q_method,
                                     double convergence, double log_convergence, bool verbose, int maxIter,
                                     double maxSeconds, QMethodInput *inputParams, double *outBestLL, double **outBestQ,
                                     Matrix **bestBootstrap, double *outBestTime, int *outFinishReason,
                                     int *outIterTotal);
#ifdef __cplusplus
}
#endif
#endif
