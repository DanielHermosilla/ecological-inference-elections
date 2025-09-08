#ifndef MAIN_H_EIM
#define MAIN_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "LP.h"
#include "MCMC.h"
#include "exact.h"
#include "globals.h"
#include "multinomial.h"
#include "multivariate-cdf.h"
#include "multivariate-pdf.h"
#include "utils_matrix.h"

    // ---...--- //
    //
    EMContext *createEMContext(Matrix *X, Matrix *W, const char *method, QMethodInput params);
    /**
     * @brief Implements the whole EM algorithm.
     *
     * Given a method for estimating "q", it calculates the EM until it converges to arbitrary parameters. As of in the
     * paper, it currently supports hnr, mult, mvn_cdf and mvn_pdf methods.
     *
     * @param[in] currentP Matrix of dimension (cxg) with the initial probabilities for the first iteration.
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
     * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
     * @param[in] LLconvergence Threshold regarding the convergence of the log-likelihood between iterations.
     * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
     * @param[in] maxSeconds Double with the value of the maximum amount of seconds to use.
     * @param[in] verbose Wether to verbose useful outputs.
     * @param[in, out] time The time that the algorithm took.
     * @param[in, out] iterTotal Total amount of iterations.
     * @param[in, out] logLLarr The loglikelihood array
     * @param[in, out] finishing_reason The reason that the algorithm has been stopped. It can either be 0, 1, 2, 3,
     * representing a normal convergence, log likelihood decrease, maximum time reached and maximum iterations reached,
     * respectively.
     *
     * @return Matrix: A matrix with the final probabilities. In case it doesn't converges, it returns the last
     * probability that was computed
     *
     * @note This is the main function that calls every other function for "q"
     *
     * @see getInitialP() for getting initial probabilities. group_proportional method is recommended.
     *
     * @warning
     * - Pointers shouldn't be NULL.
     * - `x` and `w` dimensions must be coherent.
     *
     */
    EMContext *EMAlgoritm(Matrix *X, Matrix *W, const char *p_method, const char *q_method, const double convergence,
                          const double LLconvergence, const int maxIter, const double maxSeconds, const bool verbose,
                          double *time, int *iterTotal, double *logLLarr, int *finishing_reason,
                          QMethodInput *inputParams);

    Matrix precomputeNorm(Matrix *W);

    /**
     * @brief Checks if a candidate didn't receive any votes.
     *
     * Given an array of size TOTAL_CANDIDATES, it sets to "1" the index where a possible candidate haven't received
     * any vote. It also returns a boolean indicating whether a candidate hasn't receive any vote
     *
     * @param[in,out] *canArray Array of size TOTAL_CANDIDATES full of zeroes, indicating with a "1" on the index
     * where a given candidate haven't received a vote
     *
     * @return bool: A boolean that shows if it exists a candidate with no votes
     *
     */
    void cleanup(EMContext *ctx);
#ifdef __cplusplus
}
#endif
#endif // UTIL_H
