#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "main.h"
#include "utils_matrix.h"

    /**
     *  Returns an array of col-major matrices with bootstrapped matrices.
     *
     * @param[in] xmat The original X array
     * @param[in] wmat The original W array
     * @param[in] bootiter The amount of iterations for bootstrapping
     * @param[in] p_method The method for obtaining the initial probability
     * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
     * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
     * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
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
     *
     * @return An allocated array of size bootiter * TOTAL_BALLOTS that stores matrices.
     */
    Matrix bootstrapA(const Matrix *xmat, const Matrix *wmat, int bootiter, const char *q_method, const char *p_method,
                      const double convergence, const double log_convergence, const int maxIter,
                      const double maxSeconds, const bool verbose, QMethodInput *inputParams);

    Matrix bootSingleMat(Matrix *xmat, Matrix *wmat, int bootiter, const bool verbose);

#ifdef __cplusplus
}
#endif
#endif
