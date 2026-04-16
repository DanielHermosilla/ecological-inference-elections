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
#include "saddlepoint.h"
#include "utils_matrix.h"

    typedef struct
    {
        Matrix probabilities;    // G x C
        double *q;               // B x G x C (flattened with Q_3D indexing)
        double *predicted_votes; // B x G x C
        double *phi;             // K
        Matrix responsibilities; // B x K
        double *component_prob;  // G x C x K
        Matrix probabilities_inv;    // C x G when symmetric EM_weight is used
        double *q_inv;               // B x C x G
        double *predicted_votes_inv; // B x C x G
        int mixture_h;
        int total_iterations;
        double total_time;
        double log_likelihood;
        int finish_id;
    } EMMixtureResult;

    // ---...--- //
    //
    EMContext *createEMContext(Matrix *X, Matrix *W, const char *method, QMethodInput params);
    void getInitialP(EMContext *ctx, const char *p_method, Matrix *probMatrix);
    QMethodConfig getQMethodConfig(const char *q_method, QMethodInput inputParams);
    void normalizeProbabilityRows(Matrix *P);
    void projectQ(EMContext *ctx, QMethodInput inputParams);
    void getP(EMContext *ctx);
    void getPredictedVotes(EMContext *ctx);
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
                          double *time, int *iterTotal, double *logLLarr, int *finishing_reason, Matrix *probMatrix,
                          QMethodInput *inputParams);

    EMMixtureResult *EMAlgoritmMixture(Matrix *X, Matrix *W, const char *p_method, const char *q_method,
                                       const double convergence, const double LLconvergence, const int maxIter,
                                       const double maxSeconds, const bool verbose, Matrix *probMatrix,
                                       QMethodInput *inputParams, int mixture_h);
    EMMixtureResult *EMAlgoritmRowMixture(Matrix *X, Matrix *W, const char *p_method, const char *q_method,
                                          const double convergence, const double LLconvergence, const int maxIter,
                                          const double maxSeconds, const bool verbose, Matrix *probMatrix,
                                          QMethodInput *inputParams, int row_mixture_h);

    bool hasMismatch(Matrix *X, Matrix *W);

    double computeLogLikForProbability(Matrix *X, Matrix *W, Matrix *probMatrix, const char *q_method,
                                       QMethodInput *inputParams);

    Matrix precomputeNorm(double *scale_factors, Matrix *X, Matrix *W);

    void precomputeScaleFactors(double *scale_factors, Matrix *X, Matrix *W);

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
    void cleanupMixtureResult(EMMixtureResult *res);
#ifdef __cplusplus
}
#endif
#endif // UTIL_H
