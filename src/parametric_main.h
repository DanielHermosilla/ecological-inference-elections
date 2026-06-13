#ifndef PARAMETRIC_MAIN_H_EIM
#define PARAMETRIC_MAIN_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_matrix.h"

    typedef struct
    {
        Matrix *probabilities;    // length B, each G x C aggregated over components
        Matrix *q;                // length B, each G x C aggregated over components
        Matrix *expected;         // length B, each G x C aggregated over components
        Matrix beta;              // A x (K-1), coefficients for component membership probabilities
        Matrix phi;               // B x K, prior component probabilities from V and beta
        Matrix responsibilities;  // B x K
        double *component_prob;   // G x C x K
        int mixture_h;
        int total_iterations;
        double total_time;
        double log_likelihood;
        int finish_id;
    } EMParametricMixtureResult;

    Matrix *EM_Algorithm(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                         const double maxtime, const double ll_threshold, const int maxnewton, const bool verbose,
                         double *out_elapsed, int *total_iterations, double *logLikelihood,
                         Matrix **out_q, Matrix **out_expected,
                         const char *adjust_prob_cond_method, bool adjust_prob_cond_every);
    Matrix *EM_Algorithm_Method(Matrix *X, Matrix *W, Matrix *V, Matrix *beta, Matrix *alpha, const int maxiter,
                                const double maxtime, const double ll_threshold, const int maxnewton,
                                const bool verbose, double *out_elapsed, int *total_iterations,
                                double *logLikelihood, Matrix **out_q, Matrix **out_expected,
                                const char *q_method, QMethodInput *q_params,
                                const char *adjust_prob_cond_method, bool adjust_prob_cond_every);
    EMParametricMixtureResult *EM_Algorithm_Parametric_Mixture(
        Matrix *X, Matrix *W, Matrix *V, Matrix *component_prob_init, Matrix *membership_beta_init, int mixture_h,
        const int maxiter, const int miniter, const double maxtime, const double convergence, const double ll_threshold,
        const int maxnewton, const bool verbose, const char *q_method, QMethodInput *q_params,
        const char *adjust_prob_cond_method, bool adjust_prob_cond_every);
    void cleanupParametricMixtureResult(EMParametricMixtureResult *res);
#ifdef __cplusplus
}
#endif
#endif // PARAMETRIC_MAIN_H_EIM
