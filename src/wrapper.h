#ifndef WRAPPER_H_EIM
#define WRAPPER_H_EIM

/* From CRAN guide to packages:
 *Macros defined by the compiler/OS can cause problems. Identifiers starting with an underscore followed by an
 *upper-case letter or another underscore are reserved for system macros and should not be used in portable code
 *(including not as guards in C/C++ headers). Other macros, typically upper-case, may be defined by the compiler or
 *system headers and can cause problems. Some of these can be avoided by defining _POSIX_C_SOURCE before including any
 *system headers, but it is better to only use all-upper-case names which have a unique prefix such as the package name.
 */

#ifdef __cplusplus
extern "C"
{
#endif
#include "bootstrap.h"
#include "exact.h"
#include "hitAndRun.h"
#include "main.h"
#include "utils_matrix.h"
#ifdef __cplusplus
}
#endif
#include <Rcpp.h>

/**
 * @brief Runs the Expected Maximization algorithm for every method.
 *
 * Given the stopping parameters of the EM method, it calculates an approximation of the RxG ecological inference
 * probability matrix.
 *
 * @param[in] Rcpp::String em_method The method for the EM algorithm. Options: "mvn_pdf", "mult", "exact".
 * (default: "mult")
 *
 * @param[in] Rcpp::String probability_method The method for obtaining the first probability. Options: "Group
 * proportional", "proportional", "uniform". (default: "Group proportional")
 *
 * @param[in] Rcpp::IntegerVector maximum_iterations A single integer value with the maximum iterations allowed for the
 * EM-algorithm. (default: 1000)
 *
 * @param[in] Rcpp:: maximum_seconds A single integer value with the maximum seconds to run the algorithm. (default:
 * 3600)
 *
 * @param[in] Rcpp::NumericVector stopping_threshold The absolute difference between subsequent probabilities matrices
 * to stop the algorithm. (default: 0.001)
 *
 * @param[in] Rcpp::LogicalVector verbose Boolean to determine if the algorithm will print helpful messages. (default:
 * false)
 *
 * @param[in] Rcpp::String monte_method The method to obtain an approximation of the CDF of the Normal vector.
 * The Alan Genz method are used, whereas it is heavily recommended to use the newest one, showing a faster and more
 * precise results. Options: "Genz", "Genz2". (default: "Genz2")
 *
 * @param[in] Rcpp::NumericVector monte_error The error threshold used to calculate the Montecarlo simulation
 * precition. The algorithm will do an early exit if this threshold is either accomplished or the maximum iterations are
 * done. (default: 0.000001)
 *
 * @param[in] Rcpp::IntegerVector monte_iter. The maximum amount of iterations to do in the Montecarlo
 * simulation. (default: 5000)
 *
 * @return Rcpp::List A list with the final probability ("result"), log-likelihood ("log_likelihood"), total
 * iterations that were made ("total_iterations"), time taken ("total_time"), stopping reason ("stopping_reason"),
 * finish id ("finish_id") and q value ("q").
 */
Rcpp::List EMAlgorithmFull(Rcpp::String em_method = "mult", Rcpp::String probability_method = "group_proportional",
                           Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
                           Rcpp::NumericVector maximum_seconds = Rcpp::NumericVector::create(3600),
                           Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
                           Rcpp::NumericVector log_threshold = Rcpp::NumericVector::create(-1000),
                           Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false),
                           Rcpp::IntegerVector step_size = Rcpp::IntegerVector::create(3000),
                           Rcpp::IntegerVector samples = Rcpp::IntegerVector::create(1000),
                           Rcpp::String monte_method = "genz2",
                           Rcpp::NumericVector monte_error = Rcpp::NumericVector(1e-6),
                           Rcpp::IntegerVector monte_iter = Rcpp::IntegerVector(5000));

/**
 * @brief Sets the `X` and `W` parameters on C
 *
 * Given an R's matrix, it sets the global parameters of `X` and `W` and computes all of its
 * important values (total candidates, votes per ballot, etc)
 *
 * @param Rcpp::NumericMatrix candidate_matrix A (c x b) matrix object of R that contains the votes that each
 * candidate `c` got on a ballot `b`.
 * @param Rcpp::NumericMatrix group_matrix A (b x g) matrix object of R that contains the votes that each
 * demographic group `g` did on a ballot `b`.
 *
 */
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix);

/**
 *  Returns an array of col-major matrices with bootstrapped matrices.
 *
 * @param[in] xmat The original X array
 *
 * @param[in] wmat The original W array
 *
 * @param[in] bootiter The amount of iterations for bootstrapping
 *
 * @param[in] p_method The method for obtaining the initial probability
 *
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 *
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
 *
 * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
 *
 * @param[in] maxSeconds Double with the value of the maximum amount of seconds to use.
 *
 * @param[in] verbose Wether to verbose useful outputs.
 *
 * @param[in, out] time The time that the algorithm took.
 *
 * @param[in, out] iterTotal Total amount of iterations.
 *
 * @param[in, out] logLLarr The loglikelihood array
 *
 * @param[in, out] finishing_reason The reason that the algorithm has been stopped. It can either be 0, 1, 2, 3,
 * representing a normal convergence, log likelihood decrease, maximum time reached and maximum iterations reached,
 * respectively.
 *
 *
 * @return An allocated array of size bootiter * TOTAL_BALLOTS that stores matrices.
 */

Rcpp::NumericMatrix bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                                 Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                                 Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_threshold,
                                 Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                                 Rcpp::IntegerVector samples, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);

/*
 * Returns a list with an heuristic-optimal bootstrapped matrix with an ideal group aggregation.
 *
 * @param[in] sd_statistic String indicates the statistic for the standard deviation (gxc) matrix. It can take the value
 * 'maximum', in which case computes the maximum over the standard deviation matrix, or 'average', in which case
 * computes the average.
 *
 * @param[in] String indicates the statistic for the standard deviation (gxc) matrix. It can take the value 'maximum',
 * in which case computes the maximum over the standard deviation matrix, or 'average', in which case computes the
 * average.
 *
 * @param[in] candidate_matrix The 'X' matrix of dimension (cxb).
 *
 * @param[in] group_matrix The 'W' matrix of dimension (bxg).
 *
 * @param[in] xmat The original X array
 *
 * @param[in] wmat The original W array
 *
 * @param[in] nboot The amount of iterations for bootstrapping
 *
 * @param[in] probability_method The method for obtaining the initial probability
 *
 * @param[in] em_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 *
 * @param[in] stopping_threshold Threshold value for convergence. Usually it's set to 0.001.
 *
 * @param[in] maximum_iterations Integer with a threshold of maximum iterations. Usually it's set to 100.
 *
 * @param[in] maximum_seconds Double with the value of the maximum amount of seconds to use.
 *
 * @param[in] verbose Wether to verbose useful outputs.
 *
 * @param[in] step_size The step size for the hnr method.
 *
 * @param[in] samples The amount of samples for the hnr method.
 *
 * @param[in] monte_method The method for the montecarlo simulation.
 *
 * @param[in] monte_error The error threshold for the montecarlo simulation.
 *
 * @param[in] monte_iter The amount of iterations for the montecarlo simulation
 *
 * @return A key-value list with 'bootstrap_result' and the cutting 'indices'.
 */
Rcpp::List groupAgg(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold, Rcpp::LogicalVector feasible,
                    Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Rcpp::IntegerVector nboot,
                    Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                    Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold,
                    Rcpp::NumericVector log_threshold, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                    Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                    Rcpp::IntegerVector monte_iter);

/*
 *
 * Greedy approach towards aggregating groups, evaluating with log-likelihoods
 *
 * @param[in] Rcpp::String em_method The method for the EM algorithm. Options: "mvn_pdf", "mult", "exact".
 * (default: "mult")
 *
 * @param[in] Rcpp::String probability_method The method for obtaining the first probability. Options: "Group
 * proportional", "proportional", "uniform". (default: "Group proportional")
 *
 * @param[in] Rcpp::IntegerVector maximum_iterations A single integer value with the maximum iterations allowed for the
 * EM-algorithm. (default: 1000)
 *
 * @param[in] Rcpp:: maximum_seconds A single integer value with the maximum seconds to run the algorithm. (default:
 * 3600)
 *
 * @param[in] Rcpp::NumericVector stopping_threshold The absolute difference between subsequent probabilities matrices
 * to stop the algorithm. (default: 0.001)
 *
 * @param[in] Rcpp::LogicalVector verbose Boolean to determine if the algorithm will print helpful messages. (default:
 * false)
 *
 * @param[in] Rcpp::String monte_method The method to obtain an approximation of the CDF of the Normal vector.
 * The Alan Genz method are used, whereas it is heavily recommended to use the newest one, showing a faster and more
 * precise results. Options: "Genz", "Genz2". (default: "Genz2")
 *
 * @param[in] Rcpp::NumericVector monte_error The error threshold used to calculate the Montecarlo simulation
 * precition. The algorithm will do an early exit if this threshold is either accomplished or the maximum iterations are
 * done. (default: 0.000001)
 *
 * @param[in] Rcpp::IntegerVector monte_iter. The maximum amount of iterations to do in the Montecarlo
 * simulation. (default: 5000)
 *
 */
Rcpp::List groupAggGreedy(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold,
                          Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                          Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                          Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                          Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                          Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples,
                          Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);

#endif // WRAPPER_H
