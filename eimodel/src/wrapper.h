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
#include "exact.h"
#include "globals.h"
#include "hitAndRun.h"
#include "main.h"
#include "utils_matrix.h"

#ifdef __cplusplus
}
#endif

/**
 * @brief Runs the Expected Maximization algorithm for the `MVN PDF`, `Multinomial` or `Exact` method.
 *
 * Given the stopping parameters of the EM method, it calculates an approximation of the RxG ecological inference
 * probability matrix.
 *
 * @param[in] Rcpp::String em_method The method for the EM algorithm. Options: "MVN PDF", "Multinomial", "Exact".
 * (default: "Multinomial")
 * @param[in] Rcpp::String probability_method The method for obtaining the first probability. Options: "Group
 * proportional", "Proportional", "Uniform". (default: "Group proportional")
 * @param[in] Rcpp::IntegerVector maximum_iterations A single integer value with the maximum iterations allowed for the
 * EM-algorithm. (default: 1000)
 * @param[in] Rcpp:: maximum_minutes A single integer value with the maximum minutes to run the algorithm. (default:
 * 1440)
 * @param[in] Rcpp::NumericVector stopping_threshold The absolute difference between subsequent probabilities matrices
 * to stop the algorithm. (default: 0.001)
 * @param[in] Rcpp::LogicalVector verbose Boolean to determine if the algorithm will print helpful messages. (default:
 * false)
 *
 * @note For running the `Hit and Run` and `MVN CDF` method refer to its own dedicated functions.
 *
 * @return Rcpp::List A list with the final probability ("result"), log-likelihood array ("log_likelihood"), total
 * iterations that were made ("total_iterations") and the time that was taken ("total_time").
 */
Rcpp::List EMAlgorithmAll(Rcpp::String em_method = "Multinomial",
                          Rcpp::String probability_method = "Group proportional",
                          Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
                          Rcpp::IntegerVector maximum_minutes = Rcpp::IntegerVector::create(1440),
                          Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
                          Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false));

/**
 * @brief Runs the Expected Maximization algorithm for the Multivariate CDF method.
 *
 * Given the stopping parameters of the EM method and the Multivariate method (as shown in the Alan Genz paper), it
 * calculates the Expected Maximization algorithm for the Multivariate CDF method. Note that this method depends heavily
 * on randomized sampling, hence, it is to expect some variance between iterations that could give innacurate log
 * likelihood values. Higher iterations and lower error thresholds will improve accuracy at the cost of computation
 * time.
 *
 * @param[in] Rcpp::String probability_method The method for obtaining the first probability. Options: "Group
 * proportional", "Proportional", "Uniform". (default: "Group proportional")
 * @param[in] Rcpp::IntegerVector maximum_iterations A single integer value with the maximum iterations allowed for the
 * EM-algorithm. (default: 1000)
 * @param[in] Rcpp:: maximum_minutes A single integer value with the maximum minutes to run the algorithm. (default:
 * 1440)
 * @param[in] Rcpp::NumericVector stopping_threshold The absolute difference between subsequent probabilities matrices
 * to stop the algorithm. (default: 0.001)
 * @param[in] Rcpp::LogicalVector verbose Boolean to determine if the algorithm will print helpful messages. (default:
 * false)
 * @param[in] Rcpp::String multivariate_method The method to obtain an approximation of the CDF of the Normal vector.
 * The Alan Genz method are used, whereas it is heavily recommended to use the newest one, showing a faster and more
 * precise results. Options: "Genz", "Genz2". (default: "Genz2")
 * @param[in] Rcpp::NumericVector multivariate_epsilon The error threshold used to calculate the Montecarlo simulation
 * precition. The algorithm will do an early exit if this threshold is either accomplished or the maximum iterations are
 * done. (default: 0.000001)
 * @param[in] Rcpp::IntegerVector multivariate_iterations. The maximum amount of iterations to do in the Montecarlo
 * simulation. (default: 5000)
 *
 * @note For running the other EM algorithm methods, refer to its own functions.
 *
 * @return Rcpp::List A list with the final probability ("result"), log-likelihood array ("log_likelihood"), total
 * iterations that were made ("total_iterations") and the time that was taken ("total_time").
 */
Rcpp::List EMAlgorithmCDF(Rcpp::String probability_method = "Group proportional",
                          Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
                          Rcpp::IntegerVector maximum_minutes = Rcpp::IntegerVector::create(1440),
                          Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
                          Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false),
                          Rcpp::String multivariate_method = "Genz2",
                          Rcpp::NumericVector multivariate_epsilon = Rcpp::NumericVector::create(0.000001),
                          Rcpp::IntegerVector multivariate_iterations = Rcpp::IntegerVector::create(5000));

/**
 * @brief Runs the Expected Maximization algorithm for the Hit and Run method.
 *
 * Given the stopping parameters of the EM method and the Hit and Run method, it
 * calculates the Expected Maximization algorithm. As shown in the paper, a higher step size would reduce the
 * correlation between samples, making a better sampling.
 *
 * @param[in] Rcpp::String probability_method The method for obtaining the first probability. Options: "Group
 * proportional", "Proportional", "Uniform". (default: "Group proportional")
 * @param[in] Rcpp::IntegerVector maximum_iterations A single integer value with the maximum iterations allowed for
 * the EM-algorithm. (default: 1000)
 * @param[in] Rcpp:: maximum_minutes A single integer value with the maximum minutes to run the algorithm. (default:
 * 1440)
 * @param[in] Rcpp::NumericVector stopping_threshold The absolute difference between subsequent probabilities
 * matrices to stop the algorithm. (default: 0.001)
 * @param[in] Rcpp::LogicalVector verbose Boolean to determine if the algorithm will print helpful messages.
 * (default: false)
 * @param[in] Rcpp::IntegerVector step_size The distance between sampling points, as shown in the paper. (default: 3000)
 * @param[in] Rcpp::IntegerVector samples Amount of samples. (default: 1000)
 *
 * @note Much of the time taken on this algorithm is due to the sampling generation. This can be precomputed and also
 * parallelized.
 *
 * @return Rcpp::List A list with the final probability ("result"), log-likelihood array ("log_likelihood"), total
 * iterations that were made ("total_iterations") and the time that was taken ("total_time").
 */
Rcpp::List EMAlgorithmHitAndRun(Rcpp::String probability_method = "Group proportional",
                                Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
                                Rcpp::IntegerVector maximum_minutes = Rcpp::IntegerVector::create(1440),
                                Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
                                Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false),
                                Rcpp::IntegerVector step_size = Rcpp::IntegerVector::create(3000),
                                Rcpp::IntegerVector samples = Rcpp::IntegerVector(1000));

/**
 * @brief Precomputes the values from the Hit and Run method that can be calculated before the EM algorithm.
 *
 * Computes all of the sets that are going to be used among all of the EM iterations. This values are independent
 * between each stage of the algorithm.
 *
 * @param[in] Rcpp::IntegerVector samples The amount of samples to use. (default: 1000)
 * @param[in] Rcpp::IntegerVector step_size The distance between sampling points, as shown in the paper. (default: 3000)
 *
 * @warning: The `X` and `W` matrices must be assigned before with the `setParameters` function.
 *
 * @note: This values must be freed with the `cleanHitAndRun()` function.
 *
 */
void RprecomputeHR(Rcpp::IntegerVector samples = Rcpp::IntegerVector::create(1000),
                   Rcpp::IntegerVector step_size = Rcpp::IntegerVector::create(3000));

/**
 * @brief Precomputes the values from the Exact method that can be calculated before the EM algorithm.
 *
 * Computes all of the sets that are going to be used among all of the EM iterations. This values are independent
 * between each stage of the algorithm.
 *
 * @warning: The `X` and `W` matrices must be assigned before with the `setParameters` function.
 *
 * @note: This values must be freed with the `cleanExact()` function.
 *
 */
void RprecomputeExact();

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
 * void readFromFile(Rcpp::String filename);
 */

/**
 * Cleans the precomputed values for the exact method
 *
 */
void clean_exact_precompute();

/**
 * Cleans the precomputed values for the Hit and Run method
 *
 */
void clean_hr_precompute();

/**
 * Cleans the global variables defined on C
 *
 */
void clean_everything();

#endif // WRAPPER_H
