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

#endif // WRAPPER_H
