/*
Copyright (c) 2025 Daniel Hermosilla

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "wrapper.h"
#include "Rcpp/vector/instantiation.h"
#include "bootstrap.h"
#include "main.h"
#include <Rcpp.h>
#include <vector>

Matrix convertToMatrix(const Rcpp::NumericMatrix &mat)
{
    int rows = mat.nrow(), cols = mat.ncol();
    double *data = (double *)malloc(rows * cols * sizeof(double)); // Allocate on heap
    std::memcpy(data, mat.begin(), rows * cols * sizeof(double));  // Copy data from R matrix
    return {data, rows, cols};                                     // Safe to return
}

// ---- Helper Function: Initialize QMethodInput ---- //
QMethodInput initializeQMethodInput(const std::string &EMAlg, int samples, int step_size, int monte_iter,
                                    double monte_error, const std::string &monte_method)
{
    QMethodInput inputParams = {0}; // Default initialization

    if (EMAlg == "hnr")
    {
        inputParams.S = samples;
        inputParams.M = step_size;
    }
    else if (EMAlg == "mvn_cdf")
    {
        inputParams.monteCarloIter = monte_iter;
        inputParams.errorThreshold = monte_error;
        inputParams.simulationMethod = monte_method.c_str();
    }

    return inputParams;
}

// ---- Set Parameters ---- //
// [[Rcpp::export]]
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    setParameters(&XR, &WR);
}

// ---- Run EM Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmFull(Rcpp::String em_method, Rcpp::String probability_method,
                           Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                           Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose,
                           Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                           Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;

    Matrix pIn = getInitialP(probabilityM.c_str());

    double timeIter = 0, logLLarr = 0;
    int totalIter = 0, finish = 0;
    double *qvalue = NULL;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    Matrix Pnew = EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], maximum_iterations[0], maximum_seconds[0],
                             verbose[0], &timeIter, &totalIter, &logLLarr, &qvalue, &finish, inputParams);

    // ---- Create human-readable stopping reason ---- //
    std::vector<std::string> stop_reasons = {"Convergence achieved", "Log-likelihood decrease", "Maximum time reached",
                                             "Maximum iterations reached"};
    std::string stopping_reason = (finish >= 0 && finish < 4) ? stop_reasons[finish] : "Unknown";

    if (verbose[0])
    {
        Rprintf("\nThe calculated matrix is:\n");
        printMatrix(&Pnew);
        Rprintf("\nExecution Time: %.5f seconds | Log-likelihood: %.5f.\n", timeIter, logLLarr);
    }

    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    Rcpp::NumericVector condProb(qvalue, qvalue + TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS);
    free(qvalue);
    cleanup();
    if (EMAlg == "hnr")
        cleanHitAndRun();
    else if (EMAlg == "exact")
        cleanExact();

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = logLLarr,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb);
}

// ---- Run Bootstrapping Algorithm ---- //
// [[Rcpp::export]]
Rcpp::NumericVector bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                                 Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                                 Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose,
                                 Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    double *result = bootstrapA(&XR, &WR, nboot[0], EMAlg.c_str(), probabilityM.c_str(), stopping_threshold[0],
                                maximum_iterations[0], maximum_seconds[0], verbose[0], inputParams);

    // double *bootstrapA(const Matrix *xmat, const Matrix *wmat, int bootiter, const char *q_method, const char
    // *p_method, const double convergence, const int maxIter, const double maxSeconds, const bool verbose, QMethodInput
    // inputParams);
    int result_size = nboot[0] * XR.rows * WR.cols;
    Rcpp::NumericVector output(result, result + result_size);
    free(result);

    return output;
}
