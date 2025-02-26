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
#include "Rcpp/String.h"
#include "main.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h>
#include <Rcpp.h>
#include <cstdio>
#include <iostream>
#include <vector>

// [[Rcpp::export]]
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix)
{

    // Check dimensions
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
    {
        Rcpp::stop("Error: X matrix has zero dimensions!");
    }
    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
    {
        Rcpp::stop("Error: W matrix has zero dimensions!");
    }

    // Convert to Matrix struct
    int xrows = candidate_matrix.nrow(), xcols = candidate_matrix.ncol();
    double *ptrX = (double *)R_alloc(xrows * xcols, sizeof(double));
    std::memcpy(ptrX, candidate_matrix.begin(), xrows * xcols * sizeof(double));
    Matrix XR = {ptrX, xrows, xcols};

    int wrows = group_matrix.nrow(), wcols = group_matrix.ncol();
    double *ptrW = (double *)R_alloc(wrows * wcols, sizeof(double));
    std::memcpy(ptrW, group_matrix.begin(), wrows * wcols * sizeof(double));

    Matrix WR = {ptrW, wrows, wcols};

    // Check pointers
    if (!XR.data || !WR.data)
    {
        Rcpp::stop("Error: Couldn't allocate memory for matrices. Consider removing unused objects or calling the "
                   "garbage collector (`gc()`)");
    }

    // Call C function
    setParameters(&XR, &WR);
    return;
}

// [[Rcpp::export]]
Rcpp::List EMAlgorithmFull(Rcpp::String em_method, Rcpp::String probability_method,
                           Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                           Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose,
                           Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                           Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    double timeIter = 0;
    int totalIter = 0;
    int finish = 0;
    std::string EMAlg = em_method;
    double logLLarr;
    double *qvalue = NULL;
    QMethodInput inputParams;
    if (EMAlg == "hnr")
    {
        QMethodInput inputParams = {.S = samples[0], .M = step_size[0]};
    }
    else if (EMAlg == "mvn_cdf")
    {
        static std::string monteMethod = monte_method;
        QMethodInput inputParams = {
            .monteCarloIter = monte_iter[0], .errorThreshold = monte_error[0], .simulationMethod = monteMethod.c_str()};
    }
    else
    {
        QMethodInput inputParams = {0};
    }
    // ---...--- //

    Matrix Pnew = EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], maximum_iterations[0], maximum_seconds[0],
                             verbose[0], &timeIter, &totalIter, &logLLarr, &qvalue, &finish, inputParams);
    freeMatrix(&pIn);
    std::string stopping_reason;

    switch (finish)
    {
    case 0:
        stopping_reason = "Convergence achieved";
        break;
    case 1:
        stopping_reason = "Log-likelihood decrease";
        break;
    case 2:
        stopping_reason = "Maximum time reached";
        break;
    case 3:
        stopping_reason = "Maximum iterations reached";
        break;
    default:
        stopping_reason = "Unknown";
        break;
    }

    if (verbose[0])
    {
        Rprintf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        Rprintf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr);
    }

    // ---- Return the results ---- //
    totalIter++;
    // ---- Final probability
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);
    // ---- q array
    // TODO: Maybe the error goes with TOTAL_BALLOTS, etc..
    Rcpp::NumericVector qTensor(qvalue, qvalue + TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS);
    free(qvalue);
    qvalue = NULL;
    // ---- Memory deallocation
    cleanup();
    if (EMAlg == "hnr")
        cleanHitAndRun();
    else if (EMAlg == "exact")
        cleanExact();

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = logLLarr,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = qTensor);
}

/*
Rcpp::List bootstrapAlg(Rcpp::IntegerVector nboot, Rcpp::String probability_method, Rcpp::IntegerVector
maximum_iterations, Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector
verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector
monte_error, Rcpp::IntegerVector monte_iter) {
*/
