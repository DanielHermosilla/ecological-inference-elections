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
Rcpp::List EMAlgorithmAll(Rcpp::String em_method, Rcpp::String probability_method,
                          Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                          Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    QMethodInput inputParams = {0};
    double timeIter = 0;
    int totalIter = 0;
    int finish = 0;
    double *logLLarr = (double *)R_alloc(maximum_iterations[0], sizeof(double));
    // ---...--- //

    std::string EMAlg = em_method;
    Matrix Pnew = EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], maximum_iterations[0], maximum_seconds[0],
                             verbose[0], &timeIter, &totalIter, logLLarr, &finish, inputParams);
    freeMatrix(&pIn);
    std::string stopping_reason;

    switch (finish)
    {
    case 0:
        stopping_reason = "Convergence achieved";
        break;
    case 1:
        stopping_reason = "Maximum time reached";
        break;
    case 2:
        stopping_reason = "Maximum iterations reached";
        break;
        /*
    case 1:
        stopping_reason = "Log-likelihood decrease";
        break;
        */
    default:
        stopping_reason = "Unknown";
        break;
    }

    if (verbose[0])
    {
        Rprintf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        Rprintf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
        Rprintf("\nThe reason why the algorithm ended was:\t'%s'\n", stopping_reason.c_str());
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    totalIter++;
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish);
    // ---...--- //
}

// [[Rcpp::export]]
Rcpp::List EMAlgorithmCDF(Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                          Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold,
                          Rcpp::LogicalVector verbose, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                          Rcpp::IntegerVector monte_iter)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    std::string monteMethod = monte_method;
    QMethodInput inputParams = {
        .monteCarloIter = monte_iter[0], .errorThreshold = monte_error[0], .simulationMethod = monteMethod.c_str()};

    double timeIter = 0;
    int totalIter = 0;
    int finish = 0;
    double *logLLarr = (double *)R_alloc(maximum_iterations[0], sizeof(double));
    // ---...--- //

    Matrix Pnew = EMAlgoritm(&pIn, "mvn_cdf", stopping_threshold[0], maximum_iterations[0], maximum_seconds[0],
                             verbose[0], &timeIter, &totalIter, logLLarr, &finish, inputParams);
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
        Rprintf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    totalIter++;
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish);
    // ---...--- //
}

// [[Rcpp::export]]
Rcpp::List EMAlgorithmHitAndRun(Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                                Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold,
                                Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    QMethodInput inputParams = {.S = samples[0], .M = step_size[0]};
    double timeIter = 0;
    int totalIter = 0;
    int finish = 0;
    double *logLLarr = (double *)R_alloc(maximum_iterations[0], sizeof(double));
    // ---...--- //

    Matrix Pnew = EMAlgoritm(&pIn, "hnr", stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0],
                             &timeIter, &totalIter, logLLarr, &finish, inputParams);
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
        Rprintf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    totalIter++;
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish);
    // ---...--- //
}

// [[Rcpp::export]]
void RprecomputeHR(Rcpp::IntegerVector samples, Rcpp::IntegerVector step_size)
{
    // ---- Precompute the Omega Set and leave it on the global variable ---- //
    generateOmegaSet(step_size[0], samples[0], 42); // TODO: Rely on R's seed and randomizer
    // ---- Precompute a multinomial multiplication that is constant throughout the loops ---- //
    preComputeMultinomial();
    return;
}

// [[Rcpp::export]]
void RprecomputeExact()
{
    // ---- Precompute the H set and the K set and leave it as a global variable ---- //
    generateHSets();
    generateKSets();
    return;
}

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
void clean_exact_precompute()
{
    cleanExact();
}

// [[Rcpp::export]]
void clean_hr_precompute()
{
    cleanHitAndRun();
}

// [[Rcpp::export]]
void clean_everything()
{
    cleanup();
}
