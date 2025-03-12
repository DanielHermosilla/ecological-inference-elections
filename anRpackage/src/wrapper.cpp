#include <Rcpp.h>
#include <cstdio>
#include <iostream>
// Include the corrected wrapper.h
#include "Rcpp/String.h"
#include "main.h"
#include "wrapper.h"

// [[Rcpp::export]]
Rcpp::List EMAlgorithmAll(Rcpp::String em_method, Rcpp::String probability_method,
                          Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector stopping_threshold,
                          Rcpp::LogicalVector verbose)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    QMethodInput inputParams = {0};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(maximum_iterations[0] * sizeof(double));
    // ---...--- //

    std::string EMAlg = em_method;
    Matrix Pnew = EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], maximum_iterations[0], verbose[0], &timeIter,
                             &totalIter, logLLarr, inputParams);
    freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);
    free(logLLarr);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter);
    // ---...--- //
}

// [[Rcpp::export]]
Rcpp::List EMAlgorithmCDF(Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                          Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose,
                          Rcpp::String multivariate_method, Rcpp::NumericVector multivariate_epsilon,
                          Rcpp::IntegerVector multivariate_iterations)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    std::string monteMethod = multivariate_method;
    QMethodInput inputParams = {.monteCarloIter = multivariate_iterations[0],
                                .errorThreshold = multivariate_epsilon[0],
                                .simulationMethod = monteMethod.c_str()};

    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(maximum_iterations[0] * sizeof(double));
    // ---...--- //

    Matrix Pnew = EMAlgoritm(&pIn, "MVN CDF", stopping_threshold[0], maximum_iterations[0], verbose[0], &timeIter,
                             &totalIter, logLLarr, inputParams);
    freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);
    free(logLLarr);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter);
    // ---...--- //
}

// [[Rcpp::export]]
Rcpp::List EMAlgorithmHitAndRun(Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                                Rcpp::NumericVector stopping_threshold, Rcpp::LogicalVector verbose,
                                Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = probability_method;
    Matrix pIn = getInitialP(probability_method.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    QMethodInput inputParams = {.S = samples[0], .M = step_size[0]};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(maximum_iterations[0] * sizeof(double));
    // ---...--- //

    Matrix Pnew = EMAlgoritm(&pIn, "Hit and Run", stopping_threshold[0], maximum_iterations[0], verbose[0], &timeIter,
                             &totalIter, logLLarr, inputParams);
    // freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Return the results ---- //
    // ---- Final probability ----
    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    // ---- Final log-likelihood array ----
    Rcpp::NumericVector RlogLikelihood(logLLarr, logLLarr + totalIter);
    free(logLLarr);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = RlogLikelihood,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter);
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
        std::cerr << "Error: X matrix has zero dimensions!" << std::endl;
        return;
    }
    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
    {
        std::cerr << "Error: W matrix has zero dimensions!" << std::endl;
        return;
    }

    // Convert to Matrix struct
    int xrows = candidate_matrix.nrow(), xcols = candidate_matrix.ncol();
    double *ptrX = (double *)malloc(xrows * xcols * sizeof(double));
    std::memcpy(ptrX, candidate_matrix.begin(), xrows * xcols * sizeof(double));
    Matrix XR = {ptrX, xrows, xcols};

    int wrows = group_matrix.nrow(), wcols = group_matrix.ncol();
    double *ptrW = (double *)malloc(wrows * wcols * sizeof(double));
    std::memcpy(ptrW, group_matrix.begin(), wrows * wcols * sizeof(double));

    Matrix WR = {ptrW, wrows, wcols};

    // Check pointers
    if (!XR.data || !WR.data)
    {
        std::cerr << "Error: Matrix data pointer is NULL!" << std::endl;
        return;
    }

    // Call C function
    setParameters(&XR, &WR);
    return;
}

/*
void readFromFile(Rcpp::String filename)
{
    std::string file = filename;


    // Case when it wants to read the matrices from a JSON file.
    Matrix RX, RW, RP;
    readJSONAndStoreMatrices(file.c_str(), &RW, &RX, &RP);
    setParameters(&RX, &RW);
    return;
}
*/

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
