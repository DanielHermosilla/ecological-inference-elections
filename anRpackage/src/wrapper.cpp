#include <Rcpp.h>
#include <cstdio>
#include <iostream>
// Include the corrected wrapper.h
#include "Rcpp/String.h"
#include "wrapper.h"

static bool beenPrecomputed = false;
static bool matricesPrecomputed = false;

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
    Matrix Pnew = EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], maximum_iterations[0], verbose, &timeIter,
                             &totalIter, logLLarr, inputParams);
    freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Clean the variables ---- //
    cleanup();
    matricesPrecomputed = false;
    if (EMAlg == "Exact")
    {
        cleanExact();
        beenPrecomputed = false;
    }
    // ---...--- //

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

    Matrix Pnew = EMAlgoritm(&pIn, "MVN CDF", stopping_threshold[0], maximum_iterations[0], verbose, &timeIter,
                             &totalIter, logLLarr, inputParams);
    freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Clean the variables ---- //
    cleanup();
    matricesPrecomputed = false;
    // ---...--- //

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

    Matrix Pnew = EMAlgoritm(&pIn, "Hit and Run", stopping_threshold[0], maximum_iterations[0], verbose, &timeIter,
                             &totalIter, logLLarr, inputParams);
    freeMatrix(&pIn);

    if (verbose)
    {
        printf("\nThe calculated matrix is\n");
        printMatrix(&Pnew);
        printf("\nIt took %.5f seconds to run with a log-likelihood of %.5f.\n", timeIter, logLLarr[totalIter]);
    }

    // ---- Clean the variables ---- //
    cleanup();
    matricesPrecomputed = false;
    cleanHitAndRun();
    beenPrecomputed = false;
    // ---...--- //

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
    // TODO: Enable a saving function
    if (beenPrecomputed)
    {
        std::cout << "There has been already a precomputation for Hit and Run" << std::endl;
        return;
    }
    if (!matricesPrecomputed)
    {
        std::cout << "The `X` and `W` matrices haven't been computed yet. Hint: Call RsetParameters." << std::endl;
        return;
    }

    // ---- Precompute the Omega Set and leave it on the global variable ---- //
    generateOmegaSet(step_size[0], samples[0], 42);
    // ---- Precompute a multinomial multiplication that is constant throughout the loops ---- //
    preComputeMultinomial();
    beenPrecomputed = true;
    return;
}

// [[Rcpp::export]]
void RprecomputeExact()
{
    // TODO: Enable a saving function
    if (beenPrecomputed)
    {
        std::cout << "There has been already a precomputation for the Exact method" << std::endl;
        return;
    }
    if (!matricesPrecomputed)
    {
        std::cout << "The `X` and `W` matrices haven't been computed yet. Hint: Call RsetParameters." << std::endl;
        return;
    }

    // ---- Precompute the H set and the K set and leave it as a global variable ---- //
    generateHSets();
    generateKSets();
    beenPrecomputed = true;
    return;
}

// [[Rcpp::export]]
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix)
{

    // Case when the matrices were already computed and want to set other parameters
    if (matricesPrecomputed == true)
    {
        cleanup();
        matricesPrecomputed = false;
    }
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
    matricesPrecomputed = true;
    return;
}

// [[Rcpp::export]]
void readFromFile(Rcpp::String filename)
{
    std::string file = filename;

    // Case when the matrices were already computed and want to set other parameters
    if (matricesPrecomputed == true)
    {
        cleanup();
        matricesPrecomputed = false;
    }

    // Case when it wants to read the matrices from a JSON file.
    Matrix RX, RW, RP;
    readJSONAndStoreMatrices(file.c_str(), &RW, &RX, &RP);
    setParameters(&RX, &RW);
    matricesPrecomputed = true;
    return;
}
