#include <Rcpp.h>
#include <cstdio>
#include <iostream>
// Include the corrected wrapper.h
#include "Rcpp/String.h"
#include "wrapper.h"

static bool beenPrecomputed = false;
static bool matricesPrecomputed = false;

// [[Rcpp::export]]
Rcpp::List hitAndRunCDF(Rcpp::String Pmethod = "group proportional", Rcpp::IntegerVector maxIt = 1000,
                        Rcpp::NumericVector epsilon = 0.001, Rcpp::LogicalVector verbose = false,
                        Rcpp::String monteCarloMethod = "Genz2", Rcpp::NumericVector mvnEpsilon = 0.000001,
                        Rcpp::IntegerVector mvnIter = 10000)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = Pmethod;
    Matrix pIn = getInitialP(Pmethod.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    std::string monteMethod = monteCarloMethod;
    QMethodInput inputParams = {
        .monteCarloIter = mvnIter[0], .errorThreshold = mvnEpsilon[0], .simulationMethod = monteMethod.c_str()};

    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(maxIt[0] * sizeof(double));
    // ---...--- //

    Matrix Pnew =
        EMAlgoritm(&pIn, "MVN CDF", epsilon[0], maxIt[0], verbose, &timeIter, &totalIter, logLLarr, inputParams);
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
Rcpp::List hitAndRunEM(Rcpp::String Pmethod = "group proportional", Rcpp::IntegerVector maxIt = 1000,
                       Rcpp::NumericVector epsilon = 0.001, Rcpp::LogicalVector verbose = false,
                       Rcpp::IntegerVector stepSize = 3000, Rcpp::IntegerVector samples = 1000)
{
    // ---- Get the initial probability ---- //
    std::string probabilityM = Pmethod;
    Matrix pIn = getInitialP(Pmethod.get_cstring());
    // ---...--- //

    // ---- Define the main parameters ---- //
    QMethodInput inputParams = {.S = samples[0], .M = stepSize[0]};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(maxIt[0] * sizeof(double));
    // ---...--- //

    Matrix Pnew =
        EMAlgoritm(&pIn, "Hit and Run", epsilon[0], maxIt[0], verbose, &timeIter, &totalIter, logLLarr, inputParams);
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
void RprecomputeHR(Rcpp::IntegerVector S, Rcpp::IntegerVector M)
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
    generateOmegaSet(M[0], S[0], 42);
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
void RsetParameters(Rcpp::NumericMatrix x, Rcpp::NumericMatrix w, Rcpp::String jsonFile = "")
{

    // Case when the matrices were already computed and want to set other parameters
    if (matricesPrecomputed == true)
    {
        cleanup();
        matricesPrecomputed = false;
    }

    // Case when it wants to read the matrices from a JSON file.
    if (jsonFile != "")
    {
        std::string file = jsonFile;
        Matrix X, W, P;
        readJSONAndStoreMatrices(file.c_str(), &W, &X, &P);
        setParameters(&X, &W);
        matricesPrecomputed = true;
        return;
    }

    // Check dimensions
    if (x.nrow() == 0 || x.ncol() == 0)
    {
        std::cerr << "Error: X matrix has zero dimensions!" << std::endl;
        return;
    }
    if (w.nrow() == 0 || w.ncol() == 0)
    {
        std::cerr << "Error: W matrix has zero dimensions!" << std::endl;
        return;
    }

    // Convert to Matrix struct
    int xrows = x.nrow(), xcols = x.ncol();
    double *ptrX = (double *)malloc(xrows * xcols * sizeof(double));
    std::memcpy(ptrX, x.begin(), xrows * xcols * sizeof(double));
    Matrix XR = {ptrX, xrows, xcols};

    int wrows = w.nrow(), wcols = w.ncol();
    double *ptrW = (double *)malloc(wrows * wcols * sizeof(double));
    std::memcpy(ptrW, w.begin(), wrows * wcols * sizeof(double));

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
void readFilePrint(Rcpp::String filename, Rcpp::String method)
{
    std::string file = filename;
    std::string itMethod = method;

    Matrix X, W, P;

    readJSONAndStoreMatrices(file.c_str(), &W, &X, &P);
    setParameters(&X, &W);

    Matrix pIn = getInitialP("group proportional");

    QMethodInput inputParams = {.monteCarloIter = 100000, .errorThreshold = 0.000001, .simulationMethod = "Genz2"};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(10000 * sizeof(double));

    Matrix Pnew = EMAlgoritm(&pIn, itMethod.c_str(), 0.001, 1000, false, &timeIter, &totalIter, logLLarr, inputParams);
    free(logLLarr);
    printf("\nThe calculated matrix was\n");
    printMatrix(&Pnew);
    printf("\nThe real one was:\n");
    printMatrix(&P);
    printf("\nIt took %.5f seconds to run.", timeIter);
    freeMatrix(&P);
    freeMatrix(&Pnew);
    cleanup();
    // ---- Clean precomputed variables for a given iteration ---- //
    if (strcmp(itMethod.c_str(), "Exact") == 0)
    {
        cleanExact();
    }
    // ---- Hit and Run method ----
    else if (strcmp(itMethod.c_str(), "Hit and Run") == 0)
    {
        cleanHitAndRun();
    }
}
