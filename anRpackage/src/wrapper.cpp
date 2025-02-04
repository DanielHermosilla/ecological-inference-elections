#include <Rcpp.h>
#include <cstdio>
#include <iostream>
// Include the corrected wrapper.h
#include "Rcpp/String.h"
#include "wrapper.h"

static bool beenPrecomputed = false;
static bool matricesPrecomputed = false;

// [[Rcpp::export]]
void hitAndRunEM(Rcpp::NumericMatrix X, Rcpp::NumericMatrix W, Rcpp::String Pmethod = "group proportional",
                 Rcpp::String jsonPath = "", Rcpp::IntegerVector maxIt = 1000, Rcpp::NumericVector epsilon = 0.001,
                 Rcpp::LogicalVector verbose = false, Rcpp::LogicalVector precompute = false,
                 Rcpp::IntegerVector stepSize = 1000)
{
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

// [[Rcpp::plugins(openmp)]]
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
