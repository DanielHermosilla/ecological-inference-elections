#include <Rcpp.h>
#include <cstdio>
#include <iostream>
// Include the corrected wrapper.h
#include "Rcpp/String.h"
#include "wrapper.h"

// [[Rcpp::export]]
void RsetParameters(Rcpp::NumericMatrix x, Rcpp::NumericMatrix w)
{
    std::cout << "Function RsetParameters called" << std::endl;

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

    std::cout << "X Matrix: " << xrows << "x" << xcols << std::endl;
    std::cout << "W Matrix: " << wrows << "x" << wcols << std::endl;

    // Check pointers
    if (!XR.data || !WR.data)
    {
        std::cerr << "Error: Matrix data pointer is NULL!" << std::endl;
        return;
    }
    // Call C function
    // std::cout << "This is a trial for setting parameters:\t" << std::endl;
    setParameters(&XR, &WR);
    // freeMatrix(&XR);
    // freeMatrix(&WR);
    printf("The X matrix is:\n");
    printMatrix(&XR);
    printf("The W matrix is:\n");
    printMatrix(&WR);
    printf("\nGroup proportional method:\n");
    Matrix p = getInitialP("group proportional");
    printMatrix(&p);
    printf("\nProportional method:\n");
    p = getInitialP("proportional");
    printMatrix(&p);
    printf("\nUniform method:\n");
    p = getInitialP("uniform");
    printMatrix(&p);

    // setParameters(&XR, &WR);
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
void readFilePrint(Rcpp::String filename, Rcpp::String method)
{
    std::string file = filename;
    std::string itMethod = method;

    Matrix X, W, P;

    readJSONAndStoreMatrices(file.c_str(), &W, &X, &P);
    /*
    printf("The matrices are:\nFor W:\n");
    printMatrix(&W);
    printf("\nFor X:\n");
    printMatrix(&X);
    printf("\nFor P:\n");
    printMatrix(&P);
*/
    setParameters(&X, &W);
    /*
    printf("\nGroup proportional method:\n");
    Matrix pIn = getInitialP("group proportional");
    printMatrix(&pIn);
    freeMatrix(&pIn);

    printf("\nProportional method:\n");
    pIn = getInitialP("proportional");
    printMatrix(&pIn);
    freeMatrix(&pIn);

    printf("\nUniform method:\n");
    pIn = getInitialP("uniform");
    printMatrix(&pIn);
    freeMatrix(&pIn);
*/
    Matrix pIn = getInitialP("group proportional");

    // It will run multinomial...
    QMethodInput inputParams = {.monteCarloIter = 100000, .errorThreshold = 0.00001, .simulationMethod = "Genz2"};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLarr = (double *)malloc(10000 * sizeof(double));
    // double *logLLresults = nullptr;

    Matrix Pnew = EMAlgoritm(&pIn, itMethod.c_str(), 0.001, 10000, true, &timeIter, &totalIter, logLLarr, inputParams);
    free(logLLarr);
    printf("\nThe calculated matrix was\n");
    printMatrix(&Pnew);
    printf("\nThe real one was:\n");
    printMatrix(&P);
    printf("\nIt took %.5f seconds to run.", timeIter);
    // free(&logLLresults);

    // printMatrix(&p);
    // freeMatrix(&W);
    // freeMatrix(&X);
    freeMatrix(&P);
    freeMatrix(&Pnew);
    // free(W.data);
    // free(X.data);
    // free(P.data);
    cleanup();
    // *logLLresults = NULL;
}
