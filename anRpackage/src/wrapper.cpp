#include <Rcpp.h>
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

// [[Rcpp::export]]
void readFilePrint(Rcpp::String filename)
{
    std::string file = filename;

    Matrix X, W, P;

    readJSONAndStoreMatrices(file.c_str(), &W, &X, &P);
    printf("The matrices are:\nFor W:\n");
    printMatrix(&W);
    printf("\nFor X:\n");
    printMatrix(&X);
    printf("\nFor P:\n");
    printMatrix(&P);

    setParameters(&X, &W);
    printf("\nGroup proportional method:\n");
    printf("\nProportional method:\n");
    Matrix pIn = getInitialP("proportional");
    printMatrix(&pIn);
    printf("\nUniform method:\n");
    pIn = getInitialP("uniform");
    printMatrix(&pIn);
    pIn = getInitialP("group proportional");
    printMatrix(&pIn);

    // It will run multinomial...
    QMethodInput inputParams = {0};
    double timeIter = 0;
    int totalIter = 0;
    double *logLLresults = NULL;

    Matrix Pnew = EMAlgoritm(&pIn, "Multinomial", 0.001, 1000, true, &timeIter, &totalIter, &logLLresults, inputParams);
    printf("\nThe calculated matrix was\n");
    printMatrix(&Pnew);
    printf("\nThe real one was:\n");
    printMatrix(&P);
    // printMatrix(&p);
    freeMatrix(&W);
    freeMatrix(&X);
    freeMatrix(&P);
}
