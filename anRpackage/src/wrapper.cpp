#include <Rcpp.h>
#include <iostream>

// Include the corrected wrapper.h
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
    std::cout << "Calling C function..." << std::endl;
    std::cout << "A matrix for X is" << std::endl;
    printMatrix(&XR);
    std::cout << "Now adding a column of zeroes" << std::endl;
    addColumnOfZeros(&XR, 0);
    printMatrix(&XR);
    freeMatrix(&XR);
    freeMatrix(&WR);

    // setParameters(&XR, &WR);
}
