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
    double *ptrX = x.begin();
    Matrix XR = {ptrX, xrows, xcols};

    int wrows = w.nrow(), wcols = w.ncol();
    double *ptrW = w.begin();
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
    // setParameters(&XR, &WR);
}
