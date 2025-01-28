#include "main.h"
#include "matrixUtils.h"
#include <R.h>
#include <Rinternals.h>

SEXP R_EM_algorithm(SEXP X, SEXP W, SEXP maxIterations, SEXP convergence, SEXP method, SEXP verbose, SEXP pMethod)
{

    // ---- Create the matrices on the C standard ---- //
    // ---- Create the X matrix ----
    int xcols = ncols(X); // Ballot boxes
    int xrows = nrows(X); // Candidates
    double *xData = REAL(X);
    Matrix RX = {.cols = xcols, .rows = xrows, .data = xData};

    // ---- Create the W matrix ----
    int wcols = ncols(W); // Demographic groups
    int wrows = nrows(W); // Ballot boxes
    double *wData = REAL(W);
    Matrix RW = {.cols = wcols, .rows = wrows, .data = wData};
    // ---...--- //

    // ---- Set the global parameters ---- //
    setParameters(&RX, &RW);

    // ---...--- //

    // ---- Calculate the initial probability ---- //
    const char *probabilityMethod = CHAR(pMethod);
    Matrix initialProbability = getInitialP(probabilityMethod);
    // ---...--- //

    // ---- Get the EM stopping parameters ---- //
    int maximumIterations = INTEGER(maxIterations)[0];
    double convergenceThreshold = REAL(convergence)[0];
}
