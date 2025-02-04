#ifndef WRAPPER_H
#define WRAPPER_H

#ifdef __cplusplus
extern "C"
{
#endif
#include "exact.h"
#include "globals.h"
#include "hitAndRun.h"
#include "main.h"
#include "utils_file.h"
#include "utils_matrix.h"

#ifdef __cplusplus
}
#endif
void RsetParameters(Rcpp::NumericMatrix x, Rcpp::NumericMatrix w);

void readFilePrint(Rcpp::String filename);

#endif // WRAPPER_H
