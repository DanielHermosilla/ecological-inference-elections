#ifndef WRAPPER_H
#define WRAPPER_H

#ifdef __cplusplus
extern "C"
{
#endif

#include "matrixUtils.h"

#ifdef __cplusplus
}
#endif
void RsetParameters(Rcpp::NumericMatrix x, Rcpp::NumericMatrix w);

#endif // WRAPPER_H
