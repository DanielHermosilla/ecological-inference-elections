#include "globals.h"
#include "main.h"
#include "utils_matrix.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <dirent.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/*
Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                  const double maxSeconds, const bool verbose, double *time, int *iterTotal, double *logLLarr,
                  double *qVal, int *finishing_reason, QMethodInput inputParams
*/
double *bootstrap(int bootiter, Matrix *currentP, const char *q_method, const double convergence, const int maxIter,
                  const double maxSeconds, const bool verbose, QMethodInput inputParams)
{
}
