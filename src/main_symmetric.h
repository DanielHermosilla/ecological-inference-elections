#ifndef MAIN_SYMMETRIC_H_EIM
#define MAIN_SYMMETRIC_H_EIM

#ifdef __cplusplus
extern "C"
{
#endif

#include "main.h"

    bool shouldRunSymmetricEMWeight(const QMethodInput *inputParams);

    void runSymmetricEMWeight(EMContext *ctx_forward, const char *p_method, const char *q_method,
                              const double convergence, const double LLconvergence, const int maxIter,
                              const double maxSeconds, const bool verbose, double *time, int *iterTotal,
                              double *logLLarr, int *finishing_reason, Matrix *probMatrix,
                              QMethodInput *inputParams, QMethodConfig config_forward);

#ifdef __cplusplus
}
#endif

#endif // MAIN_SYMMETRIC_H_EIM
