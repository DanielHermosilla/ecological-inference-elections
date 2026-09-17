#ifndef KL_EIM
#define KL_EIM

#ifdef __cplusplus
extern "C"
{
#endif

#include "globals.h"

    int KL_joint_symmetric(const Matrix *X, const Matrix *W, Matrix *q_forward, Matrix *q_reverse, int b);
    int KL_joint_symmetric_ctx(EMContext *ctx_forward, EMContext *ctx_reverse, int b);

#ifdef __cplusplus
}
#endif

#endif
