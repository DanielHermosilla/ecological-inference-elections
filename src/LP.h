#ifndef LP_EIM
#define LP_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"

    // ---...--- //
    int LP_NW(const Matrix *X, const Matrix *W, Matrix *q_bgc, int b);
    int LPW(const Matrix *X, const Matrix *W, Matrix *q_bgc, int b);
    int LP_NW_ctx(EMContext *ctx, int b);
    int LPW_ctx(EMContext *ctx, int b);

#ifdef __cplusplus
}
#endif
#endif
