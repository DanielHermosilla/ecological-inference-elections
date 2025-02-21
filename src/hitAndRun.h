#ifndef COMPUTE_HITANDRUN_H
#define COMPUTE_HITANDRUN_H

#include "combinations.h"
#include "globals.h"
#include "matrixUtils.h"

// ---- Define a structure to store the Omega sets ---- //
typedef struct
{
    uint32_t b;
    Matrix **data;
    size_t size;
} OmegaSet;
// ---...--- //

// ---- Macro for finding the minimum ---- //
#define MIN(a, b) ((a) < (b) ? (a) : (b))
// ---...--- //

/**
 * @brief Computes the `q` values for all the ballot boxes given a probability matrix. Uses the Hit and Run method.
 *
 * Given a probability matrix with, it returns a flattened array with estimations of the conditional probability. The
 * array can be accesed with the macro `Q_3D` (it's a flattened tensor).
 *
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 * @param[in] M. The step size used for the creation of samples.
 * @param[in] S. The total amount of samples
 *
 * @return A pointer towards the flattened tensor.
 *
 */
double *computeQHitAndRun(Matrix const *probabilities, int M, int S);

void cleanHitAndRun();
#endif
