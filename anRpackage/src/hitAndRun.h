#ifndef HITANDRUN_H
#define HITANDRUN_H

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_combinations.h"
#include "utils_matrix.h"

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
     * Given a probability matrix with, it returns a flattened array with estimations of the conditional probability.
     * The array can be accesed with the macro `Q_3D` (it's a flattened tensor).
     *
     * @param[in] *probabilities. A pointer towards the probabilities matrix.
     * @param[in] params A QMethodInput struct with the `M` (step size) and `S` (samples) parameters
     *
     * @return A pointer towards the flattened tensor.
     *
     */
    double *computeQHitAndRun(Matrix const *probabilities, QMethodInput params);

    void cleanHitAndRun();

#ifdef __cplusplus
}
#endif
#endif // HITANDRUNH
