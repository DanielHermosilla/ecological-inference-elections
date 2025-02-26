#ifndef HITANDRUN_H_EIM
#define HITANDRUN_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
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

    /*
     * @brief Precomputes the sets used for the simulation.
     *
     * Precomputes the sets that are independent from each EM iteration. It is made with parallelism towards the ballot
     * boxes and with a static assignment for ensuring reproducibility.
     *
     * @param[in] M. The step size between consecutive samples. Note that the direction is assigned randomly.
     * @param[in] S. The amount of samples for each ballot box.
     *
     * @return void. Written on the global variable.
     */
    void generateOmegaSet(int M, int S);

    void preComputeMultinomial();
    void cleanHitAndRun();

#ifdef __cplusplus
}
#endif
#endif // HITANDRUNH
