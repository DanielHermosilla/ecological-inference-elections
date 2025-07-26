#ifndef HITANDRUN_H_EIM
#define HITANDRUN_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "uthash.h"
#include "utils_hash.h"
#include "utils_matrix.h"

// ---...--- //

// ---- Macro for finding the minimum ---- //
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
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
    void computeQHitAndRun(EMContext *ctx, QMethodInput params, double *ll);

    /*
     * @brief Precomputes the sets used for the simulation.
     *
     * Precomputes the sets that are independent from each EM iteration. It is made with parallelism towards the
     * ballot boxes and with a static assignment for ensuring reproducibility.
     *
     * @param[in] M. The step size between consecutive samples. Note that the direction is assigned randomly.
     * @param[in] S. The amount of samples for each ballot box.
     *
     * @return void. Written on the global variable.
     */
    void generateOmegaSet(EMContext *ctx, int M, int S);

    double preMultinomialCoeff(EMContext *ctx, const int b, IntMatrix *currentMatrix);

    void preComputeMultinomial(EMContext *ctx);

    void encode(EMContext *ctx);

    void precomputeQConstant(EMContext *ctx, int size);

#ifdef __cplusplus
}
#endif
#endif // HITANDRUNH
