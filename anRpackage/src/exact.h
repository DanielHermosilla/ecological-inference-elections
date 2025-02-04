#ifndef COMPUTE_EXACT_H
#define COMPUTE_EXACT_H

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_hash.h"
#include "utils_matrix.h"

    // ---- Define a structure to store the sets ---- //
    typedef struct
    {
        uint32_t b;
        uint16_t g;
        size_t **data;
        size_t size;
    } Set;
    // ---...--- //

    /**
     * @brief Calculate the value of `q_{bgc}`.
     *
     * It calculates all of the values for q_{bgc} by the definition on the paper. It returns the array of type
     * `double`.
     *
     * @param[in] *probabilities A pointer to the matrix with the probabilities.
     * @param[in] params The parameters to use for the `q` probability. On this case, it should be empty.
     *
     * @return *double: A pointer toward the array.
     *
     * @note: A single pointer is used to store the array continously. This is for using cBLAS operations later.
     *
     */
    double *computeQExact(const Matrix *probabilities, QMethodInput params);

    /**
     * @brief Cleans all of the allocated memory associated with the exact method
     *
     * Given the precomputed sets of possibilities, it frees everything.
     */
    void cleanExact();
    /**
     * @brief Precalculates the `H` set for every index.
     *
     * Given that the set wasn't calculated before, it calculates the H set defined as every possible combination for a
     * given `g` group
     *
     * @return void: Results written at the global variable HSETS.
     *
     */
    void generateHSets();

    /**
     * @brief Precalculates the `K` set for every index.
     *
     * Given that the set wasn't calculated before, it calculates the K set defined as every possible combination as a
     * cummulative set given the first `f` groups.
     *
     * @return void: Results written at the global variable KSETS.
     *
     */

    void generateKSets();

#ifdef __cplusplus
}
#endif
#endif
