// ---- Avoid circular dependencies
#ifndef GLOBALS_H
#define GLOBALS_H

#ifdef __cplusplus

extern "C"
{
#endif

#include "utils_matrix.h"
#include <stdint.h>

// Macro for accessing a 3D flattened array (b x g x c)
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])

    extern uint32_t TOTAL_VOTES;
    extern uint32_t TOTAL_BALLOTS;
    extern uint16_t TOTAL_CANDIDATES;
    extern uint16_t TOTAL_GROUPS;
    extern uint16_t *BALLOTS_VOTES;    // Total votes per ballot
    extern uint32_t *CANDIDATES_VOTES; // Total votes per candidate
    extern uint32_t *GROUP_VOTES;      // Total votes per group
    extern double *inv_BALLOTS_VOTES;
    extern Matrix *X;
    extern Matrix *W;
#ifdef __cplusplus
}
#endif
#endif // GLOBALS_H
