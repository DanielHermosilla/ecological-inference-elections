#ifndef COMPUTE_EXACT_H
#define COMPUTE_EXACT_H

#include "matrixUtils.h"
#include <stdint.h>
#include <stdlib.h>

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

// Macro for accessing a 3D flattened array (b x g x c)
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])
/**
 * @brief Computes an approximate conditional probability using a Multinomial approach.
 *
 * @param[in] candidates       Matrix (b x c) - results of candidate "c" in ballot "b".
 * @param[in] groups           Matrix (b x g) - votes from demographic group "g".
 * @param[in] probabilities    Matrix (g x c) - probabilities of each group and candidate.
 * @param[in] candidatesVotes  Array (b) - total votes per ballot.
 * @param[in] candidatesTotal  Pointer to the total number of candidates.
 * @param[in] ballotsTotal     Pointer to the total number of ballots.
 * @param[in] groupsTotal      Pointer to the total number of groups.
 *
 * @return A pointer to a flattened 3D array (b x g x c) representing the probabilities.
 */
// double *computeQMultinomial(Matrix const *probabilities);

#endif
