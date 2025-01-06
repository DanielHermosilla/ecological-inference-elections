#ifndef COMPUTE_MULTINOMIAL_H
#define COMPUTE_MULTINOMIAL_H

#include "globals.h"
#include "matrixUtils.h"

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
double *computeQMultinomial(Matrix const *probabilities);

#endif
