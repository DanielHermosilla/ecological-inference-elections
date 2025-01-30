#include <globals.h>
#include <matrixUtils.h>

/**
 * @brief Main function for generating all combinations possible for a given group and ballot.
 *
 * Given that the combinations are created towards a recursion, the main function is to add as a wrapper towards the
 * recursion function. It will return a pointer with all of the possible configurations. The combinations have an upper
 * constraint given the results that a candidate have gotten.
 *
 * @param[in] b. The index of the ballot box that is going to be calculated.
 * @param[in] totalVotes. The total amount of votes that a group had. For example, if a group had `10` votes, the sum of
 * each element will be 10.
 * @param[in] numCandidates. The amount of candidates. It will also determine the dimensions of the array.
 * @param[in, out] *count. Pointer that will store the total amount of combinations. Useful for iterating over a set.
 *
 * @return size_t **: A pointer that will store arrays of arrays, having all of the possible combinations.
 *
 */
size_t **generateAllConfigurationsPerGroup(int b, int totalVotes, int numCandidates, size_t *count);

Matrix **generateAllConfigurationsPerBallot(int b, int totalBallotVotes, int numCandidates, int numGroups,
                                            size_t *count);
