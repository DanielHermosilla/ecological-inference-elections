#include <globals.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/**
 * @brief Calculate the possible configurations via a recursion.
 *
 * Given the initial parameters, it will calculate an array with all of the possible voting outcomes, having the
 * constraint on the votes that a candidate had scored. The constraint will put an upper limit towards the dimension of
 * the given candidate.
 *
 * @param[in] b. The index of the ballot box that is going to be calculated.
 * @param[in] *votes. The array that will be storing every vote.
 * @param[in] position. The current position that indicates the index of the candidate that's being calculated.
 * @param[in] remainingVotes. The remaining votes to distribute over each dimension.
 * @param[in] numCandidates. The amount of candidates. It will also determine the dimensions of the array.
 * @param[in, out] ***results. The main array for storing each possible combination.
 * @param[in, out] *count. A counter that will have the amount of combinations
 *
 * @return void.
 *
 */

void generateConfigurations(int b, size_t *votes, int position, int remainingVotes, int numCandidates,
                            size_t ***results, size_t *count)
{
    // ---- Base case: we're on the last candidate ---- //
    if (position == numCandidates - 1)
    {
        // ---- Assign remaining votes to the last candidate ----
        votes[position] = remainingVotes;

        // ---- If the last candidate actually had less votes, ditch that combination ----
        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // ---- Exit the recursion and don't save anything
            return;
        }

        // ---- Store the result ---- //
        // ---- Up to this point, the combination is valid, hence, the results will be stored.
        (*results) = realloc(*results, (*count + 1) * sizeof(size_t *));
        (*results)[*count] = malloc(numCandidates * sizeof(size_t));
        memcpy((*results)[*count], votes, numCandidates * sizeof(size_t));
        (*count)++;
        return;
        // ---...--- //
    }
    // ---...--- //

    // ---- Loop over all the remaining votes ---- //
    for (int i = 0; i <= remainingVotes; i++)
    { // ---- For each remaining vote
        // ---- Assing that amount of votes to the candidate in the given position ----
        votes[position] = i;

        // ---- If the candidate actually had less votes, ditch that combination ----
        if (votes[position] > MATRIX_AT_PTR(X, position, b))
        {
            // ---- Exit the recursion and dont save anything
            return;
        }
        // ---- Call the recursion ----
        generateConfigurations(b, votes, position + 1, remainingVotes - i, numCandidates, results, count);
    }
    // ---...--- //
}

/**
 * @brief Main function for generating all combinations possible.
 *
 * Given that the combinations are created towards a recursion, the main function is to add as a wrapper towards the
 * recursion function. It will return a pointer with all of the possible configurations. The combinations have an upper
 * constraint given the results that a candidate have gotten.
 *
 * @param[in] b. The index of the ballot box that is going to be calculated.
 * @param[in] totalVotes. The total amount of votes to handle. For example, if a group had `10` votes, the sum of each
 * element will be 10.
 * @param[in] numCandidates. The amount of candidates. It will also determine the dimensions of the array.
 * @param[in, out] *count. Pointer that will store the total amount of combinations. Useful for iterating over a set.
 *
 * @return size_t **: A pointer that will store arrays of arrays, having all of the possible combinations.
 *
 */
size_t **generateAllConfigurationsPerGroup(int b, int totalVotes, int numCandidates, size_t *count)
{
    // ---- Initialize parameters ---- //
    size_t **results = NULL;
    *count = 0;
    size_t *votes = malloc(numCandidates * sizeof(size_t));
    // --- ... --- //

    // ---- Call the recursion ---- //
    generateConfigurations(b, votes, 0, totalVotes, numCandidates, &results, count);
    // --- ... --- //
    free(votes);
    return results;
}

void generateConfigurationsPerBallotRecursion(int b, Matrix *votes, int cPosition, int gPosition, int remainingVotes,
                                              int numCandidates, int numGroups, size_t ***results, size_t *count,
                                              size_t *acumulatedCandidateVotes, size_t *acumulatedGroupVotes)
{
    // ---- Base case: we're on the last coordinate ---- //
    if (cPosition == numCandidates - 1 && gPosition == numGroups - 1)
    {
        MATRIX_AT_PTR(votes, gPosition, cPosition) = remainingVotes;

        if (MATRIX_AT_PTR(votes, gPosition, cPosition) > MATRIX_AT_PTR(X, cPosition, b) ||
            MATRIX_AT_PTR(votes, gPosition, cPosition) > MATRIX_AT_PTR(W, b, gPosition) ||
            acumulatedCandidateVotes[cPosition] + remainingVotes != MATRIX_AT_PTR(X, cPosition, b) ||
            acumulatedGroupVotes[gPosition] + remainingVotes != MATRIX_AT_PTR(W, b, gPosition))
        {
            return;
        }
        acumulatedCandidateVotes[cPosition] += remainingVotes;
        acumulatedGroupVotes[gPosition] += remainingVotes;

        // ---- Store the result ---- //
        // ---- Up to this point, the combination is valid, hence, the results will be stored.
        (*results) = realloc(*results, (*count + 1) * sizeof(size_t *));
        (*results)[*count] = malloc(sizeof(Matrix));
        memcpy((*results)[*count], votes, sizeof(Matrix));
        (*count)++;
        return;
    }

    if (cPosition == numCandidates && gPosition != numGroups - 1)
    {
        if (acumulatedGroupVotes[gPosition] != MATRIX_AT_PTR(W, b, gPosition))
        {
            return; // Skip invalid group configurations
        }
        generateConfigurationsPerBallotRecursion(b, votes, 0, gPosition + 1, remainingVotes, numCandidates, numGroups,
                                                 results, count, acumulatedCandidateVotes, acumulatedGroupVotes);
        return;
    }

    // ---- Loop over all the remaining votes ---- //
    for (int i = 0; i <= remainingVotes && i <= MATRIX_AT_PTR(W, b, gPosition) && i <= MATRIX_AT_PTR(X, cPosition, b);
         i++)
    { // ---- For each remaining vote
        // ---- Assing that amount of votes to the candidate in the given position ----
        MATRIX_AT_PTR(votes, gPosition, cPosition) = i;
        acumulatedCandidateVotes[cPosition] += i;
        acumulatedGroupVotes[gPosition] += i;

        if (acumulatedCandidateVotes[cPosition] > MATRIX_AT_PTR(X, cPosition, b) ||
            acumulatedGroupVotes[gPosition] > MATRIX_AT_PTR(W, b, gPosition))
        {
            acumulatedCandidateVotes[cPosition] -= i; // Undo addition
            acumulatedGroupVotes[gPosition] -= i;
            return;
        }
        // ---- Call the recursion ----
        generateConfigurationsPerBallotRecursion(b, votes, cPosition + 1, gPosition, remainingVotes - i, numCandidates,
                                                 numGroups, results, count, acumulatedCandidateVotes,
                                                 acumulatedGroupVotes);

        acumulatedCandidateVotes[cPosition] -= i; // Undo addition
        acumulatedGroupVotes[gPosition] -= i;
    }
    // ---...--- //
}

size_t **generateAllConfigurationsPerBallot(int b, int totalBallotVotes, int numCandidates, int numGroups,
                                            size_t *count)
{
    // ---- Initialize parameters ---- //
    size_t **results = NULL;
    *count = 0;
    Matrix votes = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    size_t *acumulatedCandidateVotes = calloc(numCandidates, sizeof(size_t));
    size_t *acumulatedGroupVotes = calloc(numGroups, sizeof(size_t));
    // --- ... --- //

    // ---- Call the recursion ---- //
    generateConfigurationsPerBallotRecursion(b, &votes, 0, 0, totalBallotVotes, numCandidates, numGroups, &results,
                                             count, acumulatedCandidateVotes, acumulatedGroupVotes);
    // --- ... --- //
    freeMatrix(&votes);
    free(acumulatedCandidateVotes);
    free(acumulatedGroupVotes);
    return results;
}
