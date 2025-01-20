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

void generateConfigurationsRecursion(int b, size_t *votes, int position, int remainingVotes, int numCandidates,
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
        generateConfigurationsRecursion(b, votes, position + 1, remainingVotes - i, numCandidates, results, count);
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
    generateConfigurationsRecursion(b, votes, 0, totalVotes, numCandidates, &results, count);
    // --- ... --- //
    free(votes);
    return results;
}

Matrix *copyMatrix2(const Matrix *src)
{
    Matrix *dest = malloc(sizeof(Matrix));
    dest->rows = src->rows;
    dest->cols = src->cols;

    // Allocate memory for the data
    dest->data = malloc(src->rows * src->cols * sizeof(size_t));

    // Copy the data
    memcpy(dest->data, src->data, src->rows * src->cols * sizeof(size_t));

    return dest;
}

void generateConfigurationsPerBallotRecursion2(int b, Matrix *votes, int cPosition, int gPosition, int remainingVotes,
                                               int numCandidates, int numGroups, Matrix ***results, size_t *count,
                                               size_t *acumulatedGroupVotes)
{
    size_t *acumulatedCandidateVotes = calloc(numCandidates, sizeof(size_t));
    for (uint16_t g = 0; g < numGroups; g++)
    {
        // Remaining votes must be the total votes that a group did
        size_t **results2 = NULL;
        *count = 0;
        size_t *votes2 = malloc(numCandidates * sizeof(size_t));

        generateConfigurationsRecursion(b, votes2, 0, remainingVotes, numCandidates, &results2, count);

        for (size_t comb = 0; comb < *count; comb++)
        {
            size_t ***validResults = NULL;
            size_t *validCount = 0;
            bool valid = true;
            for (size_t r = 0; r < TOTAL_CANDIDATES; r++)
            {
                if (*results2[r] > (size_t)MATRIX_AT_PTR(X, r, b))
                {
                    valid = false;
                    break;
                }
            }
            if (valid)
            {
                (*validResults) = realloc(*validResults, (*validCount + 1) * sizeof(size_t *));
                (*validResults)[*validCount] = malloc(numCandidates * sizeof(size_t));
                memcpy((*validResults)[*validCount], results2, numCandidates * sizeof(size_t));
                (*validCount)++;
            }
        }
    }
}

void generateConfigurationsPerBallotRecursion(int b, Matrix *votes, int cPosition, int gPosition, int remainingVotes,
                                              int numCandidates, int numGroups, Matrix ***results, size_t *count,
                                              size_t *acumulatedCandidateVotes, size_t *acumulatedGroupVotes)
{
    // ---- Base case: we're on the last coordinate of the matrix ---- //
    if (cPosition == numCandidates - 1 && gPosition == numGroups - 1)
    {
        MATRIX_AT_PTR(votes, gPosition, cPosition) = remainingVotes;

        // ---- Case where the remaining votes doesn't match the last candidate or group total amount of votes
        if (MATRIX_AT_PTR(votes, gPosition, cPosition) > MATRIX_AT_PTR(X, cPosition, b) ||
            MATRIX_AT_PTR(votes, gPosition, cPosition) > MATRIX_AT_PTR(W, b, gPosition) ||
            acumulatedCandidateVotes[cPosition] + remainingVotes != MATRIX_AT_PTR(X, cPosition, b) ||
            acumulatedGroupVotes[gPosition] + remainingVotes != MATRIX_AT_PTR(W, b, gPosition))
        {
            return;
        }
        printf("\nAdding the matrix:\n");
        printMatrix(votes);
        printf("\n");
        acumulatedCandidateVotes[cPosition] += remainingVotes;
        acumulatedGroupVotes[gPosition] += remainingVotes;

        // ---- Store the result ---- //
        // ---- Up to this point, the combination is valid, hence, the results will be stored.
        (*results) = realloc(*results, (*count + 1) * sizeof(size_t *));
        (*results)[*count] = copyMatrix2(votes);
        (*count)++;
        return;
    }

    if (cPosition == numCandidates - 1 && gPosition != numGroups - 1)
    {
        // ---- Case where the total amount of votes from the row doesn't sum up to the total amount of votes the group
        // had ----
        int toAdd = MATRIX_AT_PTR(W, b, gPosition) - acumulatedGroupVotes[gPosition];
        printf("\ntoAdd is %d\nThe current acumulated groupVotes is %d\nThe candidate upper limit is %d\n\n", toAdd,
               (int)acumulatedGroupVotes[gPosition], (int)MATRIX_AT_PTR(X, b, cPosition));
        if (toAdd < 0 || toAdd > MATRIX_AT_PTR(X, b, cPosition))
        {
            printf("\ntoAdd is invalid!\n\n");
            return; // Skip invalid group configurations
        }
        MATRIX_AT_PTR(votes, gPosition, cPosition) = toAdd;
        acumulatedGroupVotes[gPosition] += MATRIX_AT_PTR(votes, gPosition, cPosition);
        // ---- We call the recursion, but this time we move to the next group coordinate ----
        generateConfigurationsPerBallotRecursion(b, votes, 0, gPosition + 1, remainingVotes - toAdd, numCandidates,
                                                 numGroups, results, count, acumulatedCandidateVotes,
                                                 acumulatedGroupVotes);
        return;
    }

    // ---- Loop over all the remaining votes ---- //
    for (int i = 0; i <= remainingVotes; i++)
    { // ---- For each remaining vote
        // ---- Pass the amount of votes to the candidate in the given position ----
        MATRIX_AT_PTR(votes, gPosition, cPosition) = i;
        acumulatedCandidateVotes[cPosition] += i;
        acumulatedGroupVotes[gPosition] += i;

        if (i > MATRIX_AT_PTR(W, b, gPosition) || i > MATRIX_AT_PTR(X, cPosition, b))
        {
            acumulatedCandidateVotes[cPosition] -= i; // Undo addition
            acumulatedGroupVotes[gPosition] -= i;
            return;
        }
        // ---- If the sum of columns or rows are greater than the candidate total votes or group, it's an invalid
        // combination ----
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

        // ---- When the recursion is finished, the loop will advance to the next `i`, hence, we need to remove the
        // acumulated arrays addition ----
        acumulatedCandidateVotes[cPosition] -= i; // Undo addition
        acumulatedGroupVotes[gPosition] -= i;
    }
    // ---...--- //
}

Matrix **generateAllConfigurationsPerBallot(int b, int totalBallotVotes, int numCandidates, int numGroups,
                                            size_t *count)
{
    // ---- Initialize parameters ---- //
    Matrix **results = NULL;
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
