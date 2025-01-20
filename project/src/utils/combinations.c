#include <globals.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

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

double *substractVectors(double *a, double *b, size_t size)
{

    double *toReturn = calloc(size, sizeof(double));
    for (size_t i = 0; i < size; i++)
    {
        toReturn[i] = a[i] - b[i];
    }
    return toReturn;
}

void generateConfigurationsRecursion2(int b, Matrix *votes, double *restriction, int cPosition, int gPosition,
                                      double remainingVotes, int numCandidates, int numGroups, Matrix ***results,
                                      size_t *count)
{
    // sleep(1);
    //  ---- Base case: we're on the last candidate ---- //
    if (cPosition == numCandidates - 1)
    {
        // ---- Assign remaining votes to the last candidate ----
        MATRIX_AT_PTR(votes, gPosition, cPosition) = remainingVotes;

        // ---- If the last candidate actually had less votes, ditch that combination ----
        if (MATRIX_AT_PTR(votes, gPosition, cPosition) > restriction[cPosition])
        {
            // ---- Exit the recursion and don't save anything
            return;
        }
        double rowSum = 0;
        for (int i = 0; i < numCandidates; i++)
        {
            rowSum += MATRIX_AT_PTR(votes, gPosition, i);
        }
        if (rowSum != MATRIX_AT_PTR(W, b, gPosition))
        {
            return; // Invalid row configuration
        }

        // ---- Store the result ---- //
        // ---- Up to this point, the matrix should be valid, hence, the results will be stored.
        if (gPosition == numGroups - 1)
        {

            for (int k = 0; k < numCandidates; k++)
            {
                double colSum = 0;
                for (int i = 0; i < numGroups; i++)
                {
                    colSum += MATRIX_AT_PTR(votes, i, k);
                }
                if (colSum != MATRIX_AT_PTR(X, k, b))
                    return;
            }
            (*results) = realloc(*results, (*count + 1) * sizeof(Matrix *));
            (*results)[*count] = votes;
            (*count)++;
            return;
        }
        else
        {
            double *newColumnRestriction = calloc(numCandidates, sizeof(double));
            for (int i = 0; i < numCandidates; i++)
            {
                newColumnRestriction[i] = restriction[i] - MATRIX_AT_PTR(votes, gPosition, i);
            }
            generateConfigurationsRecursion2(b, votes, newColumnRestriction, 0, gPosition + 1,
                                             MATRIX_AT_PTR(W, b, gPosition + 1), numCandidates, numGroups, results,
                                             count);
            free(newColumnRestriction);
            return;
        }
        // ---...--- //
    }
    // ---...--- //

    // ---- Loop over all the remaining votes ---- //

    for (double i = 0; i <= remainingVotes; i++)
    { // ---- For each remaining vote
        // ---- Assing that amount of votes to the candidate in the given position ----
        MATRIX_AT_PTR(votes, gPosition, cPosition) = i;

        // ---- If the candidate actually had less votes, ditch that combination ----
        if (MATRIX_AT_PTR(votes, gPosition, cPosition) > restriction[cPosition])
        {
            // ---- Exit the recursion and dont save anything
            return;
        }
        // ---- Call the recursion ----
        generateConfigurationsRecursion2(b, votes, restriction, cPosition + 1, gPosition, remainingVotes - i,
                                         numCandidates, numGroups, results, count);
    }
    // ---...--- //
}

/*
void generateConfigurationsPerBallotRecursion2(int b, Matrix *votes, int cPosition, int gPosition, int remainingVotes,
                                               int numCandidates, int numGroups, Matrix ***results, size_t *count,
                                               size_t *acumulatedGroupVotes)
{

    double **overallResults = NULL;
    size_t *overallCounter = 0;

    double *pastRestriction = getColumn(X, b);
    double **groupResults = NULL;
    double *initialArray = malloc(numCandidates * sizeof(double));
    double *currentCombination = calloc(numCandidates, sizeof(double));
    size_t *amountOfCombinations = 0;
    double *groupRestriction = substractVectors(pastRestriction, currentCombination, numCandidates);

    generateConfigurationsRecursion2(b, initialArray, pastRestriction, 0, MATRIX_AT_PTR(W, b, 0), numCandidates,
                                     &groupResults, amountOfCombinations);

    for (size_t k = 0; k < *amountOfCombinations; k++)
    {
        currentCombination = groupResults[k];
        for (uint16_t g = 1; g < numGroups; g++)
        {
            groupRestriction = substractVectors(pastRestriction, currentCombination, numCandidates);
        }
    }

    double *acumulatedCandidateVotes = calloc(numCandidates, sizeof(double));
    double **results2 = NULL;
    double *votes2 = malloc(numCandidates * sizeof(double));
    *count = 0;
    double *firstRestriction = getColumn(X, b);
    generateConfigurationsRecursion2(b, votes2, firstRestriction, 0, MATRIX_AT_PTR(W, b, 0), numCandidates, &results2,
                                     count);
    for (size_t k = 0; k < *count; k++)
    {
        double *currentCombination = results2[k];
        double *currentConstraint = substractVectors() for (uint16_t g = 1; g < numGroups; g++)
        {
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
*/
Matrix **generateAllConfigurationsPerBallot(int b, int totalBallotVotes, int numCandidates, int numGroups,
                                            size_t *count)
{
    // ---- Initialize parameters ---- //
    Matrix **results = NULL;
    *count = 0;
    Matrix votes = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *initialRestriction = getColumn(X, b);
    // --- ... --- //

    // ---- Call the recursion ---- //
    generateConfigurationsRecursion2(b, &votes, initialRestriction, 0, 0, MATRIX_AT_PTR(W, b, 0), numCandidates,
                                     numGroups, &results, count);

    // --- ... --- //
    freeMatrix(&votes);
    return results;
}
