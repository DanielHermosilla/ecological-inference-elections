#include "hitAndRun.h"
#include "globals.h"
#include <stdio.h>

// ---- Define a structure to store the Omega sets ---- //
typedef struct
{
    uint32_t b;
    Matrix **data;
    size_t size;
} Set;

#define MIN(a, b) ((a) < (b) ? (a) : (b)) // Macro for finding the minimum

Set **OMEGASET = NULL; // Global pointer to store all H sets

/**
 *  @brief Yields an initial point of the polytope given a ballot
 *
 * Given a ballot box index, it returns a matrix of size (GxC) with a valid starting point for the Hit and Run
 * algorithm.
 *
 * @param[in] b The ballot box index
 *
 * @return A matrix with the starting point.
 * Description.
 */
Matrix startingPoint(int b)
{
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *groupVotes = getRow(W, b);
    double *candidateVotes = getColumn(X, b);
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            MATRIX_AT(toReturn, g, c) = MIN(groupVotes[g], candidateVotes[c]);
            groupVotes[g] -= MATRIX_AT(toReturn, g, c);
            candidateVotes[c] -= MATRIX_AT(toReturn, g, c);
        }
    }
    free(groupVotes);
    free(candidateVotes);
    return toReturn;
}

void generateOmegaSet(int M, int S)
{

    int totalSamples = M * S;
    // ---- Allocate memory for the `b` index ----
    OMEGASET = malloc(TOTAL_BALLOTS * sizeof(Set *));

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        OMEGASET[b] = malloc(sizeof(Set));
        // ---- Allocate memory for the `g` index ----
        OMEGASET[b]->b = b;
        OMEGASET[b]->size = totalSamples;
        Matrix startingZ = startingPoint(b);

        for (int s = 0; s < S; s++)
        {
            Matrix steppingZ = startingZ;
            for (int m = 0; m < M; m++)
            {

                // ---- Sample random indexes ---- //
                int groupIndex1 = rand() % TOTAL_GROUPS;
                int groupIndex2 = rand() % (TOTAL_GROUPS - 1);
                if (groupIndex2 >= groupIndex1)
                    groupIndex2++;

                int candidateIndex1 = rand() % TOTAL_CANDIDATES;
                int candidateIndex2 = rand() % (TOTAL_CANDIDATES - 1);
                if (candidateIndex2 >= candidateIndex1)
                    candidateIndex2++;
                // ---...--- //

                double firstSubstraction = MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) - 1;
                double secondSubstraction = MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) - 1;

                if (firstSubstraction < 0 || secondSubstraction < 0)
                    continue;

                MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) -= 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) -= 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex1) += 1;
                MATRIX_AT(steppingZ, groupIndex1, candidateIndex2) += 1;
            }
        }
    }
}

/**
 * @brief Computes the `q` values for all the ballot boxes given a probability matrix. Uses the Hit and Run method.
 *
 * Given a probability matrix with, it returns a flattened array with estimations of the conditional probability. The
 * array can be accesed with the macro `Q_3D` (it's a flattened tensor).
 *
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 *
 * @return A pointer towards the flattened tensor.
 *
 */
double *computeQHitAndRun(Matrix const *probabilities)
{
    generateOmegaSet();
    double *a = malloc(TOTAL_BALLOTS * sizeof(double *));

    return a;
}
