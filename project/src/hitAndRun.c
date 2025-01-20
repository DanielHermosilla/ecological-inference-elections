#include "hitAndRun.h"
#include <stdio.h>

// ---- Define a structure to store the Omega sets ---- //
typedef struct
{
    uint32_t b;
    Matrix **data;
    size_t size;
} Set;

Set **OMEGASET = NULL; // Global pointer to store all H sets

void generateOmegaSet()
{
    // ---- Allocate memory for the `b` index ----
    OMEGASET = malloc(TOTAL_BALLOTS * sizeof(Set *));

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box
        OMEGASET[b] = malloc(sizeof(Set));
        // ---- Allocate memory for the `g` index ----
        OMEGASET[b]->b = b;

        // ---- Parameters for the function ----
        size_t totalVotes = BALLOTS_VOTES[b];

        // --- ... --- //

        // ---- Compute the set combinations ---- //
        size_t count = 0;
        printf("\nRunning for ballot:\t%d\n", b);
        Matrix **configurations =
            generateAllConfigurationsPerBallot(b, totalVotes, TOTAL_CANDIDATES, TOTAL_GROUPS, &count);

        // ---- Store configurations and size ----
        OMEGASET[b]->data = configurations;
        OMEGASET[b]->size = count;
        printf("Printing the generated configurations:\n");
        for (size_t a = 0; a < count; a++)
        {
            printMatrix(OMEGASET[b]->data[a]);
            printf("\n");
        }
        // --- ... --- //
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
