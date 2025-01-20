#include "hitAndRun.h"
#include "globals.h"
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdio.h>
#include <sys/_types/_size_t.h>

// ---- Define a structure to store the Omega sets ---- //
typedef struct
{
    uint32_t b;
    Matrix **data;
    size_t size;
} Set;

// Macro for finding the minimum
#define MIN(a, b) ((a) < (b) ? (a) : (b))

Set **OMEGASET = NULL; // Global pointer to store all H sets
double **multinomialVals = NULL;

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

    // ---- Allocate memory for the `b` index ----
    OMEGASET = malloc(TOTAL_BALLOTS * sizeof(Set *));

#pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For every ballot box

        // ---- Allocate memory for the set ---- //
        OMEGASET[b] = malloc(sizeof(Set));
        OMEGASET[b]->b = b;
        OMEGASET[b]->size = S;
        OMEGASET[b]->data = malloc(S * sizeof(Matrix *));
        // ---...--- //
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

                // ---- Check non negativity condition ---- //
                double firstSubstraction = MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) - 1;
                double secondSubstraction = MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) - 1;

                if (firstSubstraction < 0 || secondSubstraction < 0)
                    continue;
                // ---...--- //
                MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) -= 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) -= 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex1) += 1;
                MATRIX_AT(steppingZ, groupIndex1, candidateIndex2) += 1;
            }
            // ---- Add the combination to the set ---- //
            Matrix *append = malloc(sizeof(Matrix));
            *append = copyMatrix(&steppingZ);
            OMEGASET[b]->data[s] = append;
            // ---...--- //
        }
    }
}

/**
 * @brief Computes the pre-computable values of the expression that doesn't depend on EM iterations
 *
 * Given a ballot box index and a matrix represent an element from the Hit and Run set, it computes the following:
 *
 * $$\Prod_{g'\in G}\binom{w_{bg'}}{z_{bg'1}\cdots z_{bg'C}}$$
 *
 * @param[in] b The ballot box index
 * @param[in] *currentMatrix A pointer towards the current matricial element, of size GxC.
 *
 * double The result of the calculation.
 */
double multinomialCoeff(const int b, Matrix *currentMatrix)
{
    double result = 0;
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        // --- Compute ln(w_bf!). When adding by one, it considers the last element too ---
        result += gsl_sf_lngamma((int)MATRIX_AT_PTR(W, b, g) + 1);

        for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate
            // ---- Divide by each h_i! ----
            result -= gsl_sf_lngamma(MATRIX_AT_PTR(currentMatrix, g, i) + 1);
        }
    }
    // ---- Return the original result by exponentiating ----
    return exp(result);
}

/**
 * @brief Calculate the chained product between probabilities as defined in `a`:
 *
 * Given an `H` element and the `f` index it computes the chained product. It will use logarithmics for reducing
 * complexity.
 *
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] setIndex The index of the set
 *
 * @return double: The result of the product
 *
 *
 */
double prod(const Matrix *probabilities, const int b, const int g, const int setIndex)
{
    double log_result = 0;
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate

        double currentElement = MATRIX_AT_PTR(OMEGASET[b]->data[setIndex], g, c);
        // ---- Get the matrix value at `f` and `c` ----
        double prob = MATRIX_AT_PTR(probabilities, g, c);

        // ---- If the multiplication gets to be zero, it would be undefined in the logarithm (and zero for all of the
        // chain). Do an early stop if that happens ----
        if (prob == 0.0 && currentElement > 0)
        {
            return 0.0; // ---- Early stop ----
        }

        // ---- Ensures the probability is greater than zero for not getting a undefined logarithm. Anyways, that
        // shouldn't happen. ----
        if (prob > 0)
        {
            // ---- Add the result to the logarithm. Remember that by logarithm properties, the multiplication of the
            // arguments goes as a summatory ----
            log_result += currentElement * log(prob);
        }
    }
    // --- ... --- //

    // ---- Exponetiate the final result ----
    return exp(log_result);
}

void preComputeMultinomial()
{

    multinomialVals = malloc(TOTAL_BALLOTS * sizeof(double **));
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        Set *currentSet = OMEGASET[b];
        multinomialVals[b] = malloc(currentSet->size * sizeof(double));
        for (size_t s = 0; s < currentSet->size; s++)
        {
            multinomialVals[b][s] = multinomialCoeff(b, currentSet->data[s]);
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
double *computeQHitAndRun(Matrix const *probabilities, int M, int S)
{
    generateOmegaSet(M, S);
    double *a = malloc(TOTAL_BALLOTS * sizeof(double *));

    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        Set *currentSet = OMEGASET[b];
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                for (size_t s = 0; s < currentSet->size; s++)
                {
                    Matrix *currentMatrix = currentSet->data[s];
                    double log_product = 1;
                    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
                    {
                    }
                }

                // ---- Exponetiate the final result ----
            }
        }
    }
    return a;
}

__attribute__((destructor)) void cleanOmega()
{
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (size_t s = 0; s < OMEGASET[b]->size; s++)
        {
            freeMatrix(OMEGASET[b]->data[s]); // Free individual matrices
            free(OMEGASET[b]->data[s]);       // Free the pointers to matrices
        }
        free(OMEGASET[b]->data); // Free the data array
        free(OMEGASET[b]);       // Free the Set struct
    }
    free(OMEGASET); // Free the OMEGASET array
}
