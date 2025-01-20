#include "hitAndRun.h"
#include "globals.h"
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdio.h>
#include <sys/_types/_size_t.h>
#include <unistd.h>

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

    // #pragma omp parallel for
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
            Matrix steppingZ = copyMatrix(&startingZ);
            for (int m = 0; m < M; m++)
            {
                // ---- Sample random indexes ---- //
                int groupIndex1 = rand() % TOTAL_GROUPS;
                int groupIndex2;
                do
                {
                    groupIndex2 = rand() % TOTAL_GROUPS;
                } while (groupIndex2 == groupIndex1);

                int candidateIndex1 = rand() % TOTAL_CANDIDATES;
                int candidateIndex2;
                do
                {
                    candidateIndex2 = rand() % TOTAL_CANDIDATES;
                } while (candidateIndex2 == candidateIndex1);
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
            freeMatrix(&steppingZ);
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
double preMultinomialCoeff(const int b, Matrix *currentMatrix)
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
 * @brief Calculates the last term of the multiplication set
 *
 * Given a probability matrix, a ballot index and a set index, it calculates:
 *
 * $$\Prod_{g\in G}\Prod_{c\in C}p_{gc}^{z_{bgc}}$$
 *
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] b The index of the ballot box
 * @param[in] setIndex The index of the set
 *
 * @return double: The result of the product
 *
 *
 */
double logarithmicProduct(const Matrix *probabilities, const int b, const int setIndex)
{
    double log_result = 0;
    Matrix *currentMatrix = OMEGASET[b]->data[setIndex];

    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            log_result += MATRIX_AT_PTR(currentMatrix, g, c) * log(MATRIX_AT_PTR(probabilities, g, c));
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
            multinomialVals[b][s] = preMultinomialCoeff(b, currentSet->data[s]);
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
    srand(42);

    if (OMEGASET == NULL)
        generateOmegaSet(M, S);
    if (multinomialVals == NULL)
        preComputeMultinomial();

    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        Set *currentSet = OMEGASET[b];
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {

                if (MATRIX_AT_PTR(W, b, g) == 0) // Handle division by zero.
                {
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                    continue;
                }
                double firstTerm = 0;
                double secondTerm = 0;
                for (size_t s = 0; s < currentSet->size; s++)
                {
                    Matrix *currentMatrix = currentSet->data[s];
                    double multiplications = logarithmicProduct(probabilities, b, s) * multinomialVals[b][s];
                    firstTerm += multiplications;
                    secondTerm += multiplications * (MATRIX_AT_PTR(currentMatrix, g, c) / MATRIX_AT_PTR(W, b, g));
                }
                // printf("\nAdding the element %.4f on iteration b=%d, c=%d and g=%d\n", (1 / firstTerm) * secondTerm,
                // b,
                //     c, g);
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = (1 / firstTerm) * secondTerm;
            }
        }
    }
    return array2;
}

__attribute__((destructor)) void cleanOmega()
{
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (size_t s = 0; s < OMEGASET[b]->size; s++)
        {
            freeMatrix(OMEGASET[b]->data[s]); // Free individual matrices
            free(OMEGASET[b]->data[s]);       // Free the pointers to matrices
            ;
        }
        free(OMEGASET[b]->data); // Free the data array
        free(OMEGASET[b]);       // Free the Set struct
        free(multinomialVals[b]);
    }
    free(multinomialVals);
    free(OMEGASET); // Free the OMEGASET array
}
