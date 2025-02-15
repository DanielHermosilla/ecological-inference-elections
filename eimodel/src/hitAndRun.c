/*
Copyright (c) 2025 Daniel Hermosilla

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "hitAndRun.h"
#include <R_ext/Memory.h>
#include <R_ext/Random.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

OmegaSet **OMEGASET = NULL; // Global pointer to store all H sets
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
    // ---- Retrieve the initial variables ---- //
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *groupVotes = getRow(W, b);
    double *candidateVotes = getColumn(X, b);
    // ---...--- //
    // ---- Compute the main loop ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // --- For each group
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // --- For each candidate given a group
            MATRIX_AT(toReturn, g, c) = MIN(groupVotes[g], candidateVotes[c]);
            groupVotes[g] -= MATRIX_AT(toReturn, g, c);
            candidateVotes[c] -= MATRIX_AT(toReturn, g, c);
        }
    }
    // ---...--- //
    Free(groupVotes);
    Free(candidateVotes);
    return toReturn;
}

/*
 * @brief Precomputes the sets used for the simulation.
 *
 * Precomputes the sets that are independent from each EM iteration. It is made with parallelism (NOT SUPPORTED) towards
 * the ballot boxes and with a static assignment for ensuring reproducibility.
 *
 * @param[in] M. The step size between consecutive samples. Note that the direction is assigned randomly.
 * @param[in] S. The amount of samples for each ballot box.
 * @param[in] seedNum. An arbitrary number to seed the process.
 *
 * @return void. Written on the global variable.
 */
void generateOmegaSet(int M, int S, unsigned int seedNum)
{
    // ---- Allocate memory for the `b` index ----
    OMEGASET = Calloc(TOTAL_BALLOTS, OmegaSet *);

    /*
     * Currently, CRAN doesn't allow to use the own system random number generator.
     * Even though it has its upsides (setting a global seed from R), it doesn't allow to use
     * parallelization with RNG's in a thread-safe way. In fact, R's seed is global and unique.
     * The C's approach back then, was to define a seed per thread;
     *
     * int seed = rand_r(&seedNum) + omp_get_thread_number
     *
     * This cannot be done anymore, so, as of now, this loop cannot be parallelized, despite its
     * benefits.
     */

    GetRNGstate();

    // ---- Perform the main iterations ---- //
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {                   // ---- For every ballot box
        if (b % 5 == 0) // Checks condition every 5 iterations
            R_CheckUserInterrupt();
        // ---- Define a seed, that will be unique per thread ----
        //    unsigned int seed = rand_r(&seedNum) + omp_get_thread_number();
        // ---- Allocate memory for the OmegaSet ---- //
        OMEGASET[b] = Calloc(1, OmegaSet);
        OMEGASET[b]->b = b;
        OMEGASET[b]->size = S;
        OMEGASET[b]->data = Calloc(S, Matrix *);
        // ---...--- //
        // ---- The `base` element used as a starting point ----
        Matrix startingZ = startingPoint(b);

        for (int s = 0; s < S; s++)
        { // --- For each sample given a ballot box
            // ---- Copy the initial matrix ----
            Matrix steppingZ = copyMatrix(&startingZ);
            for (int m = 0; m < M; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int groupIndex1 = (int)(unif_rand() * TOTAL_GROUPS);
                int groupIndex2;
                do
                {
                    groupIndex2 = (int)(unif_rand() * TOTAL_GROUPS);
                } while (groupIndex2 == groupIndex1);

                int candidateIndex1 = (int)(unif_rand() * TOTAL_CANDIDATES);

                int candidateIndex2;
                do
                {
                    candidateIndex2 = (int)(unif_rand() * TOTAL_CANDIDATES);
                } while (candidateIndex2 == candidateIndex1);
                // ---...--- //

                // ---- Check non negativity condition ---- //
                double firstSubstraction = MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) - 1;
                double secondSubstraction = MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) - 1;

                if (firstSubstraction < 0 || secondSubstraction < 0)
                    continue;
                // ---...--- //

                // ---- Asign changes on the new matrix ---- //
                MATRIX_AT(steppingZ, groupIndex1, candidateIndex1) -= 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex2) -= 1;
                MATRIX_AT(steppingZ, groupIndex1, candidateIndex2) += 1;
                MATRIX_AT(steppingZ, groupIndex2, candidateIndex1) += 1;
                // ---...--- //
            } // --- End the step size loop
            // ---- Add the combination to the OmegaSet ---- //
            Matrix *append = Calloc(1, Matrix);
            *append = copyMatrix(&steppingZ);
            OMEGASET[b]->data[s] = append;
            freeMatrix(&steppingZ);
            // ---...--- //
        } // --- End the sample loop
        freeMatrix(&startingZ);
    } // --- End the ballot box loop
    PutRNGstate();
}

/**
 * @brief Computes the pre-computable values of the expression that doesn't depend on EM iterations
 *
 * Given a ballot box index and a matrix represent an element from the Hit and Run OmegaSet, it computes the following:
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
        // result += gsl_sf_lngamma((int)MATRIX_AT_PTR(W, b, g) + 1);
        result += lgamma1p((int)MATRIX_AT_PTR(W, b, g) + 1);

        for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate
            // ---- Divide by each h_i! ----
            // result -= gsl_sf_lngamma(MATRIX_AT_PTR(currentMatrix, g, i) + 1);
            result -= lgamma1p(MATRIX_AT_PTR(currentMatrix, g, i) + 1);
        }
    }
    // ---- Return the original result by exponentiating ----
    return exp(result);
}

/**
 * @brief Calculates the last term of the multiplication OmegaSet
 *
 * Given a probability matrix, a ballot index and a OmegaSet index, it calculates:
 *
 * $$\Prod_{g\in G}\Prod_{c\in C}p_{gc}^{z_{bgc}}$$
 *
 * @param[in] *probabilities A pointer toward the probabilities Matrix.
 * @param[in] b The index of the ballot box
 * @param[in] setIndex The index of the OmegaSet
 *
 * @return double: The result of the product
 *
 */
double logarithmicProduct(const Matrix *probabilities, const int b, const int setIndex)
{
    // ---- Define initial parameters ---- //
    double log_result = 0;
    Matrix *currentMatrix = OMEGASET[b]->data[setIndex];
    // ---...--- //
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            log_result += MATRIX_AT_PTR(currentMatrix, g, c) * log(MATRIX_AT_PTR(probabilities, g, c));
        }
    }
    // --- ... --- //
    // ---- Exponetiate the final result ----
    return exp(log_result);
}

/**
 * @brief Precomputes the multinomial multiplication that is independent for each EM iteration.
 *
 * Calls the main function for computing all of the calculations related with the final result that are independent from
 * each EM call. Specifically, for each ballot box and its simulations, the following is calculated:
 *
 * $$\Prod_{g\in G}\binom{w_{bg}}{z_{bg1},\cdots, z_{bgC}}$$
 *
 * It transform the main product and the factorial to logarithmic scale for making efficient calculations.
 *
 * @return. Results written at the global variable
 */
void preComputeMultinomial()
{
    // ---- Initialize space for storing all of the simulations ---- //
    multinomialVals = Calloc(TOTAL_BALLOTS, double **);
    // ---...--- //
    // ---- Compute the simulated combinations for each OmegaSet ---- //
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        // ---- Define the current OmegaSet and allocate memory for saving its size ----
        OmegaSet *currentSet = OMEGASET[b];
        multinomialVals[b] = Calloc(currentSet->size, double);
        for (size_t s = 0; s < currentSet->size; s++)
        { // --- For each simulation given a balot box
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
 * @param[in] params A QMethodInput struct with the `M` (step size) and `S` (samples) parameters
 *
 * @return A pointer towards the flattened tensor.
 *
 */
double *computeQHitAndRun(Matrix const *probabilities, QMethodInput params)
{
    // ---- Compute the variables that can be reused ---- //
    if (OMEGASET == NULL)
    {
        generateOmegaSet(params.M, params.S, 2);
    }
    if (multinomialVals == NULL)
    {
        preComputeMultinomial();
    }
    // ---...--- //

    // ---- Compute the final values and fill the returning array ---- //
    double *array2 = (double *)Calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, double); // Array to return
    // ---- Use a static assignment since the workload is even between threads ----

#ifdef _OPENMP
#pragma omp parallel for collapse(3) schedule(static)
#endif
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        OmegaSet *currentSet = OMEGASET[b];
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box
                // ---- If there's a division by zero, assign probability 0 and skip the iteration ----
                if (MATRIX_AT_PTR(W, b, g) == 0)
                {
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                    continue;
                }
                // ---- Initialize variables for the multiplications. The terms are the outer parenthesis. ----
                double firstTerm = 0;
                double secondTerm = 0;
                for (size_t s = 0; s < currentSet->size; s++)
                { // --- For each sample given a ballot box, group and candidate
                    // ---- Perform the summatory ---- //
                    Matrix *currentMatrix = currentSet->data[s];
                    double multiplications = logarithmicProduct(probabilities, b, s) * multinomialVals[b][s];
                    firstTerm += multiplications;
                    secondTerm += multiplications * (MATRIX_AT_PTR(currentMatrix, g, c) / MATRIX_AT_PTR(W, b, g));
                    // ---...--- //
                }
                // ---- Add the element to the array ---- //
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = (1 / firstTerm) * secondTerm;
            }
        }
    }
    return array2;
}

void cleanHitAndRun()
{
    if (OMEGASET != NULL && multinomialVals != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        { // --- For each ballot box
            for (size_t s = 0; s < OMEGASET[b]->size; s++)
            {                                     // For each sample given a ballot box
                freeMatrix(OMEGASET[b]->data[s]); // Free individual matrices
                Free(OMEGASET[b]->data[s]);       // Free the pointers to matrices
            }
            Free(OMEGASET[b]->data); // Free the data array
            Free(OMEGASET[b]);       // Free the OmegaSet struct
            Free(multinomialVals[b]);
        }
        Free(multinomialVals); // Free the precomputed multinomial values
        Free(OMEGASET);        // Free the OMEGASET array
        OMEGASET = NULL;
        multinomialVals = NULL;
    }
}
//__attribute__((destructor)) void cleanEverything()
//{

//   cleanHitAndRun();
//}
