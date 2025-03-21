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
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h> // for R_CheckUserInterrupt()
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

OmegaSet **OMEGASET = NULL; // Global pointer to store all H sets
double **multinomialVals = NULL;
double *logGammaArr = NULL;
double *loglogGammaArr = NULL;

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

void allocateRandoms(int M, int S, char *array)
{
}
/*
 * @brief Precomputes the sets used for the simulation.
 *
 * Precomputes the sets that are independent from each EM iteration. It is made with parallelism (NOT SUPPORTED) towards
 * the ballot boxes and with a static assignment for ensuring reproducibility.
 *
 * @param[in] M. The step size between consecutive samples. Note that the direction is assigned randomly.
 * @param[in] S. The amount of samples for each ballot box.
 *
 * @return void. Written on the global variable.
 */
void generateOmegaSet(int M, int S)
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

    /*
     * UPDATE:
     * It may be useful but dangerous to create a RNG pool before looping. It may be much faster, however
     * the randomn numbers are dinamically asigned.
     *
     * The dynamic condition that creates a RNG call is g1 == g2 (with replacement) and c1 == c2. Both expected
     * values are 1/|G| and 1/|C|. The idea would be to generate a sample pool such that 97.5% of times it doesn't
     * crash. Hence, one loop through M uses on average 2 + (1/|C| + 1/|G| - 1/|C| * 1/|G|) -> failing probability.
     * Given that, the expected value of RNG calls  is |B| * |S| * |M| * (2 + |C|^-1 + |G|^-1 - (|C|*|G|)^-1)
     *
     * On the other side, let X and Y be the failing d.r.v of |C| and |G| respectively. Then, Var(X+Y)=(|G|*|G-1| +
     * |C|*|C-1|)/12. Then, since |B|, |S|, |M| and |2| are constants:
     *
     * Var( |B| * |S| * |M| * (2 + X + Y)) = |B|^2 * |S|^2 * |M|^2 Var(2 + X + Y)
     *  								   = |B|^2 * |S|^2 * |M|^2 * (|G|*|G-1| + |C|*|C-1|)/12.
     *  							std()  = |B| * |S| * |M| sqrt( [ |G|*|G-1| + |C|*|C-1| ] / 12 )
     *
     * So, ensuring a \mu + 2\sigma approach...,
     *
     * Pool size = |B| * |S| * |M| * (2 + |C|^-1 + |G|^-1 - (|C|*|G|)^-1) + 2|B| * |S| * |M| sqrt( [ |G|*|G-1| +
     * |C|*|C-1| ] / 12 ) = |B| * |S| * |M| * (2 + |C|^-1 + |G|^-1 - (|C|*|G|)^-1 + 2 * sqrt( [ |G|*|G-1| + |C|*|C-1| ]
     * / 12)
     *
     * Is this feasible? Let's suppose 50 ballot boxes, 1000 samples and a step size of 3000. Each of the RNG calls are
     * stored as doubles (64bits) aswell. That's 1.2GB, just by considering the first factor and not accounting the |C|
     * and |G| terms. If a pool cannot be allocated before the loop, maybe it could be possible to execute within
     * distinct pools by dividing the ballot boxes and free within each big iteration. However, that's not so practical
     * to implement, and the 2.5% failure rate still bothers me.
     *
     * UPDATE: Tried diverses takes, with chunk, pooling, etc. None of them worked since samples were too similar...
     */

    GetRNGstate();

    // ---- Perform the main iterations ---- //
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {                               // ---- For every ballot box
        if (b % 5 == 0)             // Checks condition every 5 iterations
            R_CheckUserInterrupt(); // This might be fatal, since it doesn't free global memory
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
        {                                 // --- For each sample given a ballot box
            if (S > 5000 && S % 100 == 0) // If there's a big amount of samples, check interrupts
                R_CheckUserInterrupt();
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
        result += lgamma1p((int)MATRIX_AT_PTR(W, b, g)); // TODO: Could be saved

        for (uint16_t i = 0; i < TOTAL_CANDIDATES; i++)
        { // ---- For each candidate
            // ---- Divide by each h_i! ----
            // result -= gsl_sf_lngamma(MATRIX_AT_PTR(currentMatrix, g, i) + 1);
            result -= lgamma1p(MATRIX_AT_PTR(currentMatrix, g, i));
        }
    }
    // ---- Return the original result by exponentiating ----
    return result; // TODO: This can be avoided!
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
    // TODO: logPRod
    // ---- Define initial parameters ---- //
    double log_result = 0;
    Matrix *currentMatrix = OMEGASET[b]->data[setIndex];
    // ---...--- //
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // Cambiar a log(1) if probabilidad 0
            // Producto punto
            log_result += MATRIX_AT_PTR(currentMatrix, g, c) * log(MATRIX_AT_PTR(probabilities, g, c));
        }
    }
    // --- ... --- //
    // ---- Exponetiate the final result ----
    return log_result;
}

void precomputeLogGammas()
{
    // We must get the biggest W_{bg}
    int biggestW = (int)maxElement(W);
    logGammaArr = (double *)Calloc(biggestW + 1, double);
    loglogGammaArr = (double *)Calloc(biggestW + 1, double);

    for (int i = 0; i <= biggestW; i++)
    {
        logGammaArr[i] = lgamma1p(i);
        loglogGammaArr[i] = log(logGammaArr[i]);
    }
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

double underflowSum(double *q)
{
    // If it takes too long, it can be parallelized, but result must be atomic or be merged as an array
    double result = 0;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
            {
                double currentMax = -DBL_MAX;
                int w_bg = (int)MATRIX_AT_PTR(W, b, g);
                double q_bgc = Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
                if (q_bgc == 0)
                {
                    continue;
                }
                // Obtain the maximum value and also save the sums
                double sums[w_bg + 1];

                for (int i = 1; i <= w_bg; i++)
                {
                    double bigsum1 = loglogGammaArr[i] + logGammaArr[w_bg] - logGammaArr[i] - logGammaArr[w_bg - i];
                    bigsum1 += i * log(q_bgc) + (w_bg - i) * log(1 - q_bgc);
                    sums[i] = bigsum1;
                    if (bigsum1 > currentMax)
                    {
                        currentMax = bigsum1;
                    }
                }
                // Shift every value by the sum and exp(max)
                double partialSum = 0;
                for (int i = 1; i <= w_bg; i++)
                {
                    partialSum += exp(sums[i] - currentMax);
                }
                result += exp(currentMax) * partialSum;
            }
        }
    }
    return result;
}
/*
 * Computes the big 'Q' for the log-likelihood.
 * This value needs to be aggregated to the log-likelihood
 */
double computeQ(double *q, Matrix const *probabilities)
{

    double thirdTerm = underflowSum(q);
    double total = -thirdTerm;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            int w_bg = (int)MATRIX_AT_PTR(W, b, g);
            total += logGammaArr[w_bg]; // Second term
            double qsum = 0;
            double firstTerm = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double q_bgc = Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES);
                double p_gc = MATRIX_AT_PTR(probabilities, g, c);
                firstTerm += (p_gc == 0.0 || q_bgc == 0.0) ? 0 : q_bgc * log(MATRIX_AT_PTR(probabilities, g, c));
            }
            // First term
            total += firstTerm * w_bg;
        }
    }
    return total;
}

double logsumexp(double *log_a_shift, OmegaSet *currentSet, int g, int c, double wbg, double *ll)
{
    // Compute log-sum-exp of numerator and denominator to compute q
    int size = currentSet->size;
    // b[i] should be z_bgc/w_bg
    double sum_exp_num = 0.0;
    double sum_exp_den = 0.0;
    double v;
    for (int i = 0; i < size; i++)
    {
        v = exp(log_a_shift[i]);
        sum_exp_num += v;
        Matrix *currentMatrix = currentSet->data[i];

        sum_exp_den += (MATRIX_AT_PTR(currentMatrix, g, c) / wbg) * v;
    }
    double log_sum_exp_num = log(sum_exp_num);

    // update the log-likelihood with the terms of H
    // *ll = 0;
    for (int i = 0; i < size; i++)
    {
        double val = log_a_shift[i] - log_sum_exp_num;
        *ll -= exp(val) * val;
    }

    // return q_{bgc}
    // return exp(log_sum_exp_num - log(sum_exp_den));
    return exp(log(sum_exp_den) - log_sum_exp_num);
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
double *computeQHitAndRun(Matrix const *probabilities, QMethodInput params, double *ll)
{
    // ---- Compute the variables that can be reused ---- //
    if (OMEGASET == NULL)
    {
        generateOmegaSet(params.M, params.S);
    }
    if (multinomialVals == NULL)
    {
        preComputeMultinomial();
    }
    if (logGammaArr == NULL)
    {
        precomputeLogGammas();
    }
    // ---...--- //

    // ---- Compute the final values and fill the returning array ---- //
    double *array2 = (double *)Calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, double); // Array to return
    // ---- Use a static assignment since the workload is even between threads ----

    *ll = 0;
#pragma omp parallel for collapse(2)
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        OmegaSet *currentSet = OMEGASET[b];
        double *multiplicationValues = (currentSet->size <= 10000) ? (double[10000]){1}                // Stack
                                                                   : Calloc(currentSet->size, double); // Heap

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            double W_bg = MATRIX_AT_PTR(W, b, g);
            if (W_bg == 0)
            {
                for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                }
                continue;
            }

            // --- Precompute multiplicationValues for this (b, g) combination ---
            double firstTerm = 0;
            double max = -DBL_MAX;
            for (size_t s = 0; s < currentSet->size; s++)
            { // --- For each sample given a group and a ballot box
                Matrix *currentMatrix = currentSet->data[s];
                double a_i = logarithmicProduct(probabilities, b, s) + multinomialVals[b][s];
                multiplicationValues[s] = a_i;
                firstTerm += multiplicationValues[s];
                max = multiplicationValues[s] > max ? multiplicationValues[s] : max;
            }
            for (size_t s = 0; s < currentSet->size; s++)
            {
                multiplicationValues[s] -= max;
                firstTerm -= max; // Maybe just substract max * currentSet->size
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box
                double secondTerm = 0;
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) =
                    logsumexp(multiplicationValues, currentSet, g, c, MATRIX_AT_PTR(W, b, g), ll);
            }
        }
        if (currentSet->size > 10000)
            Free(multiplicationValues);
    }
    *ll = *ll * (-1);
    *ll += computeQ(array2, probabilities);

    return array2;
}

void cleanHitAndRun()
{
    if (OMEGASET != NULL)
    {
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
        {
            if (OMEGASET[b] != NULL) // Ensure it's valid before freeing
            {
                for (size_t s = 0; s < OMEGASET[b]->size; s++)
                {
                    if (OMEGASET[b]->data[s] != NULL)
                    {
                        freeMatrix(OMEGASET[b]->data[s]); // Free individual matrices
                        OMEGASET[b]->data[s] = NULL;      // Avoid dangling pointers
                    }
                }

                Free(OMEGASET[b]->data); // Free the data array
                OMEGASET[b]->data = NULL;

                Free(OMEGASET[b]); // Free the OmegaSet struct
                OMEGASET[b] = NULL;
            }
        }

        Free(OMEGASET); // Free the OMEGASET array
        OMEGASET = NULL;
    }
    if (multinomialVals != NULL)
    {
        Free(multinomialVals);
        multinomialVals = NULL;
    }
    if (logGammaArr != NULL)
    {
        Free(logGammaArr);
        logGammaArr = NULL;
    }
    if (loglogGammaArr != NULL)
    {
        Free(loglogGammaArr);
        loglogGammaArr = NULL;
    }
}
//__attribute__((destructor)) void cleanEverything()
//{

//   cleanHitAndRun();
//}
