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

int lessThanColRow(Matrix mat, int b, int g, int c, int candidateVotes, int groupVotes)
{
    int groupSum = 0;
    int canSum = 0;
    for (uint16_t i = 0; i < TOTAL_GROUPS; i++)
    {
        canSum += MATRIX_AT(mat, i, c);
    }
    for (uint16_t j = 0; j < TOTAL_CANDIDATES; j++)
    {
        groupSum += MATRIX_AT(mat, g, j);
    }
    int slackC = candidateVotes - canSum;
    int slackG = groupVotes - groupSum;

    return MIN(slackC, slackG);
}

Matrix startingPoint3(int b)
{
    // ---- Retrieve the initial variables ---- //
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    double *groupVotes = getRow(W, b);
    double *candidateVotes = getColumn(X, b);

    // ---- Calculate the expected value ---- //
    double totalC = 0;
    double totalG = 0;
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double mult = groupVotes[g] * candidateVotes[c];

            if (g == 0)
                totalC += candidateVotes[c];

            // In case of mismatch, we divide for the maximum

            MATRIX_AT(toReturn, g, c) = mult;
        }
        totalG += groupVotes[g];
    }
    // ---...--- //

    // ---- Division for mismatchs ---- //
    double divide = MAX(BALLOTS_VOTES[b], totalC);
    divide = MAX(divide, totalG);

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            double newValue = MATRIX_AT(toReturn, g, c) / divide;
            double floored = floor(newValue);
            MATRIX_AT(toReturn, g, c) = floored;
        }
    }

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            int groupRestriction = groupVotes[g];
            int candidateRestriction = candidateVotes[c];

            int m = lessThanColRow(toReturn, b, g, c, candidateRestriction, groupRestriction);
            if (m > 0)
            {
                MATRIX_AT(toReturn, g, c) += m;
            }
        }
    }
    // ---...--- //

    // ---...--- //
    Free(groupVotes);
    Free(candidateVotes);
    return toReturn;
}

// 0 y C*(C-1), G(G-1)-1, o sea, sin los divididos por 2, la razón, es porque creo que es más sencillo encodear y
// decodear de la segunda manera
void allocateRandoms(int M, int S, uint8_t **c1, uint8_t **c2, uint8_t **g1, uint8_t **g2)
{
    uint32_t size = M * S;
    // Allocate memory correctly
    *c1 = (uint8_t *)Calloc(size, uint8_t);
    *c2 = (uint8_t *)Calloc(size, uint8_t);
    *g1 = (uint8_t *)Calloc(size, uint8_t);
    *g2 = (uint8_t *)Calloc(size, uint8_t);

    GetRNGstate(); // Ensure R's RNG is properly initialized
                   // Fill arrays with random indices
    int allow_repeat = (TOTAL_CANDIDATES <= 1 || TOTAL_GROUPS <= 1);

    for (int i = 0; i < size; i++)
    {
        if (i % 400 == 0) // Checks condition every 400 iterations
            R_CheckUserInterrupt();
        (*c1)[i] = (uint8_t)(unif_rand() * TOTAL_CANDIDATES);
        (*g1)[i] = (uint8_t)(unif_rand() * TOTAL_GROUPS);
        do
        {
            (*c2)[i] = (uint8_t)(unif_rand() * TOTAL_CANDIDATES);
            (*g2)[i] = (uint8_t)(unif_rand() * TOTAL_GROUPS);
        } while (!allow_repeat && ((*c2)[i] == (*c1)[i] || (*g2)[i] == (*g1)[i]));
    }
    PutRNGstate(); // Finalize RNG state to prevent repeatability
}
/*
 * @brief Precomputes the sets used for the simulation.
 *
 * Precomputes the sets that are independent from each EM iteration. It is made with parallelism (NOT SUPPORTED)
 * towards the ballot boxes and with a static assignment for ensuring reproducibility.
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
    uint8_t *c1 = NULL;
    uint8_t *c2 = NULL;
    uint8_t *g1 = NULL;
    uint8_t *g2 = NULL;

    uint32_t arraySize = M * S;

    allocateRandoms(M, S, &c1, &c2, &g1, &g2);
    // Compute the partition size
    int partitionSize = M / TOTAL_BALLOTS;
    if (partitionSize == 0)
        partitionSize = 1; // Prevent division by zero in extreme cases

// ---- Perform the main iterations ---- //
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
        Matrix startingZ = startingPoint3(b);
        if (b % 5 == 0)
        {
            int totalVals = 0;
            for (int g = 0; g < TOTAL_GROUPS; g++)
            {
                for (int c = 0; c < TOTAL_CANDIDATES; c++)
                {
                    totalVals += MATRIX_AT(startingZ, g, c);
                }
            }
        }
        int ballotShift = floor(((double)b / TOTAL_BALLOTS) * (M * S));

        for (int s = 0; s < S; s++)
        {                                // --- For each sample given a ballot box
            if (S > 5000 && S % 50 == 0) // If there's a big amount of samples, check interrupts
                R_CheckUserInterrupt();
            // TODO: El sampling debe hacerse de tamaño M*S
            // ---- Copy the initial matrix ----
            Matrix steppingZ = copyMatrix(&startingZ);
            for (int m = 0; m < M; m++)
            { // --- For each step size given a sample and a ballot box
                // ---- Sample random indexes ---- //
                int shiftIndex = (s * M + ballotShift + m) % (M * S);
                uint8_t randomCDraw = c1[shiftIndex];
                uint8_t randomCDraw2 = c2[shiftIndex];
                uint8_t randomGDraw = g1[shiftIndex];
                uint8_t randomGDraw2 = g2[shiftIndex];

                // decode(randomCDraw, TOTAL_CANDIDATES, &c1, &c2);
                //  decode(randomGDraw, TOTAL_GROUPS, &g1, &g2);

                // ---- Check non negativity condition ---- //
                double firstSubstraction = MATRIX_AT(steppingZ, randomGDraw, randomCDraw);
                double secondSubstraction = MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2);

                if (firstSubstraction <= 0 || secondSubstraction <= 0)
                    continue;
                // ---...--- //

                // ---- Asign changes on the new matrix ---- //
                MATRIX_AT(steppingZ, randomGDraw, randomCDraw) -= 1;
                MATRIX_AT(steppingZ, randomGDraw2, randomCDraw2) -= 1;
                MATRIX_AT(steppingZ, randomGDraw, randomCDraw2) += 1;
                MATRIX_AT(steppingZ, randomGDraw2, randomCDraw) += 1;
                //  ---...--- //
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
    Free(c1);
    Free(c2);
    Free(g1);
    Free(g2);
}

void encode()
{
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        OmegaSet *set = OMEGASET[b];
        int S = set->size;

        MatrixHash *hashTable = NULL;

        // First pass: identify and count unique matrices
        for (int s = 0; s < S; s++)
        {
            Matrix *current = set->data[s];
            unsigned int key = computeMatrixKey(current);

            MatrixHash *entry = NULL;
            HASH_FIND(hh, hashTable, &key, sizeof(unsigned int), entry);

            int found = 0;
            for (MatrixHash *iter = entry; iter != NULL; iter = iter->hh.next)
            {
                if (iter->key == key && matricesAreEqual(iter->matrix, current))
                {
                    iter->count++;
                    found = 1;
                    break;
                }
            }

            if (!found)
            {
                MatrixHash *newEntry = Calloc(1, MatrixHash);
                newEntry->key = key;
                newEntry->matrix = current;
                newEntry->count = 1;
                HASH_ADD(hh, hashTable, key, sizeof(unsigned int), newEntry);
            }
        }

        // Allocate counts and assign total count per sample
        if (set->counts == NULL)
        {
            set->counts = Calloc(S, int);
        }

        for (int s = 0; s < S; s++)
        {
            Matrix *current = set->data[s];
            unsigned int key = computeMatrixKey(current);

            MatrixHash *entry = NULL;
            HASH_FIND(hh, hashTable, &key, sizeof(unsigned int), entry);

            for (MatrixHash *iter = entry; iter != NULL; iter = iter->hh.next)
            {
                if (iter->key == key && matricesAreEqual(iter->matrix, current))
                {
                    set->counts[s] = iter->count; // All duplicates now get the same count
                    break;
                }
            }
            if (set->counts[s] == 0) // Hash collision case
                set->counts[s] = 1;
        }

        // Clean up
        MatrixHash *entry, *tmp;
        HASH_ITER(hh, hashTable, entry, tmp)
        {
            HASH_DEL(hashTable, entry);
            Free(entry);
        }
    }
}
/**
 * @brief Computes the pre-computable values of the expression that doesn't depend on EM iterations
 *
 * Given a ballot box index and a matrix represent an element from the Hit and Run OmegaSet, it computes the
 * following:
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
    Matrix *currentMatrix = OMEGASET[b]->data[setIndex]; // Z
    // ---...--- //
    // ---- Main computation ---- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    { // ---- For each candidate
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group
            // Cambiar a log(1) if probabilidad 0
            // Producto punto
            double prob = MATRIX_AT_PTR(probabilities, g, c);
            log_result += (prob > 0.0) ? MATRIX_AT_PTR(currentMatrix, g, c) * log(prob) : 0.0;
        }
    }
    // --- ... --- //
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
        loglogGammaArr[i] = logGammaArr[i] != 0 ? log(logGammaArr[i]) : 1;
    }
}

/**
 * @brief Precomputes the multinomial multiplication that is independent for each EM iteration.
 *
 * Calls the main function for computing all of the calculations related with the final result that are independent
 * from each EM call. Specifically, for each ballot box and its simulations, the following is calculated:
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
    multinomialVals = Calloc(TOTAL_BALLOTS, double *);
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
                int w_bg = (int)MATRIX_AT_PTR(W, b, g);
                double currentMax = -DBL_MAX;
                double q_bgc = Q_3D(q, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES);
                if (q_bgc <= 0 || w_bg <= 1)
                {
                    continue;
                }
                // Obtain the maximum value and also save the sums
                double sums[w_bg + 1];

                for (int i = 2; i <= w_bg; i++)
                {
                    double bigsum1 = loglogGammaArr[i] + logGammaArr[w_bg] - logGammaArr[i] - logGammaArr[w_bg - i];
                    if (q_bgc != 1 || w_bg != i)
                    {
                        bigsum1 += i * log(q_bgc) + (w_bg - i) * log(1 - q_bgc);
                    }
                    sums[i] = bigsum1;
                    if (bigsum1 > currentMax)
                    {
                        currentMax = bigsum1;
                    }
                }
                // Shift every value by the sum and exp(max)
                double partialSum = 0;
                for (int i = 2; i <= w_bg; i++)
                {
                    // partialSum += exp(sums[i] - currentMax);
                    partialSum += exp(sums[i]);
                }
                // result += exp(currentMax) * partialSum;
                result += partialSum;
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
    double borrar = 0;
    double borrar2 = 0;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            int w_bg = (int)MATRIX_AT_PTR(W, b, g);
            borrar += w_bg == 0 ? 0 : logGammaArr[w_bg];
            total += w_bg == 0 ? 0 : logGammaArr[w_bg]; // Second term
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
            borrar2 += firstTerm * w_bg;
        }
    }
    // Rprintf("----------\nQ: %f + %f + %f\n", borrar2, borrar, -thirdTerm);

    return total;
}

double *computeQHitAndRun(Matrix const *probabilities, QMethodInput params, double *ll)
{
    // ---- Compute the variables that can be reused ---- //
    if (OMEGASET == NULL)
    {
        generateOmegaSet(params.M, params.S);
        encode();
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
            double max = -DBL_MAX;
            for (size_t s = 0; s < currentSet->size; s++)
            {                                                // --- For each sample given a group and a ballot box
                Matrix *currentMatrix = currentSet->data[s]; // Z
                double a_i = logarithmicProduct(probabilities, b, s) +
                             multinomialVals[b][s]; // $\sum_{g\in G}\sum_{c\in C}z_{bgc}*\log(p_{gc})+$
                multiplicationValues[s] = a_i;

                max = multiplicationValues[s] > max ? multiplicationValues[s] : max;
            }
            // ---- Shift the values by the maximum
            for (size_t s = 0; s < currentSet->size; s++)
            {
                multiplicationValues[s] -= max;
            }
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box
                // ---- Obtain the summatory over all of the values ---- //
                double sum_exp_den = 0.0;
                double sum_exp_num = 0.0;
                double v = 0;
                for (int i = 0; i < currentSet->size; i++)
                {                                     // --- For each sample
                    v = exp(multiplicationValues[i]); // exp(ls - m)
                    sum_exp_num += v;
                    Matrix *currentMatrix = currentSet->data[i];
                    sum_exp_den += (MATRIX_AT_PTR(currentMatrix, g, c) / W_bg) * v;
                }
                double log_sum_exp_num = log(sum_exp_num);
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) =
                    exp(log(sum_exp_den) - log_sum_exp_num);
                // ---...--- //
            } // --- End candidate loop
        } // --- End group loop

        // ---- Calculate the log-likelihood ---- //
        double sum_exp_num = 0.0;
        double v = 0;
        for (int i = 0; i < currentSet->size; i++)
        { // --- For each sample
            v = exp(multiplicationValues[i]);
            sum_exp_num += v;
        }
        for (int i = 0; i < currentSet->size; i++)
        { // --- For each sample
            double val = exp(multiplicationValues[i]) / sum_exp_num;
            // Rprintf("S = %d B = %d H = %f R = %d\n", i, b, val, currentSet->counts[i]);
            *ll -= val * log(currentSet->counts[i] * val);
        }
        // ---...--- //

        // ---- Free allocated memory ---- //
        if (currentSet->size > 10000)
            Free(multiplicationValues);
        // ---...--- //
    } // --- End ballot box loop

    // Calculo Q
    double toprint = computeQ(array2, probabilities);
    // Rprintf("- Valor de Q = %f\n- Valor de H = %f\n----------\n", toprint, *ll);

    // *ll += computeQ(array2, probabilities);
    *ll += toprint;

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

            if (OMEGASET[b]->counts != NULL)
            {
                Free(OMEGASET[b]->counts);
                OMEGASET[b]->counts = NULL;
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
