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

#include "multivariate-cdf.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <R_ext/Random.h>
#include <Rmath.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/*
 * @brief Compute the Montecarlo approximation proposed by Alan genz towards the most recent method, using univariate
 * conditional with quasirandom numbers.
 *
 * Computes an heuristic from the first Multivariate CDF proposed. More details at the following paper:
 * https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities
 *
 * Note that this method gives the option to impose an error threshold. So, it stops iterating if, either the maximum
 * iterations are accomplished or the error threshold is passed.
 *
 * @param[in] *cholesky. The cholesky matrix of the current iteration.
 * @param[in] *lowerBounds. An array with the initial lower bounds.
 * @param[in] *upperBounds. An array with the initial upper bounds.
 * @param[in] epsilon. The minimum error threshold accepted.
 * @param[in] iterations. The maximum iterations
 * @param[in] mvnDim. The dimension of the multivariate normal.
 *
 * @return The value of the estimated integral
 */
double genzMontecarloNew(const Matrix *cholesky, const double *lowerBounds, const double *upperBounds, double epsilon,
                         int iterations, int mvnDim)
{
    GetRNGstate();

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    double mean = 0;
    int currentIterations = 1;
    double a[mvnDim], b[mvnDim], y[mvnDim], currentError;
    a[0] = pnorm(lowerBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 0.0, 1.0, 1, 0);
    b[0] = pnorm(upperBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 0.0, 1.0, 1, 0);
    // ---...--- //

    do
    {
        // ---- Generate a pseudoRandom number for each iteration
        double pseudoRandom[mvnDim];
        for (int i = 0; i < mvnDim; i++)
        {
            pseudoRandom[i] = unif_rand();
        }

        // ---- Compute the base case
        y[0] = qnorm(a[0] + pseudoRandom[0] * (b[0] - a[0]), 0.0, 1.0, 1, 0);
        double summatory;
        double P = b[0] - a[0];

        // ---- Do the main loop ---- //
        for (int i = 1; i < mvnDim; i++)
        {
            // ---- Note that the summatory is equivalent to $\sum_{j=1}^{i-1}c_{ij}*y_{j}$.
            summatory = 0;
            for (int j = 0; j < i; j++)
            {
                summatory += MATRIX_AT_PTR(cholesky, i, j) * y[j];
            }
            a[i] = pnorm((lowerBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 0.0, 1.0, 1, 0);
            b[i] = pnorm((upperBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 0.0, 1.0, 1, 0);
            double difference = b[i] - a[i];
            y[i] = qnorm(a[i] + pseudoRandom[i] * difference, 0.0, 1.0, 1, 0);
            P *= difference;
        }
        // ---...--- //
        // ---- Get the stopping parameters ---- //
        intsum += P;
        currentIterations += 1;
        mean = intsum / currentIterations;
        varsum += pow(P - mean, 2);
        currentError = sqrt(varsum / (currentIterations * (currentIterations - 1)));
        // ---...--- //
    } while (currentError > epsilon && currentIterations < iterations);

    PutRNGstate();
    return mean;
}

// All of the conventions used are from genz paper
/*
 * @brief Compute the Montecarlo approximation proposed by Alan genz
 *
 * Compute a `manual` Montecarlo simulation, proposed by Dr. Alan genz book "Computation of Multivariate Normal and t
 * Probabilities". All of the variable names and conventions are adopted towards his proposed algorithm aswell. Refer to
 * https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities for more
 * information.
 *
 * Note that this method gives the option to impose an error threshold. So, it stops iterating if, either the maximum
 * iterations are accomplished or the error threshold is passed.
 *
 * @param[in] *cholesky. The cholesky matrix of the current iteration.
 * @param[in] *lowerBounds. An array with the initial lower bounds.
 * @param[in] *upperBounds. An array with the initial upper bounds.
 * @param[in] epsilon. The minimum error threshold accepted.
 * @param[in] iterations. The maximum iterations
 * @param[in] mvnDim. The dimension of the multivariate normal.
 *
 * @return The value of the estimated integral
 */
double genzMontecarlo(const Matrix *cholesky, const double *lowerBounds, const double *upperBounds, double epsilon,
                      int iterations, int mvnDim)
{

    // ---- Initialize randomizer ---- //
    GetRNGstate();

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    int currentIterations = 0;
    double currentError = DBL_MAX;
    double d[mvnDim], e[mvnDim], f[mvnDim];
    double diag0 = MATRIX_AT_PTR(cholesky, 0, 0);
    if (diag0 == 0)
        return 0;
    d[0] = pnorm(lowerBounds[0] / diag0, 0.0, 1.0, 1, 0);
    e[0] = pnorm(upperBounds[0] / diag0, 0.0, 1.0, 1, 0);
    f[0] = (e[0] - d[0]);
    //  ---...--- //
    do
    {
        // ---- Initialize the loop variables ---- //
        d[0] = pnorm(lowerBounds[0] / diag0, 0.0, 1.0, 1, 0);
        e[0] = pnorm(upperBounds[0] / diag0, 0.0, 1.0, 1, 0);
        f[0] = (e[0] - d[0]);

        double randomVector[mvnDim - 1];
        for (int i = 0; i < mvnDim - 1; i++)
        // ---...--- //

        // ---- Generate the random vector ---- //
        { // --- For each dimension
            // ---- Generate random values in [0,1) ----
            randomVector[i] = unif_rand();
        }
        // ---...--- //

        // ---- Calculate the integral with their new bounds ---- //
        double summatory;
        double y[mvnDim];
        for (int i = 1; i < mvnDim; i++)
        {
            double draw = d[i - 1] + randomVector[i - 1] * (e[i - 1] - d[i - 1]);
            y[i - 1] = qnorm(draw, 0.0, 1.0, 1, 0);
            // ---- Note that the summatory is equivalent to $\sum_{j=1}^{i-1}c_{ij}*y_{j}$.

            summatory = 0;
            for (int j = 0; j < i; j++)
            {
                summatory += MATRIX_AT_PTR(cholesky, i, j) * y[j];
            }

            double diagii = MATRIX_AT_PTR(cholesky, i, i);
            d[i] = diagii != 0 ? pnorm((lowerBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0) : 0;
            e[i] = diagii != 0 ? pnorm((upperBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0) : 0;
            f[i] = (e[i] - d[i]) * f[i - 1];
        }
        // ---...--- //

        // ---- Compute the final indicators from the current loop ---- //
        intsum += f[mvnDim - 1];
        varsum += pow(f[mvnDim - 1], 2);
        currentIterations += 1;
        if (currentIterations < 2)
            continue;
        double mean = intsum / currentIterations;
        double variance = (varsum / currentIterations) - pow(mean, 2);
        // ---...--- //
        currentError = sqrt(variance / currentIterations);

    } while (currentError > epsilon && currentIterations < iterations);
    // ---...--- //
    PutRNGstate();

    return intsum / currentIterations;
}

/*
 * @brief Calls the main function to start all of the Montecarlo process
 *
 * Calls the function with all of the parameters needed to get the Montecarlo simulation of the Multivariate CDF.
 *
 * @param[in] *chol A matrix with the cholenksy values of the current group
 * @param[in] *mu An array with the average values of the current feature vector
 * @param[in] *lowerLimits An array with the lower bounds of the integral (defined by the hypercube)
 * @param[in] *upperLimits An array with the upper bounds of the integral (defined by the hypercube)
 * @param[in] mvnDim The dimensions of the multivariate normal. Usually it's C-1.
 * @param[in] maxSamples Amount of samples for the Montecarlo simulation.
 * @param[in] epsilon The error threshold used for the genz Montecarlo.
 * @param[in] *method The method for calculating the Montecarlo simulation. Currently available methods are `Plain`,
 * `Miser` and `Vegas`.
 *
 * @note Refer to https://www.gnu.org/software/gsl/doc/html/montecarlo.html
 *
 * @return The result of the approximated integral
 */

double Montecarlo(Matrix *chol, double *mu, double *lowerLimits, double *upperLimits, int mvnDim, int maxSamples,
                  double epsilon, const char *method)
{
    // ---- Set up the initial parameter ---- //
    double result;
    // ---...--- //
    // ---- Case where there are no values ---- //
    if (MATRIX_AT_PTR(chol, 0, 0) == 0 ||
        memcmp(lowerLimits, upperLimits, sizeof(double) * (TOTAL_CANDIDATES - 1)) == 0)
        return 0;

    // ---- Perform integration ---- //
    if (strcmp(method, "genz") == 0)
    {
        result = genzMontecarlo(chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else if (strcmp(method, "genz2") == 0)
    {
        result = genzMontecarloNew(chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else
    {
        error("Multivariate CDF: An invalid method was handed to the Montecarlo simulation for calculating the "
              "Multivariate CDF "
              "integral.\nThe method handed is:\t%s\nThe current available methods are `genz` or `genz2`"
              ".\n",
              method);
    }
    // ---...--- //
}

/*
 * @brief Gets the a matrix of mu values and the inverse sigma for a given ballot box
 *
 * @param[in] b The ballot box index.
 * @param[in] *probabilitiesReduced Probabilities matrix without the last candidate.
 * @param[in, out] **cholesky An array of `g` matrices that stores (`c-1`X`c-1`) matrices with the cholesky values. The
 * cholesky matrix is filled on the upper and lower triangle with the same values since it's simmetric.
 * @param[in, out] *mu A matrix of size (`g`x`c-1`) with the averages per group.
 *
 * @return void. Results to be written on cholesky and mu
 */
void getMainParameters(int b, Matrix const probabilitiesReduced, Matrix **cholesky, Matrix *mu)
{

    // --- Initialize empty array --- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        cholesky[g] = (Matrix *)Calloc(1, Matrix);
        *cholesky[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }
    // ---...--- //

    // ---- Get mu and sigma ---- //
    getAverageConditional(b, &probabilitiesReduced, mu, cholesky);

    // ---...--- //
}

/*
 * @brief Computes the `q` values from the multivariate CDF method
 *
 * @param[in] *probabilities A pointer towards the current probabilities matrix.
 * @param[in] monteCarloSamples The amount of samples to use in the Monte Carlo simulation
 * @param[in] epsilon The error threshold used for the genz Montecarlo method.
 * @param[in] *method The method for calculating the Montecarlo simulation. Currently available methods are `Plain`,
 * `Miser` and `Vegas`.
 *
 * @return A contiguos array with all the new probabilities
 */
double *computeQMultivariateCDF(Matrix const *probabilities, QMethodInput params, double *ll)
{
    *ll = 0.0;
    int monteCarloSamples = params.monteCarloIter;
    double epsilon = params.errorThreshold;
    const char *method = params.simulationMethod;

    // ---- Define initial variables ---- //
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 = (double *)Calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, double); // Array to return
    double *logArray = (double *)Calloc(TOTAL_BALLOTS, double);
    // --- ... --- //

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        // ---- Get the values of the Multivariate CDF that only depends on `b` ---- //
        // ---- Mu and inverse Sigma matrix ----
        Matrix mu = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
        Matrix **choleskyVals = (Matrix **)Calloc(TOTAL_GROUPS, Matrix *);
        getMainParameters(b, probabilitiesReduced, choleskyVals, &mu);
        // ---- Array with the results of the Xth candidate on ballot B ----
        double *feature = getColumn(X, b); // Of size C-1
                                           // ---...--- //

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            // ---- Define the current values to use that only depends on `g` ---- //
            Matrix *currentCholesky = choleskyVals[g];

            double *currentMu = getRow(&mu, g);
            // ---- Initialize empty variables to be filled ----
            double montecarloResults[TOTAL_CANDIDATES];
            double denominator = 0;
            // ---...--- //

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box

                // ---- Define the borders of the hypercube ---- //
                // ---- First, make a copy of the feature vector ----
                double *featureCopyA = (double *)Calloc((TOTAL_CANDIDATES - 1), double);
                double *featureCopyB = (double *)Calloc((TOTAL_CANDIDATES - 1), double);

                memcpy(featureCopyA, feature, (TOTAL_CANDIDATES - 1) * sizeof(double));
                memcpy(featureCopyB, feature, (TOTAL_CANDIDATES - 1) * sizeof(double));

                // ---- Substract/add and normalize ----
                // ---- Note that the bounds NEEDS to be standarized ----
                // ---- $$a=\sigma^{-1}(a-\mu)$$ ----
                for (uint16_t k = 0; k < TOTAL_CANDIDATES - 1; k++)
                { // --- For each candidate coordinate that is going to be integrated
                    featureCopyA[k] -= 0.5;
                    featureCopyB[k] += 0.5;
                    if (k == c)
                    {
                        featureCopyA[k] -= 1.0;
                        featureCopyB[k] -= 1.0;
                    }
                    featureCopyA[k] -= currentMu[k];
                    featureCopyB[k] -= currentMu[k];
                }
                // ---...--- //
                // ---- Save the results and add them to the denominator ---- //
                if (TOTAL_CANDIDATES != 2)
                    montecarloResults[c] = Montecarlo(currentCholesky, currentMu, featureCopyA, featureCopyB,
                                                      (int)TOTAL_CANDIDATES - 1, monteCarloSamples, epsilon, method);

                else
                {
                    montecarloResults[c] =
                        (pnorm(featureCopyB[0], 0.0, sqrt(MATRIX_AT_PTR(currentCholesky, 0, 0)), 1, 0) -
                         pnorm(featureCopyA[0], 0.0, sqrt(MATRIX_AT_PTR(currentCholesky, 0, 0)), 1, 0));
                }
                montecarloResults[c] *= MATRIX_AT_PTR(probabilities, g, c);
                denominator += !isnan(montecarloResults[c]) ? montecarloResults[c] : 0;
                logArray[b] += g == 0 && !isnan(montecarloResults[c]) ? montecarloResults[c] : 0;
                // TODO: Make an arena for this loop
                Free(featureCopyA);
                Free(featureCopyB);
                // ---...--- //
            } // --- End c loop
            Free(currentMu);
            freeMatrix(currentCholesky);
            freeMatrix(choleskyVals[g]);
            Free(choleskyVals[g]);

            // ---- Add the final results to the array ----//
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // --- For each candidate
                double result = montecarloResults[c] / denominator;
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) =
                    !isnan(result) && !isinf(result) ? result : 0;
            } // --- End c loop
            // ---...--- //
        } // --- End g loop
        Free(feature);
        freeMatrix(&mu);
        Free(choleskyVals);
    } // --- End b loop
    freeMatrix(&probabilitiesReduced);

    // Compute the log-likelihood
    double finalLikelihood = 0;
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        finalLikelihood += logArray[b] != 0 ? log(logArray[b]) : 0;
    }
    *ll = finalLikelihood;
    // *ll = log(*ll);
    free(logArray);
    return array2;
}
