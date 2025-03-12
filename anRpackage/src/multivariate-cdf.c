#include "multivariate-cdf.h"
#include "R_ext/Random.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// All of the conventions used are from genz paper
/*
 * @brief Compute the Montecarlo approximation proposed by Alan Genz towards the most recent method, using univariate
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
    double error = epsilon + 1;
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

// Might implement later
double genzMontecarloChunk(const Matrix *cholesky, const double *lowerBounds, const double *upperBounds, double epsilon,
                           int iterations, int mvnDim)
{

    // ---- Initialize randomizer ---- //
    GetRNGstate();

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    int currentIterations = 0;
    double currentError = 1000;
    double d[mvnDim], e[mvnDim], f[mvnDim];
    double diag0 = MATRIX_AT_PTR(cholesky, 0, 0);
    d[0] = pnorm(lowerBounds[0] / diag0, 0.0, 1.0, 1, 0);
    e[0] = pnorm(upperBounds[0] / diag0, 0.0, 1.0, 1, 0);
    f[0] = (e[0] - d[0]);
    //  ---...--- //

    // ---- Compute the main while loop ---- //
    for (int k = 0; k < 5; k++)
    {
#pragma omp parallel for schedule(static)
        for (int j = 0; j < 250; j++)
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
                d[i] = pnorm((lowerBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0);
                e[i] = pnorm((upperBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0);
                f[i] = (e[i] - d[i]) * f[i - 1];
            }
            // ---...--- //

            // ---- Compute the final indicators from the current loop ---- //
#pragma omp atomic
            intsum += f[mvnDim - 1];
#pragma omp atomic
            varsum += pow(f[mvnDim - 1], 2);
        }
        currentIterations += 250;
        double mean = intsum / currentIterations;
        double variance = (varsum / currentIterations) - pow(mean, 2);
        currentError = sqrt(variance / currentIterations);
        // ---...--- //
        if (currentError < epsilon)
            break;
    }
    // ---...--- //
    PutRNGstate();

    return intsum / currentIterations;
}

// All of the conventions used are from genz paper
/*
 * @brief Compute the Montecarlo approximation proposed by Alan Genz
 *
 * Compute a `manual` Montecarlo simulation, proposed by Dr. Alan Genz book "Computation of Multivariate Normal and t
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
    double currentError = 1000;
    double d[mvnDim], e[mvnDim], f[mvnDim];
    double diag0 = MATRIX_AT_PTR(cholesky, 0, 0);
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
            d[i] = pnorm((lowerBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0);
            e[i] = pnorm((upperBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0);
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
 * @param[in] epsilon The error threshold used for the Genz Montecarlo.
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
    // ---- Set up the initial parameters ---- //
    double result, error;
    // ---...--- //

    // ---- Perform integration ---- //
    if (strcmp(method, "Genz") == 0)
    {
        result = genzMontecarloChunk(chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else if (strcmp(method, "Genz2") == 0)
    {
        result = genzMontecarloNew(chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else
    {
        fprintf(stderr,
                "\nAn invalid method was handed to the Montecarlo simulation for calculating the Multivariate CDF "
                "integral.\nThe method handed is:\t%s\nThe current available methods are `Genz` or `Genz2`"
                ".\n",
                method);
        exit(EXIT_FAILURE);
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
        cholesky[g] = (Matrix *)malloc(sizeof(Matrix));
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
 * @param[in] epsilon The error threshold used for the Genz Montecarlo method.
 * @param[in] *method The method for calculating the Montecarlo simulation. Currently available methods are `Plain`,
 * `Miser` and `Vegas`.
 *
 * @return A contiguos array with all the new probabilities
 */
double *computeQMultivariateCDF(Matrix const *probabilities, QMethodInput params)
{
    int monteCarloSamples = params.monteCarloIter;
    double epsilon = params.errorThreshold;
    const char *method = params.simulationMethod;

    // ---- Define initial variables ---- //
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
// --- ... --- //
#pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        // ---- Get the values of the Multivariate CDF that only depends on `b` ---- //
        // ---- Mu and inverse Sigma matrix ----
        Matrix mu = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
        Matrix **choleskyVals = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));
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
                double *featureCopyA = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
                double *featureCopyB = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));

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
                                                      (int)TOTAL_CANDIDATES - 1, monteCarloSamples, epsilon, method) *
                                           MATRIX_AT_PTR(probabilities, g, c);

                else
                {
                    montecarloResults[c] = (pnorm(featureCopyA[0], 0.0, MATRIX_AT_PTR(currentCholesky, 0, 0), 1, 0) -
                                            pnorm(featureCopyB[0], 0.0, MATRIX_AT_PTR(currentCholesky, 0, 0), 1, 0)) *
                                           MATRIX_AT_PTR(probabilities, g, c);
                }
                denominator += montecarloResults[c];
                // TODO: Make an arena for this loop
                free(featureCopyA);
                free(featureCopyB);
                // ---...--- //
            } // --- End c loop
            free(currentMu);
            freeMatrix(currentCholesky);
            freeMatrix(choleskyVals[g]);
            free(choleskyVals[g]);

            // ---- Add the final results to the array ----//
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {                         // --- For each candidate
                if (denominator == 0) // Edge case
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = montecarloResults[c] / denominator;
            } // --- End c loop
            // ---...--- //
        } // --- End g loop
        free(feature);
        freeMatrix(&mu);
        free(choleskyVals);
    } // --- End b loop
    freeMatrix(&probabilitiesReduced);
    return array2;
}
