#include "multivariate-cdf.h"
#include "globals.h"
#include "utils/matrixUtils.h"
#include <cblas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h> // Random numbers
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct
{
    Matrix *chol; // Cholesky decomposition
    double *mu;   // The mu array
} IntegrationParams;

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
    gsl_qrng *q = gsl_qrng_alloc(gsl_qrng_sobol, mvnDim);

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    double mean = 0;
    int currentIterations = 0;
    double a[mvnDim], b[mvnDim], y[mvnDim], currentError;
    a[0] = gsl_cdf_gaussian_P(lowerBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 1);
    b[0] = gsl_cdf_gaussian_P(upperBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 1);
    // ---...--- //

    do
    {
        // ---- Generate a pseudoRandom number for each iteration
        double pseudoRandom[mvnDim];
        gsl_qrng_get(q, pseudoRandom);

        // ---- Compute the base case
        y[0] = gsl_cdf_gaussian_Pinv(a[0] + pseudoRandom[0] * (b[0] - a[0]), 1);
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
            a[i] = gsl_cdf_gaussian_P((lowerBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 1);
            b[i] = gsl_cdf_gaussian_P((upperBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 1);
            double difference = b[i] - a[i];

            y[i] = gsl_cdf_gaussian_Pinv(a[i] + pseudoRandom[i] * (difference), 1);
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

    gsl_qrng_free(q);
    return mean;
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
    const gsl_rng_type *T;
    gsl_rng *rng;

    // ---- Set up the GSL RNG environment ----
    gsl_rng_env_setup();
    // ---- Use the default random number generator type
    T = gsl_rng_default;
    // ---- Alocate the random number generator
    rng = gsl_rng_alloc(T);
    // ---...--- //

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    int currentIterations = 0;
    double currentError;
    double d[mvnDim], e[mvnDim], f[mvnDim];
    d[0] = gsl_cdf_gaussian_P(lowerBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 1);
    e[0] = gsl_cdf_gaussian_P(upperBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 1);
    f[0] = e[0] - d[0];
    // ---...--- //

    // ---- Compute the main while loop ---- //
    do
    {
        // ---- Initialize the loop variables ---- //
        double y[mvnDim];
        double randomVector[mvnDim - 1];
        for (int i = 0; i < mvnDim - 1; i++)
        // ---...--- //

        // ---- Generate the random vector ---- //
        { // --- For each dimension
            // ---- Generate random values in [0,1) ----
            randomVector[i] = gsl_rng_uniform(rng);
        }
        // ---...--- //

        // ---- Calculate the integral with their new bounds ---- //
        double summatory;
        for (int i = 1; i < mvnDim; i++)
        {
            y[i - 1] = gsl_cdf_gaussian_Pinv(d[i - 1] + randomVector[i - 1] * (e[i - 1] - d[i - 1]), 1);
            // ---- Note that the summatory is equivalent to $\sum_{j=1}^{i-1}c_{ij}*y_{j}$.
            // summatory += MATRIX_AT_PTR(cholesky, i - 1, i - 1) * y[i - 1];
            summatory = 0;
            for (int j = 0; j < i; j++)
            {
                summatory += MATRIX_AT_PTR(cholesky, i, j) * y[j];
            }
            d[i] = gsl_cdf_gaussian_P((lowerBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 1);
            e[i] = gsl_cdf_gaussian_P((upperBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 1);
            f[i] = (e[i] - d[i]) * f[i - 1];
        }
        // ---...--- //

        // ---- Compute the final indicators from the current loop ---- //
        intsum += f[mvnDim - 1];
        varsum += pow(f[mvnDim - 1], 2);
        currentIterations += 1;
        currentError = sqrt((varsum / (currentIterations - pow(intsum / currentIterations, 2))) / currentIterations);
        // ---...--- //

    } while (currentError > epsilon && currentIterations < iterations);
    // ---...--- //

    gsl_rng_free(rng);

    return intsum / currentIterations;
}

/*
 * @brief Evaluates a given x vector towards the PDF function
 *
 * @param[in] *x A vector with the values to be evaluated. Must be of the same dimension of the multivariate.
 * @param[in] dim The dimension of the multivariate variable
 * @param[in] *params Pointer to the structure with the parameters used to integrate
 *
 * @note Refer to https://www.gnu.org/software/gsl/doc/html/montecarlo.html
 *
 * @return The value of the function. Note that it is a scalar
 */
double integral(double *x, size_t dim, void *params)
{
    // ---- Unpack the integration parameters ---- //
    IntegrationParams *p = (IntegrationParams *)params;
    Matrix *chol = p->chol;
    double *currentMu = p->mu;
    // ---...--- //

    // ---- Obtain the mahanalobis of the multivariate ---- //
    double maha[dim];
    getMahanalobisDist(x, currentMu, chol, maha, dim, true);

    // ---- To get the scalar mahanalobis we need to sum over all contributions ----
    double totalMahanobis = 0;
    for (uint16_t c = 0; c < dim; c++)
    { // --- For each dimension of the multivariate
        totalMahanobis += maha[c];
    }
    // ---...--- //
    return exp(-0.5 * totalMahanobis);
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
    // ---- Parameters for the integral ----

    if (strcmp(method, "Genz") != 0 && strcmp(method, "Genz2") != 0)
    {
        inverseSymmetricPositiveMatrix(chol);
    }
    IntegrationParams params = {chol, mu};

    // ---- Initialize GSL Monte Carlo integration ----
    gsl_monte_function G = {&integral, (size_t)mvnDim, &params};
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    double result, error;
    // ---...--- //

    // ---- Perform integration ---- //
    if (strcmp(method, "Plain") == 0)
    {
        gsl_monte_plain_state *s = gsl_monte_plain_alloc(mvnDim);
        gsl_monte_plain_integrate(&G, lowerLimits, upperLimits, mvnDim, maxSamples, rng, s, &result, &error);
        gsl_monte_plain_free(s);
        gsl_rng_free(rng);
    }
    else if (strcmp(method, "Miser") == 0)
    {
        gsl_monte_miser_state *s = gsl_monte_miser_alloc(mvnDim);
        gsl_monte_miser_integrate(&G, lowerLimits, upperLimits, mvnDim, maxSamples, rng, s, &result, &error);
        gsl_monte_miser_free(s);
        gsl_rng_free(rng);
    }
    else if (strcmp(method, "Vegas") == 0)
    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(mvnDim);
        gsl_monte_vegas_integrate(&G, lowerLimits, upperLimits, mvnDim, maxSamples, rng, s, &result, &error);
        gsl_monte_vegas_free(s);
        gsl_rng_free(rng);
    }
    else if (strcmp(method, "Genz") == 0)
    {
        result = genzMontecarlo(chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
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
                "integral.\nThe method handed is:\t%s\nThe current available methods are `Plain`, `Miser`, `Vegas` and "
                "`Genz`.\n",
                method);
        exit(EXIT_FAILURE);
    }
    // ---...--- //

    // ---- Obtain the normalization values ---- //
    // ---- We'll get the determinant of the cholesky by multiplying its diagonals ----
    double determinant = 1;
    for (uint16_t c = 0; c < mvnDim; c++)
    {                                             // --- For each dimension of the multivariate
        determinant *= MATRIX_AT_PTR(chol, c, c); // Diagonal of the Cholensky
    }
    // ---...--- //

    // ---- Return the final value, normalized ---- //
    double denominator = sqrt(pow(2 * M_PI, mvnDim) * determinant);
    if (denominator == 0)
        return 0; // Early exit
                  // ---...--- //

    return (1 / denominator) * result;
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
double *computeQMultivariateCDF(Matrix const *probabilities, int monteCarloSamples, double epsilon, const char *method)
{

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
                    double lowerScaled = (featureCopyA[0]) / sqrt(MATRIX_AT_PTR(currentCholesky, 0, 0));
                    double upperScaled = (featureCopyB[0]) / sqrt(MATRIX_AT_PTR(currentCholesky, 0, 0));
                    montecarloResults[c] = (gsl_cdf_gaussian_P(upperScaled, 1) - gsl_cdf_gaussian_P(lowerScaled, 1)) *
                                           MATRIX_AT_PTR(probabilities, g, c);
                }
                denominator += montecarloResults[c];
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
