#include "multivariate-cdf.h"
#include "globals.h"
#include <cblas.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h> // Random numbers
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

typedef struct
{
    Matrix *chol; // Cholesky decomposition, equivalent of the inverse matrix
    double *mu;   // The mu array
} IntegrationParams;

Matrix *featureMatrix = NULL;

// double getLimits(double *a, double *b)
double integral(double *x, size_t dim, void *params)
{
    IntegrationParams *p = (IntegrationParams *)params;
    Matrix *chol = p->chol;
    double *currentMu = p->mu;

    // ---- Obtain the mahanalobis of the multivariate ---- //
    double maha[dim];
    getMahanalobisDist(x, currentMu, chol, maha, dim, true);
    // ---...--- //

    // ---- the PDF ---- //
    // ---- exp(-0.5 * maha) ----
    // ---- To get the scalar mahanalobis we need to sum over all contributions ----
    double totalMahanobis = 0;
    for (uint16_t c = 0; c < dim; c++)
    {
        totalMahanobis += maha[c];
    }
    return exp(-0.5 * totalMahanobis);
    // ---...--- //
}

// Monte Carlo CDF approximation
// Refer to https://www.gnu.org/software/gsl/doc/html/montecarlo.html
double Montecarlo(Matrix *chol, double *mu, const double *lowerLimits, const double *upperLimits, int mvnDim,
                  int maxSamples)
{
    // Set up the parameters for the integrand
    IntegrationParams params = {chol, mu};

    // Initialize GSL Monte Carlo integration
    gsl_monte_function G = {&integral, (size_t)mvnDim, &params};
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_plain_state *s = gsl_monte_plain_alloc(mvnDim);
    gsl_rng_set(rng, 42);
    // Perform integration
    double result, error;
    // sleep(5);
    gsl_monte_plain_integrate(&G, lowerLimits, upperLimits, mvnDim, maxSamples, rng, s, &result, &error);
    // Free memory
    gsl_monte_plain_free(s);
    gsl_rng_free(rng);

    // ---- Obtain the normalization values ---- //
    // ---- We'll get the determinant of the cholesky by multiplying its diagonals ----
    double determinant = 1;
    for (uint16_t c = 0; c < mvnDim; c++)
    {
        determinant *= MATRIX_AT_PTR(chol, c, c);
    }
    double denominator = sqrt(pow(2 * M_PI, mvnDim) * determinant);
    if (denominator == 0)
        return 0; // Early exit
                  // ---...--- //

    return (1 / denominator) * result;
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
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(cholesky[g]);
    }
    // ---...--- //
}
/*
 * @brief Computes the `q` values from the multivariate CDF method
 *
 * @param[in] *probabilities A pointer towards the current probabilities matrix.
 * @param[in] monteCarloSamples The amount of samples to use in the Monte Carlo simulation
 *
 * @return A contiguos array with all the new probabilities
 */
double *computeQMultivariateCDF(Matrix const *probabilities, int monteCarloSamples)
{

    // ---- Define initial variables ---- //
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
                                                                                           // --- ... --- //

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
                }
                // ---...--- //
                // ---- Save the results and add them to the denominator ---- //
                montecarloResults[c] = Montecarlo(currentCholesky, currentMu, featureCopyA, featureCopyB,
                                                  (int)TOTAL_CANDIDATES - 1, monteCarloSamples) *
                                       MATRIX_AT_PTR(probabilities, g, c);

                denominator += montecarloResults[c];
                free(featureCopyA);
                free(featureCopyB);
                // ---...--- //
            } // --- End c loop
            free(currentMu);
            freeMatrix(currentCholesky);

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
