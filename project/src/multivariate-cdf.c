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

void getMainParameters(int b, Matrix const probabilitiesReduced, Matrix **cholesky, Matrix *mu)
{

    // --- Get the mu and sigma --- //

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        cholesky[g] = (Matrix *)malloc(sizeof(Matrix));
        *cholesky[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(b, &probabilitiesReduced, mu, cholesky);
    // ---- ... ----
    // Missing sigma is

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(cholesky[g]);
    }
}

double *computeQMultivariateCDF(Matrix const *probabilities, int monteCarloSamples)
{

    /*
    if (featureMatrix == NULL)
    {
        printf("Creating feature matrix...\n");
        featureMatrix = removeLastColumn(X); // c, b
        printf("Matrix created\n");
    }
    */
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
                                                                                           // --- ... --- //

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        // ---- Get the values of the Multivariate that only depends on `b` ---- //
        Matrix mu = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
        Matrix **choleskyVals = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));
        getMainParameters(b, probabilitiesReduced, choleskyVals,
                          &mu); // TODO: FIX THIS, IT0S NOT POSSIBLE TO GET A INVERSE MATRIX IF ITS PERFECTLY LINEAL,
                                // NEED TO FIX THIS.
                                // ---...--- //
        double *feature = getColumn(X, b); // Of size C-1

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            // ---- Define the current values to use ---- //
            Matrix *currentCholesky = choleskyVals[g];
            double *currentMu = getRow(&mu, g);
            // ---...--- //

            double montecarloResults[TOTAL_CANDIDATES];
            double denominator = 0;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {

                // Define the integral

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
                {
                    featureCopyA[k] -= 0.5;
                    featureCopyB[k] += 0.5;
                    // featureCopyA[k] -= currentMu[k];
                    // featureCopyB[k] -= currentMu[k];
                    if (k == c)
                    {
                        featureCopyA[k] -= 1.0;
                        featureCopyB[k] -= 1.0;
                    }
                    // featureCopyA[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                    // featureCopyB[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                }

                // ---...--- //
                montecarloResults[c] = Montecarlo(currentCholesky, currentMu, featureCopyA, featureCopyB,
                                                  (int)TOTAL_CANDIDATES - 1, monteCarloSamples) *
                                       MATRIX_AT_PTR(probabilities, g, c);

                denominator += montecarloResults[c];

                free(featureCopyA);
                free(featureCopyB);
            } // --- End c loop
            free(currentMu);
            freeMatrix(currentCholesky);

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                if (denominator == 0) // Edge case
                    Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = 0;
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = montecarloResults[c] / denominator;
            } // --- End c loop
        } // --- End g loop
        free(feature);
        freeMatrix(&mu);
        free(choleskyVals);
    } // --- End b loop
    freeMatrix(&probabilitiesReduced);
    return array2;
}

__attribute__((destructor)) void freeFeatures()
{
    freeMatrix(featureMatrix);
}
