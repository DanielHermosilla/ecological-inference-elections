#include "multivariate-pdf.h"
#include <cblas.h>
#include <globals.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h> // Random numbers
#include <math.h>
#include <matrixUtils.h>
#include <multivariateUtils.h>
#include <stdint.h>
#include <string.h>

typedef struct
{
    Matrix *chol; // Cholesky decomposition, equivalent of the inverse matrix
    double *mu;   // The mu array
} IntegrationParams;

Matrix *featureMatrix = NULL;

/*
// Function to evaluate the integrand
double integral(double *x, size_t dim, void *params)
{
    IntegrationParams *p = (IntegrationParams *)params;
    const Matrix *chol = p->chol;
    const double *a = p->a;
    const double *b = p->b;
    const int ballotIndex = p->ballotIndex;
    const int candidateIndex = p->candidateIndex;
    int mvnDim = p->mvnDim;

    double value = 1.0;
    for (size_t i = 0; i < dim; i++)
    {
        // Transform sample x[i] to the scaled and shifted normal space
        double transformed = a[i] + MATRIX_AT_PTR(X, candidateIndex, ballotIndex) * (b[i] - a[i]);
        double cumulative = gsl_cdf_gaussian_P(transformed / MATRIX_AT_PTR(chol, i, i), 1.0);
        value *= cumulative;
    }
    return value;
}
*/
// Function to evaluate the integrand
// void getMahanalobisDist(double *x, double *mu, Matrix *inverseSigma, double *maha, int size);

double integral(double *x, size_t dim, void *params)
{
    IntegrationParams *p = (IntegrationParams *)params;
    Matrix *chol = p->chol;
    double *currentMu = p->mu;

    // ---- Obtain the mahanalobis of the multivariate ---- //
    double maha[dim];
    getMahanalobisDist(x, currentMu, chol, maha, dim, true);
    // ---...--- //

    // ---- Return the PDF ---- //
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

    // Perform integration
    double result, error;
    gsl_monte_plain_integrate(&G, lowerLimits, upperLimits, mvnDim, maxSamples, rng, s, &result, &error);

    // Free memory
    gsl_monte_plain_free(s);
    gsl_rng_free(rng);

    // ---- Obtain the normalization values ---- //
    // ---- We'll get the determinant of the cholesky by multiplying its diagonals ----
    double det = 1;
    for (uint16_t c = 0; c < mvnDim; c++)
    {
        det *= MATRIX_AT_PTR(chol, c, c);
    }

    return result;
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
void getFeatureMatrix()
{
    *featureMatrix = removeLastColumn(X); // c, b
}

double *computeQMultivariateCDF(Matrix const probabilities, int monteCarloSamples)
{

    Matrix probabilitiesReduced = removeLastColumn(&probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
                                                                                           // --- ... --- //

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        // ---- Get the values of the Multivariate that only depends on `b` ---- //
        Matrix mu = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
        Matrix **choleskyVals = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));
        getMainParameters(b, probabilitiesReduced, choleskyVals,
                          &mu); // TODO: FIX THIS, IT0S NOT POSSIBLE TO GET A INVERSE MATRIX IF ITS PERFECTLY LINEAL,
                                // NEED TO FIX THIS.
                                // ---...--- //

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            // ---- Define the current values to use ---- //
            Matrix *currentCholesky = choleskyVals[g];
            double *currentMu = getRow(&mu, TOTAL_GROUPS);
            double *feature = getColumn(featureMatrix, b); // Of size C-1
                                                           // ---...--- //

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                // ---- Define the borders of the hypercube ---- //
                // ---- First, make a copy of the feature vector ----
                double *featureCopyA = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
                double *featureCopyB = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));

                memcpy(featureCopyA, feature, TOTAL_CANDIDATES * sizeof(double));
                memcpy(featureCopyA, feature, TOTAL_CANDIDATES * sizeof(double));

                // ---- Substract/add and normalize ----
                // ---- Note that the bounds NEEDS to be standarized ----
                // ---- $$a=\sigma^{-1}(a-\mu)$$ ----
                for (uint16_t k = 0; k < TOTAL_CANDIDATES - 1; k++)
                {
                    featureCopyA[k] -= 0.5 - currentMu[k];
                    featureCopyB[k] += 0.5 - currentMu[k];
                    if (k == c)
                    {
                        featureCopyA[k] -= 1.0;
                        featureCopyB[k] -= 1.0;
                    }
                    featureCopyA[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                    featureCopyB[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                }

                // ---...--- //

                double monteCarloResult = Montecarlo(currentCholesky, currentMu, featureCopyA, featureCopyB,
                                                     (int)TOTAL_CANDIDATES - 1, monteCarloSamples);
            }
        }
    }
}
/**
 * @brief Calculates the CDF integral with Monte Carlo simulation.
 *
 * Given the hypercube parameters, it computes an approximation of the integral.
 *
 * @param[in] *chol Matrix of dimension (cxc) with the probabilities of each group and candidate.
 * @param[in] a First component of the unitary hypercube.
 * @param[in] b Second component of the unitary hypercube.
 * @param[in] mvnDim Dimension of the Multivariate Normal.
 * @param[in] epsilon Error threshold.
 * @param[in] iterations Maximum amount of iterations in the Monte Carlo method
 *
 * @return: An approximation of the CDF integral
 */
/*
double Montecarlo(const Matrix *chol, const double *a, const double *b, int mvnDim, double epsilon, int interations)
{
    // Initialize RNG
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 42); // seed

    size_t intsum = 0;
    size_t varsum = 0;
    size_t totalSim = 0;

    double *d = (double *)malloc(mvnDim * sizeof(double));
    double *e = (double *)malloc(mvnDim * sizeof(double));
    double *f = (double *)malloc(mvnDim * sizeof(double));
    double *y = (double *)malloc(mvnDim * sizeof(double));

    // Initialize d, e, and f for the first dimension
    d[0] = gsl_cdf_gaussian_P(a[0] / MATRIX_AT_PTR(chol, 0, 0), 1.0); // Standard normal CDF
    e[0] = gsl_cdf_gaussian_P(b[0] / MATRIX_AT_PTR(chol, 0, 0), 1.0);
    f[0] = e[0] - d[0];
}
*/
