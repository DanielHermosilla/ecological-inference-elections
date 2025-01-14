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
#include <matrixUtils.h>
#include <multivariateUtils.h>
#include <stdint.h>
#include <string.h>

typedef struct
{
    const Matrix *chol;       // Cholesky decomposition, equivalent of the inverse matrix
    const int ballotIndex;    // Index from the ballot
    const int candidateIndex; // Index from the candidate
    const double *a;          // Lower bounds of the hypercube
    const double *b;          // Upper bounds of the hypercube
    int mvnDim;               // Dimension of the MVN
    const double *feature;    // The feature vector
} IntegrationParams;

Matrix *featureMatrix = NULL;

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

    // ---- Obtain the PDF of the multivariate ---- //
    double *cholDotProduct = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
    cblas_dgemm(CblasRowMajor, // Matrix storage order
                CblasNoTrans,  // Chol is not transposed
                CblasNoTrans,  // Y is not transposed
                chol->rows,    // Number of rows in Chol
                feature->cols, // Number of columns in Y
                k,             // Shared dimension (columns of Chol, rows of Y)
                1.0,           // Alpha (scaling factor for Chol * Y)
                chol, k,       // Chol matrix and leading dimension
                y, n,          // Y matrix and leading dimension
                0.0,           // Beta (scaling factor for Chol_cum)
                chol_cum, n);  // Output matrix and leading dimension

    for (size_t c = 0; c < dim; c++)
    {
    }

    return value;
}
// Monte Carlo CDF approximation
// Refer to https://www.gnu.org/software/gsl/doc/html/montecarlo.html
double Montecarlo(const Matrix *chol, int ballotIndex, int candidateIndex, double *feature, const double *lowerLimits,
                  const double *upperLimits, int mvnDim, int maxSamples)
{
    /*
    // Define integration limits
    double xl[mvnDim], xu[mvnDim];
    for (int i = 0; i < mvnDim; i++)
    {
        xl[i] = 0.0; // Lower bound for Monte Carlo sampling
        xu[i] = 1.0; // Upper bound for Monte Carlo sampling
    }
    */

    // Set up the parameters for the integrand
    IntegrationParams params = {chol, ballotIndex, candidateIndex, lowerLimits, upperLimits, mvnDim, feature};

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

    getAverageConditionalFull(b, &probabilitiesReduced, mu, cholesky);
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(cholesky[g]);
    }
}
/*
void getFeatureMatrix(){
    *featureMatrix = createMatrix(X->rows, X->cols); // c, b
    memcpy(featureMatrix->data, X->data, X->rows * X->cols * sizeof(double));
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++){
        for (uint32_t b = 0; b < TOTAL_BALLOTS; b++){

        }
    }
    double *returningArray = (double *)malloc(TOTAL_CANDIDATES-1 * sizeof(double));
    returningArray = getColumn(XwithoutLast, b);
    return returningArray;
}
*/
double *computeQMultivariateCDF(Matrix const probabilities, int monteCarloSamples)
{

    Matrix probabilitiesReduced = removeLastColumn(&probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
                                                                                           // --- ... --- //

    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {

            // ---- Get the values of the Multivariate ---- //
            Matrix mu = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
            Matrix **choleskyVals = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));
            getMainParameters(b, probabilitiesReduced, choleskyVals, &mu);
            // ---...--- //
            // ---- Define the feature vector, the same for every group
            double *feature = getColumn(X, b);

            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                Matrix *currentCholesky = choleskyVals[g];
                // ---- Define the borders of the hypercube ---- //
                // ---- First, make a copy of the feature vector ----
                double *featureCopy = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
                double *featureCopyA = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
                double *featureCopyB = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));

                memcpy(featureCopy, feature, TOTAL_CANDIDATES * sizeof(double));
                memcpy(featureCopyA, feature, TOTAL_CANDIDATES * sizeof(double));
                memcpy(featureCopyA, feature, TOTAL_CANDIDATES * sizeof(double));

                // ---- Substract/add and normalize ----
                // ---- Note that the bounds NEEDS to be standarized ----
                for (uint16_t k = 0; k < TOTAL_CANDIDATES; k++)
                {
                    featureCopyA[k] -= 0.5 - MATRIX_AT(mu, g, k);
                    featureCopyB[k] += 0.5 - MATRIX_AT(mu, g, k);
                    if (k == c)
                    {
                        featureCopyA[k] -= 1.0;
                        featureCopyB[k] -= 1.0;
                    }
                    featureCopyA[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                    featureCopyB[k] *= MATRIX_AT_PTR(currentCholesky, k, k);
                }

                // ---- Substract the candidate index if it's from the `C-1` candidates
                if (c < TOTAL_CANDIDATES - 1)
                {
                    featureCopyA[c] -= 1;
                    featureCopyB[c] -= 1;
                }
                // ---...--- //

                Matrix *currentCholesky = choleskyVals[g];
                double monteCarloResult = Montecarlo(currentCholesky, b, c, featureCopy, featureCopyA, featureCopyB,
                                                     (int)TOTAL_CANDIDATES, monteCarloSamples);
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
