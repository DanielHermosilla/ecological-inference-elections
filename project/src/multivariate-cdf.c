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

typedef struct
{
    const Matrix *chol;       // Cholesky decomposition, equivalent of the inverse matrix
    const int ballotIndex;    // Index from the ballot
    const int candidateIndex; // Index from the candidate
    const double *a;          // Lower bounds of the hypercube
    const double *b;          // Upper bounds of the hypercube
    int mvnDim;               // Dimension of the MVN
} IntegrationParams;

// Function to evaluate the integrand
double mvn_integrand(double *x, size_t dim, void *params)
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

// Monte Carlo CDF approximation
double Montecarlo(const Matrix *chol, int ballotIndex, int candidateIndex, const double *a, const double *b, int mvnDim,
                  double epsilon, int maxSamples)
{
    // Define integration limits
    double xl[mvnDim], xu[mvnDim];
    for (int i = 0; i < mvnDim; i++)
    {
        xl[i] = 0.0; // Lower bound for Monte Carlo sampling
        xu[i] = 1.0; // Upper bound for Monte Carlo sampling
    }

    // Set up the parameters for the integrand
    IntegrationParams params = {
        chol, ballotIndex, candidateIndex, a, b, mvnDim,
    };

    // Initialize GSL Monte Carlo integration
    gsl_monte_function G = {&mvn_integrand, (size_t)mvnDim, &params};
    // {&mvn_integrand, (size_t)mvnDim, &params};
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_plain_state *s = gsl_monte_plain_alloc(mvnDim);

    // Perform integration
    double result, error;
    gsl_monte_plain_integrate(&G, xl, xu, mvnDim, maxSamples, rng, s, &result, &error);

    // Free memory
    gsl_monte_plain_free(s);
    gsl_rng_free(rng);

    return result;
}

void getMainParameters(int b, Matrix const probabilitiesReduced)
{
    // --- Get the mu and sigma --- //
    Matrix muR = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
    Matrix **sigma = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        sigma[g] = (Matrix *)malloc(sizeof(Matrix));
        *sigma[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(b, &probabilitiesReduced, &muR, sigma);
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(sigma[g]);
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
