#include <gsl/gsl_rng.h> // Random numbers
#include <matrixUtils.h>
#include <multivariateUtils.h>

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
double Montecarlo(const Matrix *chol, int a, int b, int mvnDim, double epsilon, int interations)
{
    // Initialize RNG
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 42); // seed
}
