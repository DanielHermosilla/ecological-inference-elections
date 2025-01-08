#include "matrixUtils.h"
#include <cblas.h>
// Calculate Mahalanobis Distance
// https://en.wikipedia.org/wiki/Mahalanobis_distance
// Used to check similarity between multicovariate normals.
// x; observation
// mu; average
// cov_inv; inverse of the covariate matrix
//
// $$d = \sqrt{(\vec{x}-\vec{mu})\sigma^{-1}(\vec{x}-\vec{mu})}$$
//
void mahalanobis(double *x, double *mu, double *cov_inv, double *maha, int size)
{
    double diff[size];
    double temp[size];

    // Computes (x-mu)
    for (int i = 0; i < size; i++)
    {
        diff[i] = x[i] - mu[i];
    }

    // Computes in temp[size] the multiplication of (\sigma^{-1})*diff. It's a multiplication of a matrix with a vector.
    cblas_dsymv(CblasRowMajor, CblasUpper, size, 1.0, cov_inv, size, diff, 1, 0.0, temp, 1);

    // Computes the final product, that being (x-u)^T*temp. That's a dot product.
    double result = 0.0;
    for (int i = 0; i < size; i++)
    {
        result += diff[i] * temp[i];
    }
    *maha = result; // The result is squared since the approximate uses the squared result.
}

void compute_qm_mvn_pdf(double *n_trunc, double *p, double *p_trunc, double *b_m, double *diag_p, double *p_g_squared,
                        double *q_m)
{
    int I = TOTAL_CANDIDATES - 1;
    double mu[I - 1] = {0};
    double cov[(I - 1) * (I - 1)] = {0};
    double covs_U[G * (I - 1) * (I - 1)] = {0};
    double mus_U[G * (I - 1)] = {0};
    double vals_U[(I - 1)] = {0};           // Eigenvalues
    double vecs_U[(I - 1) * (I - 1)] = {0}; // Eigenvectors
    double maha[G * I] = {0};

    // Step 1: Calculate mu = b_m @ p_trunc
    cblas_dgemv(CblasRowMajor, CblasNoTrans, G, I - 1, 1.0, p_trunc, I - 1, b_m, 1, 0.0, mu, 1);

    // Step 2: Calculate Covariance Matrix
    for (int i = 0; i < I - 1; i++)
    {
        cov[i * (I - 1) + i] = mu[i] - p_trunc[i];
    }

    // Step 3: Eigen Decomposition
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', I - 1, cov, I - 1, vals_U);

    // Step 4: Mahalanobis Distance
    for (int g = 0; g < G; g++)
    {
        mahalanobis_distance(n_trunc, &mus_U[g * (I - 1)], cov, &maha[g * I], I - 1);
    }

    // Step 5: Probability Calculation
    for (int g = 0; g < G; g++)
    {
        for (int i = 0; i < I; i++)
        {
            q_m[g * I + i] = exp(-0.5 * maha[g * I + i]) * p[g * I + i];
        }
    }

    // Step 6: Normalize
    for (int g = 0; g < G; g++)
    {
        double sum = 0.0;
        for (int i = 0; i < I; i++)
        {
            sum += q_m[g * I + i];
        }
        for (int i = 0; i < I; i++)
        {
            q_m[g * I + i] /= sum;
        }
    }
}
