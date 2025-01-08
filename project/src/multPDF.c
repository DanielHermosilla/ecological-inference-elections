#include "matrixUtils.h"

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
