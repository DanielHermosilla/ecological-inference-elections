#include "multPDF.h"
#include "globals.h"
#include "matrixUtils.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
// Calculate Mahalanobis Distance
// https://en.wikipedia.org/wiki/Mahalanobis_distance
// Used to check similarity between multicovariate normals.
// x; observation
// mu; average
// inverseSigma; inverse of the covariate matrix
//
// $$d = \sqrt{(\vec{x}-\vec{mu})\sigma^{-1}(\vec{x}-\vec{mu})}$$
//

/**
 * @brief Computes the Mahalanobis distance with last candidate adjustment.
 *
 * @param[in] x Pointer to the input feature vector (size C-1).
 * @param[in] mu Pointer to the mean vector (size C-1).
 * @param[in] inverseSigma Pointer to the inverse covariance matrix (size (C-1) x (C-1)).
 * @param[out] maha Pointer to the resulting Mahalanobis distances (size C).
 * @param[in] size Size of the truncated candidate space (C-1).
 */
void mahanalobis(double *x, double *mu, Matrix *inverseSigma, double *maha, int size)
{
    double diff[size];
    double temp[size];
    double invs_devs[size];
    double diagonalInverse[size];

    // Computes (x - mu)
    for (int i = 0; i < size; i++)
    {
        diff[i] = x[i] - mu[i];
    }

    // Compute invs_devs (inverseSigma * diff)
    cblas_dsymv(CblasRowMajor, CblasUpper, size, 1.0, inverseSigma->data, size, diff, 1, 0.0, temp, 1);

    // Compute Mahalanobis Distance (truncated)
    double mahanobisTruncated = 0.0;
    for (int i = 0; i < size; i++)
    {
        // The second parenthesis
        mahanobisTruncated += diff[i] * temp[i];
        invs_devs[i] = temp[i]; // Store intermediate results
    }

    maha[size] = mahanobisTruncated; // Last element is used as a reference

    // Extract diagonal inverses
    for (int i = 0; i < size; i++)
    {
        diagonalInverse[i] = MATRIX_AT_PTR(inverseSigma, i, i);
    }

    // Correct Mahalanobis Distance for all candidates
    for (int i = 0; i < size; i++)
    {
        maha[i] = maha[size] - 2 * invs_devs[i] + diagonalInverse[i];
    }
}

double *computeQforG(int b, int g, const Matrix *probabilities, const Matrix *probabilitiesReduced)
{

    double *candidateVotesPerBallot = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
    double *mahalanobisArray = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));

    // --- Get the mu and sigma --- All of this depends of "g" and "b"
    double *mu = (double *)malloc((TOTAL_CANDIDATES - 1) * sizeof(double));
    Matrix sigma = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1);
    getAverageConditional(b, g, probabilitiesReduced, mu, &sigma);
    // printf("The sigma matrix before inverse is\n");
    // printMatrix(&sigma);
    inverseMatrixLU(&sigma); // Cholensky decomposition for getting sigma ^{-1}
    // printf("After inverse is\n");
    // printMatrix(&sigma);
    //  --- ... ----

    // --- Calculate the mahanalobis distance --- //
    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        candidateVotesPerBallot[c] = MATRIX_AT_PTR(X, c, b);
    }
    mahanalobis(candidateVotesPerBallot, mu, &sigma, mahalanobisArray,
                TOTAL_CANDIDATES - 1); // Vector of size "c" with all the distances
    // --- ... ---

    double *toReturn = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));

    double den = 0;
    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        toReturn[c] = exp(-0.5 * mahalanobisArray[c]) * MATRIX_AT_PTR(probabilities, g, c);
        den += toReturn[c];
    }

    for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
    {
        toReturn[c] /= den;
    }

    free(candidateVotesPerBallot);
    free(mahalanobisArray);
    return toReturn;
}

double *computeQMultivariatePDF(Matrix const *probabilities)
{
    printf("The most original probability matrix is");
    printMatrix(probabilities);
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return

    // #pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    {

        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            double *toInsert = computeQforG(b, g, probabilities, &probabilitiesReduced);
            return array2;
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                // printf("\nThe value to insert is:\t%.3f\n", toInsert[c]);
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = toInsert[c];
            }
            free(toInsert);
        }
    }

    freeMatrix(&probabilitiesReduced);
    return array2;
}
