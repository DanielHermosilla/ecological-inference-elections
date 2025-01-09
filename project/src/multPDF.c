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

Matrix computeQforABallot(int b, const Matrix *probabilities, const Matrix *probabilitiesReduced)
{

    // --- Get the mu and sigma --- Arrays and matrices since it depends on g
    Matrix muR = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
    Matrix **sigma = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        sigma[g] = (Matrix *)malloc(sizeof(Matrix));
        *sigma[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(b, probabilitiesReduced, &muR, sigma);
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        inverseMatrixLU(sigma[g]);
    }
    // ---- ... ----
    // printf("After inverse is\n");
    // printMatrix(&sigma);
    //  --- ... ----

    // --- Calculate the mahanalobis distance --- //
    double **mahanalobisDistances = (double **)malloc(TOTAL_GROUPS * sizeof(double *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        mahanalobisDistances[g] = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
        double *feature = getColumn(X, b);
        double *muG = getRow(&muR, g);
        mahanalobis(feature, muG, sigma[g], mahanalobisDistances[g], TOTAL_CANDIDATES - 1);
        freeMatrix(sigma[g]);
        free(feature);
        free(muG);
    }
    free(sigma);
    freeMatrix(&muR);
    // --- .... ---

    // --- Calculate the returning values ---
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    {
        double den = 0;
        double *QC = (double *)calloc(TOTAL_CANDIDATES, sizeof(double)); // Value of Q on candidate C

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {

            QC[c] = exp(-0.5 * mahanalobisDistances[g][c]) * MATRIX_AT_PTR(probabilities, g, c);
            den += QC[c];
        }
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        {
            MATRIX_AT(toReturn, g, c) = QC[c] / den;
        }
        free(mahanalobisDistances[g]);
        free(QC);
    }
    free(mahanalobisDistances); // Might have to remove

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

        Matrix resultsForB = computeQforABallot((int)b, probabilities, &probabilitiesReduced);
        printf("\nThe results for %d ballot is:\n", b);
        printMatrix(&resultsForB);
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            {
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = MATRIX_AT(resultsForB, g, c);
            }
        }
        freeMatrix(&resultsForB);
    }

    freeMatrix(&probabilitiesReduced);
    printf("The array to return is:\n");
    printf("[");
    for (int a = 0; a < (int)TOTAL_BALLOTS * (int)TOTAL_CANDIDATES * (int)TOTAL_GROUPS; a++)
    {
        printf("%.4f, ", array2[a]);
    }
    printf("]");
    return array2;
}
