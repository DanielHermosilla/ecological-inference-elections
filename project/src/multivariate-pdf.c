#include "multivariate-pdf.h"
#include "globals.h"
#include "matrixUtils.h"
#include <cblas.h>
#include <lapacke.h>
#include <math.h>
#include <unistd.h>

/**
 * @brief Computes the `q` values for a given ballot box.
 *
 * Given a ballot box index, probabilities and the reduced version (with C-1 candidates) of the probabilities matrix, it
 * calculates the `q` values in a flattened way
 *
 * @param[in] b. The index of the ballot box
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 * @param[in] *probabilitiesReduced. A pointer towards the reduced probabilities matrix.
 *
 * @return A (g x c) matrix with the values of `q` according the candidate and group index.
 *
 */

Matrix computeQforABallot(int b, const Matrix *probabilities, const Matrix *probabilitiesReduced)
{

    // --- Get the mu and sigma --- //
    Matrix muR = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
    Matrix **sigma = (Matrix **)malloc(TOTAL_GROUPS * sizeof(Matrix *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        sigma[g] = (Matrix *)malloc(sizeof(Matrix));
        *sigma[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(b, probabilitiesReduced, &muR, sigma);
    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        // ---- Calculates the Cholensky matrix ---- //
        inverseSymmetricPositiveMatrix(sigma[g]);
    }
    // ---- ... ----

    // --- Calculate the mahanalobis distance --- //
    double **mahanalobisDistances = (double **)malloc(TOTAL_GROUPS * sizeof(double *));

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        mahanalobisDistances[g] = (double *)malloc(TOTAL_CANDIDATES * sizeof(double));
        // ---- Get the feature vector (the candidate results) ----
        double *feature = getColumn(X, b);
        // ---- Get the average values for the candidate results ----
        double *muG = getRow(&muR, g);
        // ---- Call the mahanalobis function ----
        getMahanalobisDist(feature, muG, sigma[g], mahanalobisDistances[g], TOTAL_CANDIDATES - 1, false);
        // ---- Free allocated and temporary values ----
        freeMatrix(sigma[g]);
        free(sigma[g]);
        free(feature);
        free(muG);
    }
    free(sigma);
    freeMatrix(&muR);
    // --- .... --- //

    // --- Calculate the returning values --- //
    // ---- Create the matrix to return ----
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        // ---- Initialize variables ----
        double den = 0;
        double *QC = (double *)calloc(TOTAL_CANDIDATES, sizeof(double)); // Value of Q on candidate C

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            // ---- The `q` value is calculated as exp(-0.5 * mahanalobis) * probabilities ----
            QC[c] = exp(-0.5 * mahanalobisDistances[g][c]) * MATRIX_AT_PTR(probabilities, g, c);
            // ---- Add the values towards the denominator to later divide by it ----
            den += QC[c];
        }
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            // ---- Store each value, divided by the denominator ----
            MATRIX_AT(toReturn, g, c) = QC[c] / den;
        }
        // ---- Free allocated memory ----
        free(mahanalobisDistances[g]);
        free(QC);
    }
    free(mahanalobisDistances); // Might have to remove

    return toReturn;
    // --- ... --- //
}

/**
 * @brief Computes the `q` values for all the ballot boxes given a probability matrix. Uses the Multivariate PDF method.
 *
 * Given a probability matrix with, it returns a flattened array with estimations of the conditional probability. The
 * array can be accesed with the macro `Q_3D` (it's a flattened tensor).
 *
 * @param[in] *probabilities. A pointer towards the probabilities matrix.
 *
 * @return A pointer towards the flattened tensor.
 *
 */
double *computeQMultivariatePDF(Matrix const *probabilities)
{
    // ---- Initialize values ---- //
    // ---- The probabilities without the last column will be used for each iteration ----
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return
                                                                                           // --- ... --- //

    // ---- Fill the array with the results ---- //
#pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot
        // ---- Call the function for calculating the `q` results for a given ballot
        Matrix resultsForB = computeQforABallot((int)b, probabilities, &probabilitiesReduced);
        // #pragma omp parallel for collapse(2)
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group given a ballot box
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate given a ballot box and a group
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = MATRIX_AT(resultsForB, g, c);
            }
        }
        // ---- Frees allocated space ----
        freeMatrix(&resultsForB);
    }

    freeMatrix(&probabilitiesReduced);
    return array2;
    // --- ... --- //
}
