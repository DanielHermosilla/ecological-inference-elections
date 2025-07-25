/*
Copyright (c) 2025 fastei team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "multivariate-pdf.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

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

Matrix computeQforABallot(EMContext *ctx, int b, const Matrix *probabilities, const Matrix *probabilitiesReduced,
                          double *ll)
{
    Matrix *X = &ctx->X;

    // --- Get the mu and sigma --- //
    Matrix muR = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES - 1);
    Matrix **sigma = (Matrix **)Calloc(TOTAL_GROUPS, Matrix *);

    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        sigma[g] = (Matrix *)Calloc(1, Matrix);
        *sigma[g] = createMatrix(TOTAL_CANDIDATES - 1, TOTAL_CANDIDATES - 1); // Initialize
    }

    getAverageConditional(ctx, b, probabilitiesReduced, &muR, sigma);

    // ---- ... ----

    // ---- Get the inverse matrix for each sigma ---- //
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        inverseSymmetricPositiveMatrix(sigma[g]);
    }
    // ---- ... ----
    // ---- Get the determinant FOR THE LOG-LIKELIHOOD ---- //
    double det = 1;
    for (uint16_t c = 0; c < TOTAL_CANDIDATES - 1; c++)
    {
        det *= MATRIX_AT_PTR(sigma[0], c, c);
    }
    det = 1.0 / (det * det);
    double normalizeConstant = R_pow(R_pow_di(M_2_PI, (int)(TOTAL_CANDIDATES - 1)) * det, 0.5);

    // --- Calculate the mahanalobis distance --- //
    double **mahanalobisDistances = (double **)Calloc(TOTAL_GROUPS, double *);

    // // #pragma omp parallel for
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group ----
        mahanalobisDistances[g] = (double *)Calloc(TOTAL_CANDIDATES, double);
        // ---- Get the feature vector (the candidate results) ----
        double *feature = getColumn(X, b);
        // ---- Get the average values for the candidate results ----
        double *muG = getRow(&muR, g);
        // ---- Call the mahanalobis function ----
        getMahanalobisDist(feature, muG, sigma[g], mahanalobisDistances[g], TOTAL_CANDIDATES - 1, false);
        // ---- Free allocated and temporary values ----
        freeMatrix(sigma[g]);
        Free(sigma[g]);
        Free(feature);
        Free(muG);
    }
    Free(sigma);
    freeMatrix(&muR);
    // --- .... --- //

    // --- Calculate the returning values --- //
    // ---- Create the matrix to return ----
    Matrix toReturn = createMatrix(TOTAL_GROUPS, TOTAL_CANDIDATES);
    for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
    { // ---- For each group
        // ---- Initialize variables ----
        double den = 0;
        double *QC = (double *)Calloc(TOTAL_CANDIDATES, double); // Value of Q on candidate C

        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            // ---- The `q` value is calculated as exp(-0.5 * mahanalobis) * probabilities ----
            QC[c] = exp(-0.5 * mahanalobisDistances[g][c]) * MATRIX_AT_PTR(probabilities, g, c);
            // ---- Add the values towards the denominator to later divide by it ----
            den += QC[c];
        }
        *ll += g == 0 && den > 0 ? log(den) * normalizeConstant : 0;
        for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
        { // ---- For each candidate given a group
            // ---- Store each value, divided by the denominator ----
            if (den != 0)
                MATRIX_AT(toReturn, g, c) = QC[c] / den;
            else
                MATRIX_AT(toReturn, g, c) = 0;
        }
        // ---- Free allocated memory ----
        Free(mahanalobisDistances[g]);
        Free(QC);
    }
    Free(mahanalobisDistances); // Might have to remove

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
void computeQMultivariatePDF(EMContext *ctx, QMethodInput params, double *ll)
{
    // ---- Initialize values ---- //
    *ll = 0;
    Matrix *X = &ctx->X;
    Matrix *W = &ctx->W;
    IntMatrix *intX = &ctx->intX;
    IntMatrix *intW = &ctx->intW;
    double *q = ctx->q;
    Matrix *probabilities = &ctx->probabilities;
    // ---- The probabilities without the last column will be used for each iteration ----

    // ---- The idea is to remove the column with the most votes, so it'll swap it to the last column and later reswap
    // it. ----
    Matrix probabilitiesReduced = removeLastColumn(probabilities);
    // ---- Fill the array with the results ---- //
    // #pragma omp parallel for
    for (uint32_t b = 0; b < TOTAL_BALLOTS; b++)
    { // ---- For each ballot
        // ---- Call the function for calculating the `q` results for a given ballot
        Matrix resultsForB = computeQforABallot(ctx, (int)b, probabilities, &probabilitiesReduced, ll);
        // TODO: Optimize the return
        // #pragma omp parallel for collapse(2)
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        { // ---- For each group given a ballot box
            for (uint16_t c = 0; c < TOTAL_CANDIDATES; c++)
            { // ---- For each candidate given a ballot box and a group
                double results = MATRIX_AT(resultsForB, g, c);

                Q_3D(q, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) =
                    !isnan(results) && !isinf(results) ? results : 0;
            }
        }
        // ---- Frees allocated space ----
        freeMatrix(&resultsForB);
    }

    freeMatrix(&probabilitiesReduced);
    if (isnan(*ll) || isinf(*ll))
        *ll = 0;
    // --- ... --- //
}
