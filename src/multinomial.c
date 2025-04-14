/*
Copyright (c) 2025 Daniel Hermosilla

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

#include "multinomial.h"
#include "globals.h"
#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <Rmath.h>
#include <stdint.h>

#ifndef BLAS_INT
#define BLAS_INT int
#endif

double *logGammaArr2 = NULL;

/**
 * @brief Computes the value of `r` without the denominator.
 *
 * Given that the probabilities cancel the denominator, it just computes the numerator of the `r` definition.
 *
 * @param[in] *probabilities Matrix of dimension (gxc) with the probabilities of each group and candidate.
 * @param[in] *mult Matrix of dimension (bxc) with the matricial multiplication of w * p.
 * @param[in] b The index `b`
 * @param[in] c The index `c`
 * @param[in] g The index `g`
 *
 * @return: The value for `r` at the position given.
 */
double computeR(Matrix const *probabilities, Matrix const *mult, int const b, int const c, int const g)
{

    return MATRIX_AT_PTR(mult, b, c) - MATRIX_AT_PTR(probabilities, g, c);
}

/*
 * Computes the WHOLE log-likelihood using the Multinomial method.
 *
 * @note: The *ll of the original algorithm will take care of the \Delta log-likelihood, not needing
 * to compute the lgamma function
 *
 * This would be called outside the computing function, so it is not designed to "save" some calculations.
 */
static void precomputeLogGammas()
{
    // We must get the biggest W_{bg}
    int biggestX = (int)maxElement(X);
    logGammaArr2 = (double *)Calloc(biggestX + 1, double);

    for (int i = 0; i <= biggestX; i++)
    {
        logGammaArr2[i] = lgamma1p(i);
    }
}

/**
 * @brief Computes an approximate of the conditional probability by using a Multinomial approach.
 *
 * Given the observables parameters and the probability matrix, it computes an approximation of `q` with the Multinomial
 * approach.
 *
 * @param[in] *probabilities Matrix of dimension (gxc) with the probabilities of each group and candidate.
 * @param[in] params The optional parameters to the function. For this specific method, there's no supported optional
 * parameters yet.
 *
 * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices, but
 * it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
 *
 */
double *computeQMultinomial(Matrix const *probabilities, QMethodInput params, double *ll)
{

    if (logGammaArr2 == NULL)
    {
        precomputeLogGammas();
    }
    double *array2 = (double *)Calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, double); // Array to return
    *ll = 0;
    // -- Summatory calculation for g --
    // This is a simple matrix calculation, to be computed once.
    Matrix WP = createMatrix((int)TOTAL_BALLOTS, (int)TOTAL_CANDIDATES);

    double alpha = 1.0;
    double beta = 0.0;
    BLAS_INT m = (BLAS_INT)TOTAL_BALLOTS;
    BLAS_INT k = (BLAS_INT)TOTAL_GROUPS;
    BLAS_INT n = (BLAS_INT)TOTAL_CANDIDATES;
    char noTranspose = 'N';

    // WP = alpha * W * probabilities + beta * WP
    F77_CALL(dgemm)
    (&noTranspose, &noTranspose,    // transA = 'N', transB = 'N'
     &m, &n, &k,                    // M, N, K
     &alpha, W->data, &m,           // A, LDA = m
     probabilities->data, &k,       // B, LDB = k
     &beta, WP.data, &m FCONE FCONE // C, LDC = m
                                    // string lengths for 'N', 'N'
    );

    // ---- Do not parallelize ----
    double totalWP[TOTAL_BALLOTS];

    for (int b = 0; b < (int)TOTAL_BALLOTS; b++)
    {
        totalWP[b] = 0;
        for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            // Because it cannot be initialized
            totalWP[b] += MATRIX_AT(WP, b, c);
    }

    for (int b = 0; b < (int)TOTAL_BALLOTS; b++)
    { // --- For each ballot box
        for (int g = 0; g < (int)TOTAL_GROUPS; g++)
        { // --- For each group given a ballot box
            // ---- Create temporal variables ----
            double tempSum = 0.0;
            double finalNumerator[TOTAL_CANDIDATES];

            for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            { // --- For each candidate given a group and a ballot box

                // ---- Compute x*p*r^{-1} ---- //
                double numerator = MATRIX_AT_PTR(probabilities, g, c) * (int)MATRIX_AT_PTR(X, c, b);
                double denominator = computeR(probabilities, &WP, b, c, g);

                finalNumerator[c] = denominator != 0 ? numerator / denominator : 0;

                // ---- Store the value for reusing it later ----
                tempSum += finalNumerator[c];
                // ---...--- //
                // Add the log-likelihood
                if (g == 0)
                {
                    *ll += MATRIX_AT(WP, b, c) != 0 ? MATRIX_AT_PTR(X, c, b) * log(MATRIX_AT(WP, b, c) / totalWP[b]) -
                                                          lgamma1p((int)MATRIX_AT_PTR(X, c, b))
                                                    : 0;
                    Rprintf("Término S: %f Término de X_cb: %f\n",
                            MATRIX_AT_PTR(X, c, b) * log(MATRIX_AT(WP, b, c) / totalWP[b]),
                            lgamma1p((int)MATRIX_AT_PTR(X, c, b)));
                }
            }

            for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            { // ---- For each candidate given a group and a ballot box
                // ---- Store the value ----
                Q_3D(array2, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = tempSum != 0 ? finalNumerator[c] / tempSum : 0;
            }
        }
        Rprintf("Sumando el término 'n': %f\n------- Terminando iteración B = %d -------\n", lgamma1p(BALLOTS_VOTES[b]),
                b);
        *ll += lgamma1p(BALLOTS_VOTES[b]);
    }
    // *ll -= TOTAL_BALLOTS * TOTAL_CANDIDATES * log(totalWP);
    freeMatrix(&WP);
    return array2;
}
