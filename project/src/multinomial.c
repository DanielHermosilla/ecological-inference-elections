#include "multinomial.h"
#include "matrixUtils.h"
#include <cblas.h>

/**
 * Access the element (b, g, c) from a row-major 3D array `q` with dimensions:
 *   b-dimension  = total number of ballot boxes
 *   g-dimension  = total number of groups
 *   c-dimension  = total number of candidates
 *
 * \param q     A pointer to the double array of size (b * g * c)
 * \param bIdx  The ballot index
 * \param gIdx  The group  index
 * \param cIdx  The candidate index
 * \param G     The total number of groups
 * \param C     The total number of candidates
 */
#define Q_3D(q, bIdx, gIdx, cIdx, G, C) ((q)[((bIdx) * (G) * (C)) + ((gIdx) * (C)) + (cIdx)])

/**
 * @brief Computes an approximate of the conditional probability by using a Multinomial approach.
 *
 * Given the observables parameters and the probability matrix, it computes an approximation of `q` with the Multinomial
 * approach.
 *
 * @param[in] *candidates Matrix of dimension (cxb) that stores the results of candidate "c" on ballot box "b".
 * @param[in] *groups Matrix of dimension (bxg) that stores the amount of votes from the demographic group "g".
 * @param[in] *probabilities Matrix of dimension (gxc) with the probabilities of each group and candidate.
 * @param[in] *candidatesVotes Array of size (c) that contains the total votes of each candidate.
 * @param[in] *ballotsVotes Array of size (b) that contains the total amount of votes per ballot.
 * @param[in] *candidatesTotal Integer that contains the total amount of candidates.
 * @param[in] *ballotsTotal Integer that contains the total amount of ballots.
 * @param[in] *groupTotal Integer that contains the total amount of groups.
 *
 * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices, but
 * it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
 *
 */

double *computeQMultinomial(Matrix const *probabilities)
{
    double *array =
        (double *)calloc((size_t)TOTAL_BALLOTS, (size_t)TOTAL_CANDIDATES * (size_t)TOTAL_GROUPS); // Array to return

    // We have to calculate the following: $$\frac{x_{bc}\cdot p_{gc}}{WP_{bc}}\cdot\left(\sum_{\forall
    // c}\frac{WP_{bc}}{x_{bc}\cdot p_{gc}}-\sum_{\forall c}\underbrace{\frac{1}{x_{bc}}}_\text{Suma sobre
    // enteros}\right)$$.

    // -- Summatory calculation for g --
    // This is a simple matrix calculation, to be computed once.
    Matrix WP = createMatrix((int)TOTAL_BALLOTS, (int)TOTAL_GROUPS);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)*ballotsTotal, // Rows in W
                (int)TOTAL_CANDIDATES,                                         // Columns in P
                (int)TOTAL_GROUPS,                                             // Shared dimension
                1.0, W->data, (int)TOTAL_GROUPS, probabilities->data, (int)TOTAL_CANDIDATES, 0.0, WP.data,
                (int)TOTAL_CANDIDATES);

    // --Summatory calculation for c --
    // First term NUMERATOR
    double numSum[TOTAL_BALLOTS];
    colSum(&WP, numSum);

    // First term DENOMINATOR; is x_{cb}^T*p_{gc}^T
    Matrix XP = createMatrix((int)TOTAL_BALLOTS, (int)TOTAL_GROUPS);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, X->rows, // Rows in X
                probabilities->rows,                            // Columns in P
                X->cols,                                        // Shared dimension
                1.0, X->data, X->cols, probabilities->data, probabilities->cols, 0.0, XP.data, probabilities->rows);

    // Second term
    if (inv_BALLOTS_VOTES == NULL || *inv_BALLOTS_VOTES == 0.0)
    {
        inv_BALLOTS_VOTES = (double *)calloc(TOTAL_BALLOTS, sizeof(double));
        for (int b = 0; b < TOTAL_BALLOTS; b++)
        {
            inv_BALLOTS_VOTES[b] = 1.0 / (double)BALLOTS_VOTES[b];
        }
    }

#pragma omp parallel for
    for (int b = 0; b < TOTAL_BALLOTS; b++)
    {
        for (int g = 0; g < TOTAL_GROUPS; g++)
        {
            for (int c = 0; c < TOTAL_CANDIDATES; c++)
            {
                double alpha = (MATRIX_AT_PTR(X, c, b) * MATRIX_AT_PTR(probabilities, g, c)) /
                               (MATRIX_AT(WP, b, c)); // First parenthesis
                double beta = (numSum[b] / MATRIX_AT(XP, b, c)) - inv_BALLOTS_VOTES[b];
                Q_3D(array, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES);
            }
        }
    }
    return array;
}
