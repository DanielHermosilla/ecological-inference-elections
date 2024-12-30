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
 * @param[in] *candidatesTotal Integer that contains the total amount of candidates.
 * @param[in] *ballotsTotal Integer that contains the total amount of ballots.
 * @param[in] *groupTotal Integer that contains the total amount of groups.
 *
 * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices, but
 * it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
 *
 */

double *computeQMultinomial(Matrix const *candidates, Matrix const *groups, Matrix const *probabilities,
                            uint32_t const *candidatesVotes, uint16_t const *candidatesTotal,
                            uint16_t const *ballotsTotal, uint16_t const *groupsTotal)
{

    // -- Summatory calculation for g calculation --
    // This is a simple matrix calculation, to be computed once.
    Matrix WP = createMatrix((int)*ballotsTotal, (int)*groupsTotal);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)*ballotsTotal, // Rows in W
                (int)*candidatesTotal,                                         // Columns in P
                (int)*groupsTotal,                                             // Shared dimension
                1.0, candidates->data, (int)*groupsTotal, probabilities->data, (int)*candidatesTotal, 0.0, WP.data,
                (int)*candidatesTotal);

    // -- Summatory calculations for c--
    // I would like to calculate by the following approach: $\sum_{\forall c}\frac{WP_{bc}}{p_{gc}}=WP\times (\text{inv}
    // P^T)$
    Matrix invPT = createMatrix((int)*candidatesTotal, (int)*groupsTotal); // The inverse of P and transposed.
    for (int g = 0; g < *groupsTotal; g++)
    {
        for (int c = 0; c < *candidatesTotal; c++)
        {
            MATRIX_AT(invPT, c, g) = 1.0 / MATRIX_AT_PTR(probabilities, g, c);
        }
    }
    Matrix S = createMatrix((int)*ballotsTotal, (int)*groupsTotal); // Let S be the summatory for all c'
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                (int)*ballotsTotal,    // M = b
                (int)*groupsTotal,     // N = g
                (int)*candidatesTotal, // K = c
                1.0,
                WP.data,               // (b x c)
                (int)*candidatesTotal, // leading dimension = c
                invPT.data,            // (c x g)
                (int)*groupsTotal,     // leading dimension = g
                0.0,
                S.data, // (b x g)
                (int)*groupsTotal);
    // -- End summatory calculations --

    double *array = (double *)calloc((size_t)ballotsTotal * (size_t)candidatesTotal * (size_t)groupsTotal,
                                     sizeof(double)); // Array to return

    // Let:
    // -alpha = $\frac{x_{bc}\cdot p_{gc}}{\hat{x_b}}$
    // -beta = $S_{bg}-\lvert C\rvert$
    // -gamma = $\frac{1}{(W\cdot P)_{bc}-p_{gc}}\right$

#pragma omp parallel for
    for (int b = 0; b < *ballotsTotal; b++)
    {
        for (int g = 0; g < *groupsTotal; g++)
        {
            double beta = MATRIX_AT(S, b, g) - *candidatesTotal;
            for (int c = 0; c < *candidatesTotal; c++)
            {
                double alpha =
                    MATRIX_AT_PTR(candidates, b, c) * MATRIX_AT_PTR(probabilities, g, c) / candidatesVotes[b];
                double gamma = 1.0 / (MATRIX_AT(WP, b, c) - MATRIX_AT_PTR(probabilities, g, c));
                Q_3D(array, b, g, c, *groupsTotal, *candidatesTotal) = beta * alpha * gamma;
            }
        }
    }
    return array;
}

int main()
{
    return 1;
}
