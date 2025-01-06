#include "multinomial.h"
#include <cblas.h>

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

/**
 * @brief Computes an approximate of the conditional probability by using a Multinomial approach.
 *
 * Given the observables parameters and the probability matrix, it computes an approximation of `q` with the Multinomial
 * approach.
 *
 * @param[in] *probabilities Matrix of dimension (gxc) with the probabilities of each group and candidate.
 *
 * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices, but
 * it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
 *
 */

double *computeQMultinomial(Matrix const *probabilities)
{
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return

    // -- Summatory calculation for g --
    // This is a simple matrix calculation, to be computed once.
    Matrix WP = createMatrix((int)TOTAL_BALLOTS, (int)TOTAL_CANDIDATES);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, (int)TOTAL_BALLOTS, // Rows in W
                (int)TOTAL_CANDIDATES,                                         // Columns in P
                (int)TOTAL_GROUPS,                                             // Shared dimension
                1.0, W->data, (int)TOTAL_GROUPS, probabilities->data, (int)TOTAL_CANDIDATES, 0.0, WP.data,
                (int)TOTAL_CANDIDATES);

    Matrix temp = createMatrix((int)TOTAL_BALLOTS, (int)TOTAL_GROUPS);
#pragma omp parallel for
    for (int b = 0; b < (int)TOTAL_BALLOTS; b++)
    {
        for (int g = 0; g < (int)TOTAL_GROUPS; g++)
        {
            for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            {

                double result = (1.0 / computeR(probabilities, &WP, b, c, g)) * MATRIX_AT_PTR(probabilities, g, c) *
                                MATRIX_AT_PTR(X, c, b);
                MATRIX_AT(temp, b, g) += result;
            }
            for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            {
                double result2 = (1.0 / computeR(probabilities, &WP, b, c, g)) * MATRIX_AT_PTR(probabilities, g, c) *
                                 MATRIX_AT_PTR(X, c, b);
                Q_3D(array2, b, g, c, (int)TOTAL_GROUPS, (int)TOTAL_CANDIDATES) = result2 / MATRIX_AT(temp, b, g);
            }
        }
    }

    freeMatrix(&WP);
    freeMatrix(&temp);
    return array2;
}
