#include "multinomial.h"
#include "globals.h"
#include <R_ext/BLAS.h>

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
 * @param[in] params The optional parameters to the function. For this specific method, there's no supported optional
 * parameters yet.
 *
 * @return A (bxgxc) continuos array with the values of each probability. Understand it as a tensor with matrices, but
 * it's fundamental to be a continuos array in memory for simplificating the posteriors calculations.
 *
 */
double *computeQMultinomial(Matrix const *probabilities, QMethodInput params)
{
    double *array2 =
        (double *)calloc(TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS, sizeof(double)); // Array to return

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
    (&noTranspose, &noTranspose, // transA = 'N', transB = 'N'
     &m, &n, &k,                 // M, N, K
     &alpha, W->data, &m,        // A, LDA = m
     probabilities->data, &k,    // B, LDB = k
     &beta, WP.data, &m,         // C, LDC = m
     (BLAS_INT)1, (BLAS_INT)1    // string lengths for 'N', 'N'
    );

    // ---- Do not parallelize ----
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
                if (denominator != 0) // --- Edge case
                    finalNumerator[c] = numerator / denominator;
                else
                    finalNumerator[c] = 0;
                // ---- Store the value for reusing it later ----
                tempSum += finalNumerator[c];
                // ---...--- //
            }

            for (int c = 0; c < (int)TOTAL_CANDIDATES; c++)
            {                     // ---- For each candidate given a group and a ballot box
                if (tempSum != 0) // --- Edge case
                {
                    Q_3D(array2, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = finalNumerator[c] / tempSum;
                    printf("Adding to q:\t%.4f\n", finalNumerator[c] / tempSum);
                }
                else
                    Q_3D(array2, b, g, c, TOTAL_GROUPS, TOTAL_CANDIDATES) = 0;
                // ---- Store the value ----
            }
        }
    }

    freeMatrix(&WP);
    return array2;
}
