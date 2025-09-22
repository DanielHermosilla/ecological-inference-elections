#ifndef MULTIVARIATE_UTILS_H_EIM
#define MULTIVARIATE_UTILS_H_EIM

#ifdef __cplusplus

extern "C"
{
#endif

#include "globals.h"
#include "utils_matrix.h"

    /**
     * @brief Computes the parameters of the unconditional probability
     *
     * Computes the first and second moments of an approximated Multivariate Normal distribution.
     *
     * @param[in] b The index of the ballot box
     * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
     * except the last one. Consider that the information of the last candidate is redundant. The probability matrix
     * could be reduced with the "removeRows()" function from matrixUtils.
     * @param[in, out] *mu An array of size c-1 to store the results of the average.
     * @param[in, out] *sigma A matrix of size (c-1, c-1) to store the sigma matrix.
     *
     * @warning Remember to eliminate one dimension for candidates.
     *
     * @return void. Results to be written on mu and sigma.
     *
     */

    void getParams(EMContext *ctx, int b, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma);

    /**
     * @brief Computes the parameters of the conditional probability
     *
     * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the
     * results of a ballot box.
     *
     * @param[in] b The index of the ballot box
     * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
     * except the last one. Consider that the information of the last candidate is redundant. The probability matrix
     * could be reduced with the "removeCols()" function from matrixUtils.
     * @param[in, out] *conditionalMu A matrix of size (gxc-1) that stores the average of each candidate given a group.
     * @param[in, out] *newSigma An array of matrices of size `g` that stores matrices of size (c-1, c-1) that
     * represents the sigma of each group.
     *
     * @warning Remember to eliminate one dimension for candidates.
     *
     * @return void. Results to be written on mu and sigma.
     *
     */
    void getAverageConditional(EMContext *ctx, int b, const Matrix *probabilitiesReduced, Matrix *conditionalMu,
                               Matrix **conditionalSigma);

    Matrix getBallotPDF(int b, const Matrix *probabilitiesReduced);

    /**
     * @brief Computes the Mahalanobis distance with last candidate adjustment.
     *
     * @param[in] x Pointer to the input feature vector (size C-1).
     * @param[in] mu Pointer to the mean vector (size C-1).
     * @param[in] inverseSigma Pointer to the inverse covariance matrix (size (C-1) x (C-1)).
     * @param[out] maha Pointer to the resulting Mahalanobis distance array (could be C-1 or C)
     * @param[in] size Size of the truncated candidate space (C-1).
     * @param[in] truncate Boolean value to see if the *maha parameter fixes for a lineally dependent sistem and returns
     * a size `C` array.
     */
    void getMahanalobisDist(double *x, double *mu, Matrix *inverseSigma, double *maha, int size, bool reduced);

    /*
     * @brief Obtains the Mahanalobis distance using the Cholesky factorization of Sigma, without having to invert it.
     *
     */
    double getMahanalobisDist2(const Matrix *sigmaL, const double *diff, double *y, double *z, double *ec,
                               double *Sdiag, int n, int need_z, int need_diag);

#ifdef __cplusplus
}
#endif
#endif
