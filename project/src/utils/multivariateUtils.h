#include "globals.h"
#include "matrixUtils.h"

#ifndef MULTIVARIATE_UTILS_H
#define MULTIVARIATE_UTILS_H

/**
 * @brief Computes the parameters of the unconditional probability
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution.
 *
 * @param[in] b The index of the ballot box
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeRows()" function from matrixUtils.
 * @param[in, out] *mu An array of size c-1 to store the results of the average.
 * @param[in, out] *sigma A matrix of size (c-1, c-1) to store the sigma matrix.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getParams(int b, const Matrix *probabilitiesReduced, double *mu, Matrix *sigma);

/**
 * @brief Computes the parameters of the unconditional probability WITH the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesNotReduced Matrix of dimension (gxc) with the probabilities of each group and candidate,
 * @param[in, out] *mu An array of size c to store the results of the average.
 * @param[in, out] *sigma A matrix of size (c, c) to store the sigma matrix.
 *
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getParamsFull(int b, const Matrix *probabilitiesNotReduced, double *mu, Matrix *sigma);

/**
 * @brief Computes the parameters of the conditional probability
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the results
 * of a ballot box.
 *
 * @param[in] b The index of the ballot box
 * @param[in] *probabilitiesReduce Matrix of dimension (gxc-1) with the probabilities of each group and candidate,
 * except the last one. Consider that the information of the last candidate is redundant. The probability matrix could
 * be reduced with the "removeCols()" function from matrixUtils.
 * @param[in, out] *conditionalMu A matrix of size (gxc-1) that stores the average of each candidate given a group.
 * @param[in, out] *newSigma An array of matrices of size `g` that stores matrices of size (c-1, c-1) that represents
 * the sigma of each group.
 *
 * @warning Remember to eliminate one dimension for candidates.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */
void getAverageConditional(int b, const Matrix *probabilitiesReduced, Matrix *conditionalMu, Matrix **conditionalSigma);

/**
 * @brief Computes the parameters of the conditional probability WITH the last candidate
 *
 * Computes the first and second moments of an approximated Multivariate Normal distribution conditional to the results
 * of a group.
 *
 * @param[in] b The index of the ballot box
 * @param[in] g The index of the group
 * @param[in] *probabilitiesNotReduced Matrix of dimension (gxc) with the probabilities of each group and candidate,
 * @param[in, out] *newMu An array of size c to store the results of the average.
 * @param[in, out] *newSigma A matrix of size (c, c) to store the sigma matrix.
 *
 * @return void. Results to be written on mu and sigma.
 *
 */

void getAverageConditionalFull(int b, const Matrix *probabilitiesNotReduced, Matrix *conditionalMu,
                               Matrix **conditionalSigma);

#endif
