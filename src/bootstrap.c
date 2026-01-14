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
#include "bootstrap.h"
#include "parametric_main.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h> // for R_CheckUserInterrupt()
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

void iterMat(const Matrix *originalX, const Matrix *originalW, Matrix *newX, Matrix *newW, const int *indexArr,
             int indexStart)
{
    // The amount of ballot boxes
    int ballotBoxes = originalW->rows;
    for (int b = 0; b < ballotBoxes; b++)
    {
        int sampledIndex = indexArr[indexStart + b];
        // For the 'w' matrix
        for (int g = 0; g < originalW->cols; g++)
        { // --- For each group given a ballot box
            MATRIX_AT_PTR(newW, b, g) = MATRIX_AT_PTR(originalW, sampledIndex, g);
        }
        // For the 'x' matrix
        for (int c = 0; c < originalX->rows; c++)
        { // --- For each candidate given a ballot box
            MATRIX_AT_PTR(newX, c, b) = MATRIX_AT_PTR(originalX, c, sampledIndex);
        }
    }
}

void iterMatParametric(const Matrix *originalX, const Matrix *originalW, const Matrix *originalV, Matrix *newX,
                       Matrix *newW, Matrix *newV, const int *indexArr, int indexStart)
{
    // The amount of ballot boxes
    int ballotBoxes = originalW->rows;
    for (int b = 0; b < ballotBoxes; b++)
    {
        int sampledIndex = indexArr[indexStart + b];
        // For the 'w' matrix
        for (int g = 0; g < originalW->cols; g++)
        { // --- For each group given a ballot box
            MATRIX_AT_PTR(newW, b, g) = MATRIX_AT_PTR(originalW, sampledIndex, g);
        }
        // For the 'x' matrix
        for (int c = 0; c < originalX->cols; c++)
        { // --- For each candidate given a ballot box
            MATRIX_AT_PTR(newX, b, c) = MATRIX_AT_PTR(originalX, sampledIndex, c);
        }
        // For the 'v' matrix
        for (int a = 0; a < originalV->cols; a++)
        { // --- For each attribute given a ballot box
            MATRIX_AT_PTR(newV, b, a) = MATRIX_AT_PTR(originalV, sampledIndex, a);
        }
    }
}

Matrix standardDeviations(Matrix *bootstrapResults, Matrix *sumMatrix, int totalIter)
{

    // Get the mean for each component
    for (int i = 0; i < sumMatrix->rows; i++)
    {
        for (int j = 0; j < sumMatrix->cols; j++)
        {
            MATRIX_AT_PTR(sumMatrix, i, j) /= totalIter;
        }
    }
    // printMatrix(sumMatrix);

    Matrix sdMatrix = createMatrix(sumMatrix->rows, sumMatrix->cols);

    // Get the summatory (x_i - \mu)^2
    for (int h = 0; h < totalIter; h++)
    {
        // Yields the summatory for each dimension
        for (int i = 0; i < sdMatrix.rows; i++)
        {
            for (int j = 0; j < sdMatrix.cols; j++)
            {
                double diff = MATRIX_AT(bootstrapResults[h], i, j) - MATRIX_AT_PTR(sumMatrix, i, j);
                MATRIX_AT(sdMatrix, i, j) += diff * diff;
            }
        }
        // freeMatrix(&bootstrapResults[h]);
    }

    // Make the division and get the square root
    for (int i = 0; i < sdMatrix.rows; i++)
    {
        for (int j = 0; j < sdMatrix.cols; j++)
        {
            double val = sqrt(MATRIX_AT(sdMatrix, i, j) / (totalIter - 1));
            MATRIX_AT(sdMatrix, i, j) = val == 0 ? NAN : val;
        }
    }
    return sdMatrix;
}

void bootstrapParametric(const Matrix *X, const Matrix *W, const Matrix *V, Matrix *beta, Matrix *alpha, int bootiter,
                         const int maxiter, const double maxtime, const double ll_threshold, const int maxnewton,
                         const bool verbose, Matrix *sdBeta, Matrix *sdAlpha,
                         const char *adjust_prob_cond_method, bool adjust_prob_cond_every)
{
    // ---- Initial variables
    int bdim = W->rows;
    int samples = bdim * bootiter;

    // ---- Border case
    if (bdim == 1)
    {
        Matrix infBetas = createMatrix(beta->rows, beta->cols);
        Matrix infAlpha = createMatrix(alpha->rows, alpha->cols);
        fillMatrix(&infBetas, 9999);
        fillMatrix(&infAlpha, 9999);
        *sdBeta = infBetas;
        *sdAlpha = infAlpha;
        return;
    }

    // ---- Generate the indices for bootstrap ---- //
    int *indices = Calloc(bdim * bootiter, int);
    // For each bootstrap replicate i
sampling:
    for (int i = 0; i < bdim * bootiter; i++)
    {
        indices[i] = (int)(unif_rand() * bdim);
    }
    // Check that every index is not the same
    for (int i = 1; i < bdim * bootiter; i++)
    {
        if (indices[i] != indices[i - 1])
            break;
        if (i == bdim * bootiter - 1)
        {
            goto sampling;
        }
    }
    // We want to avoid the case where the same ballot box is drawn FOR EACH placement
    // This has a probability of 1/b^b. Maybe this calculation could be avoided at 6 > ballot boxes,
    // since then it becomes practically 0
    // ---...--- //

    // ---- Execute the bootstrap algorithm ---- //
    Matrix sumMatBeta = createMatrix(beta->rows, beta->cols);
    Matrix *resultsBeta = Calloc(bootiter, Matrix);
    Matrix sumMatAlpha = createMatrix(alpha->rows, alpha->cols);
    Matrix *resultsAlpha = Calloc(bootiter, Matrix);

    for (int i = 0; i < bootiter; i++)
    {
        if (verbose && bootiter > 20 && (i % (bootiter / 20) == 0)) // Print every 5% (20 intervals)
        {
            double progress = (double)i / bootiter * 100;
            Rprintf("%.0f%% of iterations completed.\n", progress);
        }
        // ---- Declare variables for the current iteration
        Matrix iterX = createMatrix(X->rows, X->cols);
        Matrix iterW = createMatrix(W->rows, W->cols);
        Matrix iterV = createMatrix(V->rows, V->cols);
        Matrix BetaR = copMatrix(beta);
        Matrix AlphaR = copMatrix(alpha);

        iterMatParametric(X, W, V, &iterX, &iterW, &iterV, indices, i * bdim);

        // Declare EM variables, they're not used in this case...
        // It could be useful to yield a mean if the user wants to (logLL mean?)
        double elapsed = 0.0;
        int total_iter = 0;
        double log_likelihood = 0.0;
        Matrix *finalProb =
            EM_Algorithm(&iterX, &iterW, &iterV, &BetaR, &AlphaR, maxiter, maxtime, ll_threshold, maxnewton, false,
                         &elapsed, &total_iter, &log_likelihood, adjust_prob_cond_method, adjust_prob_cond_every);

        // Sum each value so later we can get the mean
        for (int g = 0; g < BetaR.rows; g++)
        {
            for (int c = 0; c < BetaR.cols; c++)
            {
                MATRIX_AT(sumMatBeta, g, c) += MATRIX_AT(BetaR, g, c);
            }
        }
        resultsBeta[i] = BetaR;

        for (int c = 0; c < AlphaR.rows; c++)
        {
            for (int a = 0; a < AlphaR.cols; a++)
            {
                MATRIX_AT(sumMatAlpha, c, a) += MATRIX_AT(AlphaR, c, a);
            }
        }
        resultsAlpha[i] = AlphaR;

        // ---- Release loop allocated variables ---- //
        if (finalProb != NULL)
        {
            for (int b = 0; b < iterV.rows; b++)
            {
                freeMatrix(&finalProb[b]);
            }
            Free(finalProb);
        }
        freeMatrix(&iterX);
        freeMatrix(&iterW);
        freeMatrix(&iterV);
        // ---...--- //
    }

    *sdBeta = standardDeviations(resultsBeta, &sumMatBeta, bootiter);
    *sdAlpha = standardDeviations(resultsAlpha, &sumMatAlpha, bootiter);

    if (verbose)
    {
        Rprintf("Bootstrapping finished!\nThe estimated standard deviation matrix for beta is:\n");
        printMatrix(sdBeta);
        Rprintf("The estimated standard deviation matrix for alpha is:\n");
        printMatrix(sdAlpha);
    }

    Free(indices);
    freeMatrix(&sumMatBeta);
    freeMatrix(&sumMatAlpha);

    for (int i = 0; i < bootiter; i++)
    {
        freeMatrix(&resultsBeta[i]);
        freeMatrix(&resultsAlpha[i]);
    }
    Free(resultsBeta);
    Free(resultsAlpha);
}
/**
 *  Returns an array of col-major matrices with bootstrapped matrices.
 *
 * @param[in] xmat The original X array
 * @param[in] wmat The original W array
 * @param[in] bootiter The amount of iterations for bootstrapping
 * @param[in] p_method The method for obtaining the initial probability
 * @param[in] q_method Pointer to a string that indicates the method or calculating "q". Currently it supports "Hit
 * and Run", "mult", "mvn_cdf", "mvn_pdf" and "exact" methods.
 * @param[in] convergence Threshold value for convergence. Usually it's set to 0.001.
 * @param[in] maxIter Integer with a threshold of maximum iterations. Usually it's set to 100.
 * @param[in] maxSeconds Double with the value of the maximum amount of seconds to use.
 * @param[in] verbose Wether to verbose useful outputs.
 * @param[in, out] time The time that the algorithm took.
 * @param[in, out] iterTotal Total amount of iterations.
 * @param[in, out] logLLarr The loglikelihood array
 * @param[in, out] finishing_reason The reason that the algorithm has been stopped. It can either be 0, 1, 2, 3,
 * representing a normal convergence, log likelihood decrease, maximum time reached and maximum iterations reached,
 * respectively.
 * @param[out] avgProbOut Optional pointer to store the average probability matrix from the bootstrap iterations.
 *
 *
 * @return Matrix with the bootstrap standard deviations for each probability component.
 */
Matrix bootstrapA(const Matrix *xmat, const Matrix *wmat, int bootiter, const char *q_method, const char *p_method,
                  const double convergence, const double log_convergence, const int maxIter, const double maxSeconds,
                  const bool verbose, Matrix *probMatrix, QMethodInput *inputParams, Matrix *avgProbOut)
{

    // ---- Initial variables
    int bdim = wmat->rows;
    int samples = bdim * bootiter;
    int matsize = wmat->cols * xmat->rows;

    // ---- Border case
    if (bdim == 1)
    {
        Matrix infMat = createMatrix(wmat->cols, xmat->rows);
        fillMatrix(&infMat, 9999);
        if (avgProbOut != NULL)
        {
            Matrix avgMat = createMatrix(wmat->cols, xmat->rows);
            fillMatrix(&avgMat, 9999);
            *avgProbOut = avgMat;
        }
        return infMat;
    }
    // ---- Generate the indices for bootstrap ---- //
    int *indices = Calloc(bdim * bootiter, int);
    // For each bootstrap replicate i
sampling:
    for (int i = 0; i < bdim * bootiter; i++)
    {
        indices[i] = (int)(unif_rand() * bdim);
    }
    // Check that every index is not the same
    for (int i = 1; i < bdim * bootiter; i++)
    {
        if (indices[i] != indices[i - 1])
            break;
        if (i == bdim * bootiter - 1)
        {
            goto sampling;
        }
    }
    // We want to avoid the case where the same ballot box is drawn FOR EACH placement
    // This has a probability of 1/b^b. Maybe this calculation could be avoided at 6 > ballot boxes,
    // since then it becomes practically 0
    // ---...--- //

    // ---- Execute the bootstrap algorithm ---- //
    Matrix sumMat = createMatrix(wmat->cols, xmat->rows);
    Matrix *results = Calloc(bootiter, Matrix);
    for (int i = 0; i < bootiter; i++)
    {
        if (verbose && bootiter > 20 && (i % (bootiter / 20) == 0)) // Print every 5% (20 intervals)
        {
            double progress = (double)i / bootiter * 100;
            Rprintf("%.0f%% of iterations completed.\n", progress);
        }
        // ---- Declare variables for the current iteration
        Matrix iterX = createMatrix(xmat->rows, xmat->cols);
        Matrix iterW = createMatrix(wmat->rows, wmat->cols);
        iterMat(xmat, wmat, &iterX, &iterW, indices, i * bdim);
        // setParameters(&iterX, &iterW);
        // Matrix iterP = getInitialP(p_method);

        // Declare EM variables, they're not used in this case...
        // It could be useful to yield a mean if the user wants to (logLL mean?)
        double time;
        double logLLarr = 0;
        int finishing_reason, iterTotal;
        EMContext *ctx =
            EMAlgoritm(&iterX, &iterW, p_method, q_method, convergence, log_convergence, maxIter, maxSeconds, false,
                       &time, &iterTotal, &logLLarr, &finishing_reason, probMatrix, inputParams);
        Matrix *resultP = &ctx->probabilities;
        Matrix copyP = createMatrix(resultP->rows, resultP->cols);
        memcpy(copyP.data, resultP->data, sizeof(double) * resultP->rows * resultP->cols);
        // Sum each value so later we can get the mean
        for (int j = 0; j < wmat->cols; j++)
        {
            for (int k = 0; k < xmat->rows; k++)
            {
                MATRIX_AT(sumMat, j, k) += MATRIX_AT_PTR(resultP, j, k);
            }
        }

        results[i] = copyP;
        // memcpy(&results[i * matsize], resultP.data, matsize * sizeof(double));

        cleanup(ctx); // Cleanup the context, it frees the memory allocated for the matrices
        // ---- Release loop allocated variables ---- //
        // freeMatrix(&iterP);
        /*
        cleanup();
        if (strcmp(q_method, "exact") == 0)
        {
            cleanExact();
        }
        else if (strcmp(q_method, "mcmc") == 0)
        {
            cleanHitAndRun();
        }
         */
        // else if (strcmp(q_method, "mult"))
        //{
        //    cleanMultinomial();
        //}
        freeMatrix(&iterX);
        freeMatrix(&iterW);
        // ---...--- //
    }
    Matrix sdReturn = standardDeviations(results, &sumMat, bootiter);
    if (avgProbOut != NULL)
    {
        *avgProbOut = sumMat;
    }
    if (verbose)
    {
        Rprintf("Bootstrapping finished!\nThe estimated standard deviation matrix (g x c) is:\n");
        printMatrix(&sdReturn);
    }

    Free(indices);
    if (avgProbOut == NULL)
    {
        freeMatrix(&sumMat);
    }
    for (int i = 0; i < bootiter; i++)
    {

        freeMatrix(&results[i]);
    }
    Free(results);

    return sdReturn;
}

Matrix bootSingleMat(Matrix *xmat, Matrix *wmat, int bootiter, const bool verbose)
{

    // ---- Initial variables
    int bdim = wmat->rows;
    int samples = bdim * bootiter;
    int matsize = xmat->rows;
    // ---- Generate the indices for bootstrap ---- //
    int *indices = Calloc(bdim * bootiter, int);
    for (int j = 0; j < samples; j++)
    {
        indices[j] = (int)(unif_rand() * bdim);
    }
    // ---...--- //

    // ---- Execute the bootstrap algorithm ---- //
    Matrix sumMat = createMatrix(1, xmat->rows);
    Matrix *results = Calloc(bootiter, Matrix);
    for (int i = 0; i < bootiter; i++)
    {
        if (verbose && (i % (bootiter / 20) == 0)) // Print every 5% (20 intervals)
        {
            double progress = (double)i / bootiter * 100;
            Rprintf("An %.0f%% of iterations have been done.\n", progress);
        }
        // ---- Declare variables for the current iteration
        Matrix iterX = createMatrix(xmat->rows, xmat->cols);
        Matrix iterW = createMatrix(wmat->rows, 1);
        // iterMat(xmat, wmat, &iterX, &iterW, indices, i * bdim);
        // setParameters(&iterX, &iterW);
        EMContext *ctxTemp = createEMContext(xmat, wmat, "group_proportional", (QMethodInput){0});

        Matrix *resultP = &ctxTemp->probabilities;
        Matrix copyP = createMatrix(resultP->rows, resultP->cols);
        memcpy(copyP.data, resultP->data, sizeof(double) * resultP->rows * resultP->cols);
        // Sum each value so later we can get the mean
        for (int k = 0; k < xmat->rows; k++)
        {
            MATRIX_AT_PTR(resultP, 0, k) = (double)ctxTemp->candidates_votes[k] / (double)TOTAL_VOTES;
            MATRIX_AT(sumMat, 0, k) += MATRIX_AT_PTR(resultP, 0, k);
        }

        results[i] = copyP;

        // memcpy(&results[i * matsize], resultP.data, matsize * sizeof(double));

        // ---- Release loop allocated variables ---- //
        cleanup(ctxTemp);
        freeMatrix(&iterX);
        freeMatrix(&iterW);
        // ---...--- //
    }
    Matrix sdReturn = standardDeviations(results, &sumMat, bootiter);
    if (verbose)
    {
        Rprintf("Bootstrapping finished!\nThe estimated standard deviation matrix (g x c) is:\n");
        printMatrix(&sdReturn);
    }

    Free(indices);
    freeMatrix(&sumMat);
    for (int i = 0; i < bootiter; i++)
    {
        freeMatrix(&results[i]);
    }
    Free(results);

    return sdReturn;
}
