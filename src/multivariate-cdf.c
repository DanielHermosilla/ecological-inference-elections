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

#ifdef _OPENMP
#include <omp.h>
#endif
#include "multivariate-cdf.h"
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/RS.h> /* for R_Calloc/R_Free, F77_CALL */
#include <R_ext/Random.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

// ---- Arena with reusable buffers ---- //
typedef struct
{
    double *feature;    // size C
    double *mu_row;     // size C-1
    double *featureA;   // size C-1
    double *featureB;   // size C-1
    double *mc_results; // size C
} ArenaCDF;

/**
 * @brief Create and allocate an ArenaCDF for dimension C.
 */
static ArenaCDF ArenaCDF_init(int C)
{
    ArenaCDF A = (ArenaCDF){0};
    A.feature = (double *)Calloc(C, double);
    A.mu_row = (double *)Calloc(C - 1, double);
    A.featureA = (double *)Calloc(C - 1, double);
    A.featureB = (double *)Calloc(C - 1, double);
    A.mc_results = (double *)Calloc(C, double);
    return A;
}

/**
 * @brief Free the buffers of the ArenaCDF.
 */
static void ArenaCDF_free(ArenaCDF *A)
{
    if (!A)
        return;
    if (A->feature)
        Free(A->feature);
    if (A->mu_row)
        Free(A->mu_row);
    if (A->featureA)
        Free(A->featureA);
    if (A->featureB)
        Free(A->featureB);
    if (A->mc_results)
        Free(A->mc_results);
    memset(A, 0, sizeof(*A));
}

void allocateSeed(EMContext *ctx, QMethodInput params)
{
    GetRNGstate();
    int total = (ctx->C - 1) * params.monteCarloIter;

    ctx->cdf_seeds = (double *)Calloc(total, double);
    for (int i = 0; i < total; i++)
    {
        ctx->cdf_seeds[i] = unif_rand();
    }
    PutRNGstate();
}

// ---- Helpers without dynamic allocation (avoid Calloc/Free in loops) ---- //

/**
 * @brief Copy column `col` of M into a preallocated buffer.
 *
 * @param[in]  *M   Input matrix
 * @param[in]  col  Column index
 * @param[out] *out Buffer of size M->rows
 */
static inline void getColumn_into(const Matrix *M, int col, double *out)
{
    for (int r = 0; r < M->rows; r++)
        out[r] = MATRIX_AT_PTR(M, r, col);
}

/**
 * @brief Copy row `row` of M into a preallocated buffer.
 *
 * @param[in]  *M   Input matrix
 * @param[in]  row  Row index
 * @param[out] *out Buffer of size M->cols
 */
static inline void getRow_into(const Matrix *M, int row, double *out)
{
    for (int c = 0; c < M->cols; c++)
        out[c] = MATRIX_AT_PTR(M, row, c);
}

static double pdf_midpoint(const Matrix *cholesky, const double *mu, const double *lower, const double *upper, int d)
{
    // midpoint
    double mid[d];
    for (int i = 0; i < d; i++)
        mid[i] = 0.5 * (lower[i] + upper[i]);

    // diff = mid - mu
    double diff[d];
    for (int i = 0; i < d; i++)
        diff[i] = mid[i] - mu[i];

    // Solve L z = diff  (L lower-triangular Cholesky s.t. L L^T = Σ)
    double z[d];
    for (int i = 0; i < d; i++)
    {
        double lii = MATRIX_AT_PTR(cholesky, i, i);
        if (!(lii > 0.0) || !isfinite(lii))
            return DBL_MIN; // ill-conditioned, avoid exact zero
        double s = diff[i];
        for (int j = 0; j < i; j++)
            s -= MATRIX_AT_PTR(cholesky, i, j) * z[j];
        z[i] = s / lii;
    }

    // z^T z
    double z2 = 0.0;
    for (int i = 0; i < d; i++)
        z2 += z[i] * z[i];

    // log|\Sigma|^{1/2} = sum(log L_ii)
    double log_det_sqrt = 0.0;
    for (int i = 0; i < d; i++)
    {
        double lii = MATRIX_AT_PTR(cholesky, i, i);
        log_det_sqrt += log(lii);
    }

    // log pdf = -0.5*d*log(2\pi) - log|\Sigma|^{1/2} - 0.5*z^T z
    double log_pdf = -0.5 * d * M_LN_2PI - log_det_sqrt - 0.5 * z2;

    // If caller needs a number in (0, ∞), clamp to avoid exact 0
    if (log_pdf < log(DBL_MIN))
        return DBL_MIN; // underflow guard
    if (log_pdf > log(DBL_MAX))
        return DBL_MAX;        // overflow guard (unlikely)
    double pdf = exp(log_pdf); // exponentiate for getting the final result
    return (isfinite(pdf) && pdf > 0.0) ? pdf : DBL_MIN;
}
/*
 * @brief Compute the Montecarlo approximation proposed by Alan genz towards the
 * most recent method, using univariate conditional with quasirandom numbers.
 *
 * Computes an heuristic from the first Multivariate CDF proposed. More details
 * at the following paper:
 * https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities
 *
 * Note that this method gives the option to impose an error threshold. So, it
 * stops iterating if, either the maximum iterations are accomplished or the
 * error threshold is passed.
 *
 * @param[in] *cholesky. The cholesky matrix of the current iteration.
 * @param[in] *lowerBounds. An array with the initial lower bounds.
 * @param[in] *upperBounds. An array with the initial upper bounds.
 * @param[in] epsilon. The minimum error threshold accepted.
 * @param[in] iterations. The maximum iterations
 * @param[in] mvnDim. The dimension of the multivariate normal.
 *
 * @return The value of the estimated integral
 */
double genzMontecarloNew(EMContext *ctx, const Matrix *cholesky, const double *lowerBounds, const double *upperBounds,
                         double epsilon, int iterations, int mvnDim)
{
    // GetRNGstate();

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    double mean = 0;
    int currentIterations = 1;
    double a[mvnDim], b[mvnDim], y[mvnDim], currentError;
    a[0] = pnorm(lowerBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 0.0, 1.0, 1, 0);
    b[0] = pnorm(upperBounds[0] / MATRIX_AT_PTR(cholesky, 0, 0), 0.0, 1.0, 1, 0);
    // ---...--- //

    do
    {
        // ---- Generate a pseudoRandom number for each iteration
        double pseudoRandom[mvnDim];
        for (int i = 0; i < mvnDim; i++)
        {
            /* pseudoRandom[i] = unif_rand(); */
            // pseudoRandom[i] = ctx->cdf_seeds[currentIterations * (mvnDim - 1) + i];
            int mvnDim = ctx->C - 1;
            pseudoRandom[i] = ctx->cdf_seeds[currentIterations * mvnDim + i];
            /* start currentIter = 0 and ensure allocateSeed() actually ran */
        }

        // ---- Compute the base case
        y[0] = qnorm(a[0] + pseudoRandom[0] * (b[0] - a[0]), 0.0, 1.0, 1, 0);
        double summatory;
        double P = b[0] - a[0];

        // ---- Do the main loop ---- //
        for (int i = 1; i < mvnDim; i++)
        {
            // ---- Note that the summatory is equivalent to
            // $\sum_{j=1}^{i-1}c_{ij}*y_{j}$.
            summatory = 0;
            for (int j = 0; j < i; j++)
            {
                summatory += MATRIX_AT_PTR(cholesky, i, j) * y[j];
            }
            a[i] = pnorm((lowerBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 0.0, 1.0, 1, 0);
            b[i] = pnorm((upperBounds[i] - summatory) / MATRIX_AT_PTR(cholesky, i, i), 0.0, 1.0, 1, 0);
            double difference = b[i] - a[i];
            y[i] = qnorm(a[i] + pseudoRandom[i] * difference, 0.0, 1.0, 1, 0);
            P *= difference;
        }
        // ---...--- //
        // ---- Get the stopping parameters ---- //
        intsum += P;
        currentIterations += 1;
        mean = intsum / currentIterations;
        varsum += pow(P - mean, 2);
        currentError = sqrt(varsum / (currentIterations * (currentIterations - 1)));
        // ---...--- //
    } while (currentError > epsilon && currentIterations < iterations);

    /* PutRNGstate(); */
    return mean;
}

// All of the conventions used are from genz paper
/*
 * @brief Compute the Montecarlo approximation proposed by Alan genz
 *
 * Compute a `manual` Montecarlo simulation, proposed by Dr. Alan genz book
 * "Computation of Multivariate Normal and t Probabilities". All of the variable
 * names and conventions are adopted towards his proposed algorithm aswell.
 * Refer to
 * https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities
 * for more information.
 *
 * Note that this method gives the option to impose an error threshold. So, it
 * stops iterating if, either the maximum iterations are accomplished or the
 * error threshold is passed.
 *
 * @param[in] *cholesky. The cholesky matrix of the current iteration.
 * @param[in] *lowerBounds. An array with the initial lower bounds.
 * @param[in] *upperBounds. An array with the initial upper bounds.
 * @param[in] epsilon. The minimum error threshold accepted.
 * @param[in] iterations. The maximum iterations
 * @param[in] mvnDim. The dimension of the multivariate normal.
 *
 * @return The value of the estimated integral
 */
double genzMontecarlo(EMContext *ctx, const Matrix *cholesky, const double *lowerBounds, const double *upperBounds,
                      double epsilon, int iterations, int mvnDim)
{

    // ---- Initialize randomizer ---- //
    GetRNGstate();
    // Rprintf("Using the cholesky matrix of:\n");
    // printMatrix(cholesky);

    // ---- Initialize Montecarlo variables ---- //
    double intsum = 0;
    double varsum = 0;
    int currentIterations = 0;
    double currentError = DBL_MAX;
    double d[mvnDim], e[mvnDim], f[mvnDim];
    double diag0 = MATRIX_AT_PTR(cholesky, 0, 0);
    if (diag0 == 0)
        return 0;
    d[0] = pnorm(lowerBounds[0] / diag0, 0.0, 1.0, 1, 0);
    e[0] = pnorm(upperBounds[0] / diag0, 0.0, 1.0, 1, 0);
    f[0] = (e[0] - d[0]);
    //  ---...--- //
    do
    {
        // ---- Initialize the loop variables ---- //
        d[0] = pnorm(lowerBounds[0] / diag0, 0.0, 1.0, 1, 0);
        e[0] = pnorm(upperBounds[0] / diag0, 0.0, 1.0, 1, 0);
        f[0] = (e[0] - d[0]);

        double randomVector[mvnDim - 1];
        for (int i = 0; i < mvnDim - 1; i++)
        // ---...--- //

        // ---- Generate the random vector ---- //
        { // --- For each dimension
            // ---- Generate random values in [0,1) ----
            // randomVector[i] = unif_rand();
            randomVector[i] = ctx->cdf_seeds[currentIterations * (mvnDim - 1) + i];
            // ---...--- //
        }
        // ---...--- //

        // ---- Calculate the integral with their new bounds ---- //
        double summatory;
        double y[mvnDim];
        for (int i = 1; i < mvnDim; i++)
        {
            double draw = d[i - 1] + randomVector[i - 1] * (e[i - 1] - d[i - 1]);
            y[i - 1] = qnorm(draw, 0.0, 1.0, 1, 0);
            // ---- Note that the summatory is equivalent to
            // $\sum_{j=1}^{i-1}c_{ij}*y_{j}$.

            summatory = 0;
            for (int j = 0; j < i; j++)
            {
                summatory += MATRIX_AT_PTR(cholesky, i, j) * y[j];
            }

            double diagii = MATRIX_AT_PTR(cholesky, i, i);
            d[i] = diagii != 0 ? pnorm((lowerBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0) : 0;
            e[i] = diagii != 0 ? pnorm((upperBounds[i] - summatory) / diagii, 0.0, 1.0, 1, 0) : 0;
            f[i] = (e[i] - d[i]) * f[i - 1];
        }
        // ---...--- //

        // ---- Compute the final indicators from the current loop ---- //
        intsum += f[mvnDim - 1];
        varsum += pow(f[mvnDim - 1], 2);
        currentIterations += 1;
        if (currentIterations < 2)
            continue;
        double mean = intsum / currentIterations;
        double variance = (varsum / currentIterations) - pow(mean, 2);
        // ---...--- //
        currentError = sqrt(variance / currentIterations);

    } while ((currentError > epsilon && currentIterations < iterations));
    // Rprintf("Iterations: %d, Error: %.10f\n", currentIterations, currentError);
    // ---...--- //
    PutRNGstate();

    return intsum / currentIterations;
}

/*
 * @brief Calls the main function to start all of the Montecarlo process
 *
 * Calls the function with all of the parameters needed to get the Montecarlo
 * simulation of the Multivariate CDF.
 *
 * @param[in] *chol A matrix with the cholenksy values of the current group
 * @param[in] *mu An array with the average values of the current feature vector
 * @param[in] *lowerLimits An array with the lower bounds of the integral
 * (defined by the hypercube)
 * @param[in] *upperLimits An array with the upper bounds of the integral
 * (defined by the hypercube)
 * @param[in] mvnDim The dimensions of the multivariate normal. Usually it's
 * C-1.
 * @param[in] maxSamples Amount of samples for the Montecarlo simulation.
 * @param[in] epsilon The error threshold used for the genz Montecarlo.
 * @param[in] *method The method for calculating the Montecarlo simulation.
 * Currently available methods are `Plain`, `Miser` and `Vegas`.
 *
 * @note Refer to https://www.gnu.org/software/gsl/doc/html/montecarlo.html
 *
 * @return The result of the approximated integral
 */

double Montecarlo(EMContext *ctx, Matrix *chol, double *lowerLimits, double *upperLimits, int mvnDim, int maxSamples,
                  double epsilon, const char *method)
{
    // ---- Set up the initial parameter ---- //
    double result;
    // ---...--- //
    // ---- Case where there are no values ---- //
    if (MATRIX_AT_PTR(chol, 0, 0) == 0)
    {
        return 0;
    }

    // ---- Perform integration ---- //
    if (strcmp(method, "genz") == 0)
    {
        result = genzMontecarlo(ctx, chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else if (strcmp(method, "genz2") == 0)
    {
        result = genzMontecarloNew(ctx, chol, lowerLimits, upperLimits, epsilon, maxSamples, mvnDim);
        return result;
    }
    else
    {
        error("Multivariate CDF: An invalid method was handed to the Montecarlo "
              "simulation for calculating the "
              "Multivariate CDF "
              "integral.\nThe method handed is:\t%s\nThe current available methods "
              "are `genz` or `genz2`"
              ".\n",
              method);
    }
    // ---...--- //
}

/*
 * @brief Gets the a matrix of mu values and the inverse sigma for a given
 * ballot box
 *
 * @param[in] b The ballot box index.
 * @param[in] *probabilitiesReduced Probabilities matrix without the last
 * candidate.
 * @param[in, out] **cholesky An array of `g` matrices that stores (`c-1`X`c-1`)
 * matrices with the cholesky values. The cholesky matrix is filled on the upper
 * and lower triangle with the same values since it's simmetric.
 * @param[in, out] *mu A matrix of size (`g`x`c-1`) with the averages per group.
 *
 * @return void. Results to be written on cholesky and mu
 */
static void getMainParameters(EMContext *ctx, int b, const Matrix probabilitiesReduced, Matrix **cholesky, Matrix *mu)
{
    // ---- Compute mu and sigma ---- //
    getAverageConditional(ctx, b, &probabilitiesReduced, mu, cholesky);

    // ---- Cholesky factorization (if applicable) ---- //
    if (TOTAL_CANDIDATES != 2)
    {
        for (uint16_t g = 0; g < TOTAL_GROUPS; g++)
        {
            choleskyMat(cholesky[g]);
        }
    }
}

/*
 * @brief Computes the `q` values from the multivariate CDF method
 *
 * @param[in] *probabilities A pointer towards the current probabilities matrix.
 * @param[in] monteCarloSamples The amount of samples to use in the Monte Carlo
 * simulation
 * @param[in] epsilon The error threshold used for the genz Montecarlo method.
 * @param[in] *method The method for calculating the Montecarlo simulation.
 * Currently available methods are `Genz` and `Genz1`.
 *
 * @return A contiguos array with all the new probabilities
 */
void computeQMultivariateCDF(EMContext *ctx, QMethodInput params, double *ll)
{
    *ll = 0.0;

    Matrix *X = &ctx->X;
    double *q = ctx->q;
    Matrix *P = &ctx->probabilities;

    const int B = (int)ctx->B;
    const int G = (int)ctx->G;
    const int C = (int)ctx->C;
    const int d = C - 1;

    int monteCarloSamples = params.monteCarloIter;
    double epsilon = params.errorThreshold;
    const char *method = params.simulationMethod;

    // ---- Remove last column (needed for getAverageConditional) ---- //
    Matrix P_red = removeLastColumn(P);
    // allocateSeed(ctx, params);

    // ---- Optional log-likelihood array ---- //
    double *logArray = NULL;
    if (params.computeLL)
        logArray = (double *)Calloc(B, double);

    // ---- Arena reused across all ballots ---- //
    ArenaCDF A = ArenaCDF_init(C);

    // #ifdef _OPENMP
    // #pragma omp parallel for schedule(static) if (C >= 3 && G >= 3)
    // #endif
    for (uint32_t b = 0; b < (uint32_t)B; b++)
    {
        // ---- Create per-ballot structures ---- //
        Matrix **L = (Matrix **)Calloc(G, Matrix *);
        for (int g = 0; g < G; g++)
        {
            L[g] = (Matrix *)Calloc(1, Matrix);
            *L[g] = createMatrix(d, d); // The cholesky matrix of dimension C-1 x C - 1
        }
        Matrix mu = createMatrix(G, d); // Mu matrix of dimension G x C-1

        // ---- Fill mu and cholesky (no alloc inside) ---- //
        getMainParameters(ctx, b, P_red, L, &mu);

        // ---- feature = X[, b] into Arena buffer ---- //
        getColumn_into(X, b, A.feature);

        for (uint16_t g = 0; g < G; g++)
        {
            Matrix *Lg = L[g];
            getRow_into(&mu, g, A.mu_row);

            double denom = 0.0;

            for (uint16_t c = 0; c < C; c++)
            {
                // ---- Copy and adjust hypercube bounds ---- //
                memcpy(A.featureA, A.feature, d * sizeof(double));
                memcpy(A.featureB, A.feature, d * sizeof(double));
                // Rprintf("Los valores de featureA y featureB son:\n");
                for (uint16_t k = 0; k < d; k++)
                {
                    A.featureA[k] -= 0.5;
                    A.featureB[k] += 0.5;
                    if (k == c)
                    {
                        A.featureA[k] -= 1.0;
                        A.featureB[k] -= 1.0;
                    }
                    A.featureA[k] -= A.mu_row[k];
                    A.featureB[k] -= A.mu_row[k];
                }

                // ---- Monte Carlo integral or closed form ---- //
                double val;
                if (C != 2)
                {
                    val = Montecarlo(ctx, Lg, A.featureA, A.featureB, d, monteCarloSamples, epsilon, method);
                }
                else
                {
                    double s = sqrt(MATRIX_AT_PTR(Lg, 0, 0));
                    val = pnorm(A.featureB[0], 0.0, s, 1, 0) - pnorm(A.featureA[0], 0.0, s, 1, 0);
                }

                // ---- Fallback with PDF midpoint ---- //
                if (!(val > 0.0) || !isfinite(val))
                {
                    val = pdf_midpoint(Lg, A.mu_row, A.featureA, A.featureB, d);
                    // Rprintf("Fallback!\n");
                }

                // ---- Weight by prior probability ---- //
                val *= MATRIX_AT_PTR(P, g, c);

                A.mc_results[c] = val;
                denom += isfinite(val) ? val : 0.0;

                if (params.computeLL && g == 0 && isfinite(val))
                    logArray[b] += val;
            }

            // ---- Normalize and store q ---- //
            for (uint16_t c = 0; c < C; c++)
            {
                double result = A.mc_results[c] / denom;
                Q_3D(q, b, g, c, G, C) = (isfinite(result) ? result : 0.0);
            }
        }

        // ---- Free per-ballot matrices ---- //
        for (int g = 0; g < G; g++)
        {
            freeMatrix(L[g]);
            Free(L[g]);
        }
        Free(L);
        freeMatrix(&mu);
    }

    // ---- Cleanup global ---- //
    freeMatrix(&P_red);
    ArenaCDF_free(&A);

    if (params.computeLL)
    {
        double finalLL = 0.0;
        for (uint32_t b = 0; b < (uint32_t)B; b++)
            finalLL += (logArray[b] != 0.0) ? log(logArray[b]) : 0.0;
        *ll = finalLL;
        Free(logArray);
    }
}
