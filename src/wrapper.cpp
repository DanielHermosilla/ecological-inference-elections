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

#include "wrapper.h"
#include "bootstrap.h"
#include "dynamic_program.h"
#include "exact.h"
#include "main.h"
#include <R.h>
#include <R_ext/Random.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <vector>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

Matrix convertToMatrix(const Rcpp::NumericMatrix &mat)
{
    int rows = mat.nrow(), cols = mat.ncol();
    double *data = (double *)malloc(rows * cols * sizeof(double)); // Allocate on heap
    std::memcpy(data, mat.begin(), rows * cols * sizeof(double));  // Copy data from R matrix
    return {data, rows, cols};                                     // Safe to return
}

// ---- Helper Function: Initialize QMethodInput ---- //
QMethodInput initializeQMethodInput(const std::string &EMAlg, int samples, int step_size, int monte_iter,
                                    double monte_error, int miniterations, const std::string &monte_method,
                                    bool compute_ll, const std::string &LP_method, bool project_every)
{
    QMethodInput inputParams = {0}; // Default initialization

    if (EMAlg == "mcmc")
    {
        inputParams.S = samples;
        inputParams.M = step_size;
    }
    else if (EMAlg == "mvn_cdf")
    {
        inputParams.monteCarloIter = monte_iter;
        inputParams.errorThreshold = monte_error;
        inputParams.simulationMethod = strdup(monte_method.c_str());
    }
    inputParams.miniter = miniterations;
    inputParams.computeLL = compute_ll;
    inputParams.prob_cond = strdup(LP_method.c_str());
    inputParams.prob_cond_every = project_every; // Weights by default

    return inputParams;
}

// ---- Set Parameters ---- //
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Matrix *X, Matrix *W)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    *X = convertToMatrix(candidate_matrix);
    *W = convertToMatrix(group_matrix);
}

// ---- Run EM Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmFull(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                           Rcpp::String em_method, Rcpp::String probability_method,
                           Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                           Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                           Rcpp::LogicalVector compute_ll, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                           Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                           Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                           Rcpp::LogicalVector project_every, Rcpp::NumericMatrix initial_probabilities)
{
    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;

    double timeIter = 0, logLLarr = 0;
    int totalIter = 0, finish = 0;
    Matrix X;
    Matrix W;
    Matrix P = convertToMatrix(initial_probabilities);
    RsetParameters(candidate_matrix, group_matrix, &X, &W);

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0]);

    EMContext *ctx = EMAlgoritm(&X, &W, probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0],
                                log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0],
                                &timeIter, &totalIter, &logLLarr, &finish, &P, &inputParams);

    Matrix *Pnew = &ctx->probabilities;
    double *qvalue = ctx->q;
    double *expected = ctx->predicted_votes;

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }

    // ---- Create human-readable stopping reason ---- //
    std::vector<std::string> stop_reasons = {"Converged", "Maximum time reached", "Maximum iterations reached"};
    std::string stopping_reason = (finish >= 0 && finish < 3) ? stop_reasons[finish] : "Unknown";

    Rcpp::NumericMatrix RfinalProbability(Pnew->rows, Pnew->cols, Pnew->data);

    std::size_t N = std::size_t(TOTAL_BALLOTS) * TOTAL_GROUPS * TOTAL_CANDIDATES;
    Rcpp::NumericVector condProb(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        condProb[i] = qvalue[i];
    }
    Rcpp::NumericVector expectedOut(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        expectedOut[i] = expected[i];
    }

    condProb.attr("dim") = Rcpp::IntegerVector::create(TOTAL_GROUPS, TOTAL_CANDIDATES, TOTAL_BALLOTS);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(TOTAL_GROUPS, TOTAL_CANDIDATES, TOTAL_BALLOTS);

    cleanup(ctx);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = logLLarr,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb, Rcpp::_["expected_outcome"] = expectedOut);
}

// ---- Run Bootstrapping Algorithm ---- //
// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                                 Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                                 Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                                 Rcpp::LogicalVector compute_ll, Rcpp::LogicalVector verbose,
                                 Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter,
                                 Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                                 Rcpp::LogicalVector project_every, Rcpp::NumericMatrix initial_probabilities)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    Matrix P = convertToMatrix(initial_probabilities);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    // cleanGlobals(EMAlg, false); // Cleans leftovers

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0]);

    Matrix sdResult =
        bootstrapA(&XR, &WR, nboot[0], EMAlg.c_str(), probabilityM.c_str(), stopping_threshold[0],
                   log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0], &P, &inputParams);
    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }

    // Convert to R's matrix
    Rcpp::NumericMatrix output(sdResult.rows, sdResult.cols);

    std::memcpy(output.begin(), // where to copy
                sdResult.data,  // source
                sdResult.rows * sdResult.cols * sizeof(double));

    freeMatrix(&sdResult);

    return output;
}

// ---- Run Group Aggregation Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List groupAgg(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold, Rcpp::LogicalVector feasible,
                    Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Rcpp::IntegerVector nboot,
                    Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations,
                    Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold,
                    Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector compute_ll,
                    Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples,
                    Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter,
                    Rcpp::IntegerVector miniterations, Rcpp::String LP_method, Rcpp::LogicalVector project_every,
                    Rcpp::NumericMatrix initial_probabilities)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    // Matrix P = createMatrix(WR.cols, XR.rows); // Empty initial probabilities
    Matrix P = convertToMatrix(initial_probabilities);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    std::string aggMet = sd_statistic;
    // cleanGlobals(EMAlg, false); // Cleans leftovers

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0]);

    // We'll hold the boundary indices here
    int G = WR.cols;
    int *cuttingBuffer = new int[G];
    int usedCuts = 0; // how many boundaries we actually use
    bool bestResult = false;

    EMContext *ctx = createEMContext(&XR, &WR, EMAlg.c_str(), inputParams);
    Matrix sdResult =
        aggregateGroups(ctx, cuttingBuffer, &usedCuts, &bestResult, sd_threshold[0], aggMet.c_str(), feasible[0],
                        nboot[0], probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0],
                        maximum_iterations[0], maximum_seconds[0], verbose[0], &P, &inputParams);
    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    // Convert to R's matrix
    Rcpp::NumericMatrix output(sdResult.rows, sdResult.cols);

    std::memcpy(output.begin(), // where to copy
                sdResult.data,  // source
                sdResult.rows * sdResult.cols * sizeof(double));
    // Convert to R's integer vector
    Rcpp::IntegerVector result(usedCuts);
    for (int i = 0; i < usedCuts; i++)
    {
        result[i] = cuttingBuffer[i];
    }

    // Free native memory
    freeMatrix(&sdResult);
    // free(cuttingBuffer);
    delete[] cuttingBuffer;

    return Rcpp::List::create(Rcpp::_["bootstrap_result"] = output, Rcpp::_["indices"] = result,
                              Rcpp::_["best_result"] = bestResult);
}

// ---- Run Greedy Group Aggregation Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List groupAggGreedy(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold,
                          Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                          Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                          Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                          Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                          Rcpp::NumericVector compute_ll, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                          Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                          Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                          Rcpp::LogicalVector project_every)
{

    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    Matrix P = createMatrix(WR.cols, XR.rows); // Empty initial probabilities

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    std::string set_method = sd_statistic;
    // cleanGlobals(EMAlg, false); // Cleans leftovers

    // Prepare the out-parameters as C++ local variables
    double bestLogLL = 0.0;
    double bestTime = 0.0;
    double *bestQ = NULL;
    double *bestExpected = NULL;
    int finishReason = 0;
    int totalIter = 0;

    // For storing the final partition boundaries
    // 'G' can be up to WR.cols. You might want a vector of size G.
    int G = WR.cols;
    int *boundaries = new int[G - 1];
    int numCuts = 0;
    Matrix *bestBootstrap = NULL;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0]);

    Matrix greedyP =
        aggregateGroupsExhaustive(&XR, &WR, boundaries, &numCuts, set_method.c_str(), nboot[0], sd_threshold[0],
                                  probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0],
                                  verbose[0], maximum_iterations[0], maximum_seconds[0], &P, &inputParams, &bestLogLL,
                                  &bestQ, &bestExpected, &bestBootstrap, &bestTime, &finishReason, &totalIter);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }

    if (numCuts == 0) // Case where there's not any match
    {
        freeMatrix(&greedyP);
        freeMatrix(&XR);
        freeMatrix(&WR);
        return Rcpp::List::create(Rcpp::_["indices"] = Rcpp::IntegerVector::create(-1));
    }

    // ---- Create human-readable stopping reason ---- //
    std::vector<std::string> stop_reasons = {"Converged", "Log-likelihood decrease", "Maximum time reached",
                                             "Maximum iterations reached"};
    std::string stopping_reason = (finishReason >= 0 && finishReason < 4) ? stop_reasons[finishReason] : "Unknown";

    Rcpp::NumericMatrix probabilities(greedyP.rows, greedyP.cols);

    std::memcpy(probabilities.begin(), // where to copy
                greedyP.data,          // source
                greedyP.rows * greedyP.cols * sizeof(double));

    Rcpp::NumericMatrix bootstrapSol(bestBootstrap->rows, bestBootstrap->cols);

    std::memcpy(bootstrapSol.begin(), // where to copy
                bestBootstrap->data,  // source
                bestBootstrap->rows * bestBootstrap->cols * sizeof(double));

    std::size_t N = std::size_t(WR.rows) * greedyP.rows * greedyP.cols;
    Rcpp::NumericVector condProb(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        condProb[i] = bestQ[i];
    }
    Rcpp::NumericVector expectedOut(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        expectedOut[i] = bestExpected[i];
    }
    condProb.attr("dim") = Rcpp::IntegerVector::create(greedyP.rows, greedyP.cols, WR.rows);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(greedyP.rows, greedyP.cols, WR.rows);

    // condProb.attr("dim") = Rcpp::IntegerVector::create(WR.rows, greedyP.rows, XR.rows); // (b, A, c)
    free(bestQ);
    free(bestExpected);
    freeMatrix(&greedyP);
    freeMatrix(bestBootstrap);
    Free(bestBootstrap);
    freeMatrix(&XR);
    freeMatrix(&WR);

    // Convert to R's integer vector
    Rcpp::IntegerVector result(numCuts);
    for (int i = 0; i < numCuts; i++)
    {
        result[i] = boundaries[i];
    }

    delete[] boundaries;

    return Rcpp::List::create(Rcpp::_["probabilities"] = probabilities, Rcpp::_["log_likelihood"] = bestLogLL,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = bestTime,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finishReason,
                              Rcpp::_["q"] = condProb, Rcpp::_["expected_outcome"] = expectedOut,
                              Rcpp::_["indices"] = result, Rcpp::_["bootstrap_sol"] = bootstrapSol);
}

// ---- Computes the exact log-likelihood ---- //
// Rcpp::NumericVector computeExactLL(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
//                                    Rcpp::NumericMatrix prob_matrix)
// {
//     if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
//         Rcpp::stop("Error: X matrix has zero dimensions!");
//
//     if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
//         Rcpp::stop("Error: W matrix has zero dimensions!");
//
//     Matrix XR = convertToMatrix(candidate_matrix);
//     Matrix WR = convertToMatrix(group_matrix);
//     Matrix PR = convertToMatrix(prob_matrix);
//
//     EMContext *ctx = createEMContext(&XR, &WR, "exact", (QMethodInput){0});
//     ctx->probabilities = PR;
//
//     double ll = computeExactLoglikelihood(ctx);
//
//     cleanup(ctx);
//
//     Rcpp::NumericVector result(1);
//     result[0] = ll;
//
//     return result;
// }
