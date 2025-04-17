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

#include "wrapper.h"
#include "bootstrap.h"
#include "dynamic_program.h"
#include "main.h"
#include <Rcpp.h>
#include <vector>

Matrix convertToMatrix(const Rcpp::NumericMatrix &mat)
{
    int rows = mat.nrow(), cols = mat.ncol();
    double *data = (double *)malloc(rows * cols * sizeof(double)); // Allocate on heap
    std::memcpy(data, mat.begin(), rows * cols * sizeof(double));  // Copy data from R matrix
    return {data, rows, cols};                                     // Safe to return
}

// ---- Helper Function: Initialize QMethodInput ---- //
QMethodInput initializeQMethodInput(const std::string &EMAlg, int samples, int step_size, int monte_iter,
                                    double monte_error, const std::string &monte_method)
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
        inputParams.simulationMethod = monte_method.c_str();
    }

    return inputParams;
}

void cleanGlobals(const std::string &EMAlg, bool everything)
{
    // Everything alludes to clean RsetParameters, if it's false, cleans leftovers
    if (everything)
        cleanup();
    if (EMAlg == "mcmc")
        cleanHitAndRun();
    else if (EMAlg == "exact")
        cleanExact();
}

// ---- Set Parameters ---- //
// [[Rcpp::export]]
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    setParameters(&XR, &WR);
}

// ---- Run EM Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmFull(Rcpp::String em_method, Rcpp::String probability_method,
                           Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                           Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                           Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples,
                           Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    cleanGlobals(EMAlg, false); // Cleans leftovers from previous iterations, in case there was an abort.

    Matrix pIn = getInitialP(probabilityM.c_str());

    double timeIter = 0;
    int totalIter = 0, finish = 0;
    double *qvalue = NULL;
    double logLLarr = 0;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    Matrix Pnew =
        EMAlgoritm(&pIn, EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0], maximum_iterations[0],
                   maximum_seconds[0], verbose[0], &timeIter, &totalIter, &logLLarr, &qvalue, &finish, inputParams);

    // ---- Create human-readable stopping reason ---- //
    std::vector<std::string> stop_reasons = {"Converged", "Log-likelihood decrease", "Maximum time reached",
                                             "Maximum iterations reached"};
    std::string stopping_reason = (finish >= 0 && finish < 4) ? stop_reasons[finish] : "Unknown";

    Rcpp::NumericMatrix RfinalProbability(Pnew.rows, Pnew.cols, Pnew.data);
    freeMatrix(&Pnew);

    Rcpp::NumericVector condProb(qvalue, qvalue + TOTAL_BALLOTS * TOTAL_CANDIDATES * TOTAL_GROUPS);
    condProb.attr("dim") = Rcpp::IntegerVector::create(TOTAL_BALLOTS, TOTAL_GROUPS, TOTAL_CANDIDATES);
    // Rcpp::NumericVector logArray(logLLarr, logLLarr + totalIter);
    free(qvalue);
    cleanGlobals(EMAlg, true);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = logLLarr,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb);
}

// ---- Run Bootstrapping Algorithm ---- //
// [[Rcpp::export]]
Rcpp::NumericMatrix bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                                 Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                                 Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                                 Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                                 Rcpp::IntegerVector samples, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    cleanGlobals(EMAlg, false); // Cleans leftovers

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    Matrix sdResult =
        bootstrapA(&XR, &WR, nboot[0], EMAlg.c_str(), probabilityM.c_str(), stopping_threshold[0],
                   log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0], inputParams);

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
                    Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector verbose,
                    Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                    Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    std::string aggMet = sd_statistic;
    cleanGlobals(EMAlg, false); // Cleans leftovers

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    // We'll hold the boundary indices here
    int G = WR.cols;
    int *cuttingBuffer = new int[G - 1];
    int usedCuts = 0; // how many boundaries we actually use
    bool bestResult = false;
    Matrix sdResult =
        aggregateGroups(&XR, &WR, cuttingBuffer, &usedCuts, &bestResult, sd_threshold[0], aggMet.c_str(), feasible[0],
                        nboot[0], probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0],
                        maximum_iterations[0], maximum_seconds[0], verbose[0], inputParams);

    // Convert to R's matrix
    Rcpp::NumericMatrix output(sdResult.rows, sdResult.cols);

    std::memcpy(output.begin(), // where to copy
                sdResult.data,  // source
                sdResult.rows * sdResult.cols * sizeof(double));

    /*
    if (usedCuts == -2)
    {
        delete[] cuttingBuffer;
        return Rcpp::List::create(Rcpp::_["bootstrap_result"] = output, Rcpp::_["indices"] = usedCuts,
                                  Rcpp::_["best_result"] = bestResult);
    }
    */

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
                          Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples,
                          Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{

    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    std::string set_method = sd_statistic;
    cleanGlobals(EMAlg, false); // Cleans leftovers

    // Prepare the out-parameters as C++ local variables
    double bestLogLL = 0.0;
    double bestTime = 0.0;
    double *bestQ = NULL;
    int finishReason = 0;
    int totalIter = 0;

    // For storing the final partition boundaries
    // 'G' can be up to WR.cols. You might want a vector of size G.
    int G = WR.cols;
    int *boundaries = new int[G - 1];
    int numCuts = 0;
    Matrix *bestBootstrap = NULL;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], monte_method);

    Matrix greedyP = aggregateGroupsExhaustive(
        &XR, &WR, boundaries, &numCuts, set_method.c_str(), nboot[0], sd_threshold[0], probabilityM.c_str(),
        EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0], verbose[0], maximum_iterations[0],
        maximum_seconds[0], inputParams, &bestLogLL, &bestQ, &bestBootstrap, &bestTime, &finishReason, &totalIter);

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

    Rcpp::NumericVector condProb(WR.rows * greedyP.rows * XR.rows); // b x g x c
    std::memcpy(condProb.begin(), bestQ, WR.rows * greedyP.rows * XR.rows * sizeof(double));

    condProb.attr("dim") = Rcpp::IntegerVector::create(WR.rows, greedyP.rows, XR.rows); // (b, A, c)
    free(bestQ);
    freeMatrix(&greedyP);
    freeMatrix(bestBootstrap);
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
                              Rcpp::_["q"] = condProb, Rcpp::_["indices"] = result,
                              Rcpp::_["bootstrap_sol"] = bootstrapSol);
}
