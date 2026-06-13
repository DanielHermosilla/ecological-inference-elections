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
#include "parametric_main.h"
#include <R.h>
#include <R_ext/Random.h>
#include <Rcpp.h>
#include <Rinternals.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

#ifndef Calloc
#define Calloc(n, type) ((type *)R_chk_calloc((size_t)(n), sizeof(type)))
#endif

#ifndef Free
#define Free(p) R_chk_free((void *)(p))
#endif

// ---- Helper Function: Convert R matrix to C matrix ---- //
static Matrix convertToMatrix(const Rcpp::NumericMatrix &mat)
{
    int rows = mat.nrow();
    int cols = mat.ncol();
    double *data = Calloc(rows * cols, double);
    if (data == NULL)
    {
        Rcpp::stop("Error: Failed to allocate memory for matrix conversion.");
    }
    std::memcpy(data, mat.begin(), rows * cols * sizeof(double));
    return {data, rows, cols};
}

// ---- Helper Function: Release matrix data ---- //
static void releaseMatrix(Matrix *mat)
{
    if (mat == NULL)
        return;
    if (mat->data != NULL)
    {
        Free(mat->data);
        mat->data = NULL;
    }
    mat->rows = 0;
    mat->cols = 0;
}

static Matrix *convertNumericArrayToMatrixArray(const Rcpp::NumericVector &arr, int rows, int cols, int K,
                                                const char *name)
{
    Rcpp::IntegerVector dims = arr.attr("dim");
    if (dims.size() != 3 || dims[0] != rows || dims[1] != cols || dims[2] != K)
    {
        Rcpp::stop("%s must have dimensions (%d x %d x %d).", name, rows, cols, K);
    }

    Matrix *out = Calloc(K, Matrix);
    const double *src = arr.begin();
    size_t slice = (size_t)rows * (size_t)cols;
    for (int k = 0; k < K; ++k)
    {
        out[k].rows = rows;
        out[k].cols = cols;
        out[k].data = Calloc(slice, double);
        std::memcpy(out[k].data, src + slice * (size_t)k, slice * sizeof(double));
    }
    return out;
}

static void releaseMatrixArray(Matrix *arr, int K)
{
    if (arr == NULL)
        return;
    for (int k = 0; k < K; ++k)
        releaseMatrix(&arr[k]);
    Free(arr);
}

// ---- Helper Function: Initialize QMethodInput ---- //
QMethodInput initializeQMethodInput(const std::string &EMAlg, int samples, int step_size, int monte_iter,
                                    double monte_error, int miniterations, const std::string &monte_method,
                                    bool compute_ll, const std::string &LP_method, bool project_every, bool symmetric,
                                    const std::string &symmetric_weight_method)
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
    inputParams.symmetric = symmetric;
    inputParams.symmetric_weight_method = strdup(symmetric_weight_method.c_str());

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
                               monte_method, compute_ll[0], LP_method, project_every[0], false, "average");

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
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    // ---- Create human-readable stopping reason ---- //
    const char *stop_reasons[] = {"Converged", "Maximum time reached", "Maximum iterations reached"};
    const char *stopping_reason = (finish >= 0 && finish < 3) ? stop_reasons[finish] : "Unknown";

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
    releaseMatrix(&X);
    releaseMatrix(&W);
    releaseMatrix(&P);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = logLLarr,
                              Rcpp::_["total_iterations"] = totalIter, Rcpp::_["total_time"] = timeIter,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb, Rcpp::_["expected_outcome"] = expectedOut);
}

// [[Rcpp::export]]
double EMLogLikFromProb(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                        Rcpp::NumericMatrix probability_matrix, Rcpp::String em_method,
                        Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                        Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter,
                        Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                        Rcpp::LogicalVector project_every)
{
    std::string EMAlg = em_method;

    Matrix X;
    Matrix W;
    Matrix P = convertToMatrix(probability_matrix);
    RsetParameters(candidate_matrix, group_matrix, &X, &W);

    QMethodInput inputParams = initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0],
                                                      miniterations[0], monte_method, true, LP_method,
                                                      project_every[0], false, "average");

    double ll = computeLogLikForProbability(&X, &W, &P, EMAlg.c_str(), &inputParams);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    releaseMatrix(&X);
    releaseMatrix(&W);
    releaseMatrix(&P);

    return ll;
}

// [[Rcpp::export]]
double EMLogLikFromProbMixture(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                               Rcpp::NumericVector probability_array, Rcpp::String em_method,
                               Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples,
                               Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                               Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations,
                               Rcpp::String LP_method, Rcpp::LogicalVector project_every)
{
    std::string EMAlg = em_method;

    Matrix X;
    Matrix W;
    RsetParameters(candidate_matrix, group_matrix, &X, &W);

    Rcpp::IntegerVector dims = probability_array.attr("dim");
    if (dims.size() != 3)
    {
        releaseMatrix(&X);
        releaseMatrix(&W);
        Rcpp::stop("probability_array must be a 3d array with dimensions (g x c x K).");
    }
    int K = dims[2];
    Matrix *componentProb =
        convertNumericArrayToMatrixArray(probability_array, W.cols, candidate_matrix.nrow(), K, "probability_array");

    QMethodInput inputParams = initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0],
                                                      miniterations[0], monte_method, true, LP_method,
                                                      project_every[0], false, "average");

    double ll = computeLogLikForMixtureProbability(&X, &W, componentProb, K, EMAlg.c_str(), &inputParams);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    releaseMatrixArray(componentProb, K);
    releaseMatrix(&X);
    releaseMatrix(&W);

    return ll;
}

// ---- Run Finite-Mixture EM (non-parametric only) ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmMixture(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                              Rcpp::String em_method, Rcpp::String probability_method,
                              Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                              Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                              Rcpp::LogicalVector compute_ll, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                              Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                              Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                              Rcpp::LogicalVector project_every, Rcpp::RObject initial_probabilities,
                              Rcpp::IntegerVector mixture_h, Rcpp::LogicalVector symmetric,
                              Rcpp::String symmetric_weight_method)
{
    int K = (mixture_h.size() > 0) ? mixture_h[0] : 1;
    if (K < 1)
    {
        Rcpp::stop("EMAlgorithmMixture: 'mixture_h' must be a positive integer.");
    }

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;
    std::string symmetricWeightMethod = symmetric_weight_method;

    Matrix X;
    Matrix W;
    Matrix P = {NULL, 0, 0};
    Matrix *componentProbInit = NULL;
    RsetParameters(candidate_matrix, group_matrix, &X, &W);

    Rcpp::IntegerVector initialDims = initial_probabilities.attr("dim");
    if (initialDims.size() == 3)
    {
        componentProbInit =
            convertNumericArrayToMatrixArray(Rcpp::as<Rcpp::NumericVector>(initial_probabilities), W.cols,
                                             candidate_matrix.nrow(), K, "initial_prob");
    }
    else
    {
        Rcpp::NumericMatrix initialMatrix(initial_probabilities);
        P = convertToMatrix(initialMatrix);
    }

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0], symmetric[0],
                               symmetricWeightMethod);

    EMMixtureResult *res =
        EMAlgoritmMixture(&X, &W, probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0],
                          log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0], &P,
                          componentProbInit,
                          &inputParams, K);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    if (res == NULL)
    {
        releaseMatrix(&X);
        releaseMatrix(&W);
        releaseMatrix(&P);
        releaseMatrixArray(componentProbInit, K);
        Rcpp::stop("EMAlgorithmMixture: backend returned NULL.");
    }

    const char *stop_reasons[] = {"Converged", "Maximum time reached", "Maximum iterations reached"};
    const char *stopping_reason =
        (res->finish_id >= 0 && res->finish_id < 3) ? stop_reasons[res->finish_id] : "Unknown";

    const bool do_ll = (compute_ll.size() > 0) ? static_cast<bool>(compute_ll[0]) : true;
    const double out_ll = do_ll ? res->log_likelihood : NA_REAL;

    int C = candidate_matrix.nrow();
    int B = candidate_matrix.ncol();
    int G = group_matrix.ncol();
    int outK = res->mixture_h;

    Rcpp::NumericMatrix RfinalProbability(res->probabilities.rows, res->probabilities.cols);
    std::memcpy(RfinalProbability.begin(), res->probabilities.data,
                (size_t)res->probabilities.rows * (size_t)res->probabilities.cols * sizeof(double));

    std::size_t N = std::size_t(B) * G * C;
    Rcpp::NumericVector condProb(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        condProb[i] = res->q[i];
    }
    Rcpp::NumericVector expectedOut(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        expectedOut[i] = res->predicted_votes[i];
    }
    condProb.attr("dim") = Rcpp::IntegerVector::create(G, C, B);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(G, C, B);

    Rcpp::NumericVector phi(outK);
    for (int k = 0; k < outK; ++k)
    {
        phi[k] = res->phi[k];
    }

    Rcpp::NumericMatrix responsibilities(B, outK);
    std::memcpy(responsibilities.begin(), res->responsibilities.data, (size_t)B * (size_t)outK * sizeof(double));

    Rcpp::NumericVector component_prob((R_xlen_t)G * C * outK);
    for (std::size_t i = 0; i < (std::size_t)G * (std::size_t)C * (std::size_t)outK; ++i)
    {
        component_prob[i] = res->component_prob[i];
    }
    component_prob.attr("dim") = Rcpp::IntegerVector::create(G, C, outK);

    Rcpp::RObject probInv = R_NilValue;
    Rcpp::RObject condProbInv = R_NilValue;
    Rcpp::RObject expectedOutInv = R_NilValue;
    if (res->probabilities_inv.data != NULL && res->q_inv != NULL && res->predicted_votes_inv != NULL)
    {
        Rcpp::NumericMatrix reverseProbability(res->probabilities_inv.rows, res->probabilities_inv.cols);
        std::memcpy(reverseProbability.begin(), res->probabilities_inv.data,
                    (size_t)res->probabilities_inv.rows * (size_t)res->probabilities_inv.cols * sizeof(double));
        probInv = reverseProbability;

        std::size_t Ninv = std::size_t(B) * C * G;
        Rcpp::NumericVector reverseCondProb(Ninv);
        Rcpp::NumericVector reverseExpectedOut(Ninv);
        for (std::size_t i = 0; i < Ninv; ++i)
        {
            reverseCondProb[i] = res->q_inv[i];
            reverseExpectedOut[i] = res->predicted_votes_inv[i];
        }
        reverseCondProb.attr("dim") = Rcpp::IntegerVector::create(C, G, B);
        reverseExpectedOut.attr("dim") = Rcpp::IntegerVector::create(C, G, B);
        condProbInv = reverseCondProb;
        expectedOutInv = reverseExpectedOut;
    }

    int iter = res->total_iterations;
    double elapsed = res->total_time;
    int finish = res->finish_id;

    cleanupMixtureResult(res);
    releaseMatrix(&X);
    releaseMatrix(&W);
    releaseMatrix(&P);
    releaseMatrixArray(componentProbInit, K);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = out_ll,
                              Rcpp::_["total_iterations"] = iter, Rcpp::_["total_time"] = elapsed,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb, Rcpp::_["expected_outcome"] = expectedOut,
                              Rcpp::_["phi"] = phi, Rcpp::_["responsibilities"] = responsibilities,
                              Rcpp::_["component_prob"] = component_prob, Rcpp::_["prob_inv"] = probInv,
                              Rcpp::_["q_inv"] = condProbInv, Rcpp::_["expected_outcome_inv"] = expectedOutInv);
}

// ---- Run Row-Level Finite-Mixture EM (non-parametric only) ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmRowMixture(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::String em_method, Rcpp::String probability_method,
                                 Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                                 Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                                 Rcpp::LogicalVector compute_ll, Rcpp::LogicalVector verbose,
                                 Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter,
                                 Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
                                 Rcpp::LogicalVector project_every, Rcpp::NumericMatrix initial_probabilities,
                                 Rcpp::IntegerVector row_mixture_h)
{
    int H = (row_mixture_h.size() > 0) ? row_mixture_h[0] : 1;
    if (H < 1)
    {
        Rcpp::stop("EMAlgorithmRowMixture: 'row_mixture_h' must be a positive integer.");
    }

    std::string probabilityM = probability_method;
    std::string EMAlg = em_method;

    Matrix X;
    Matrix W;
    Matrix P = convertToMatrix(initial_probabilities);
    RsetParameters(candidate_matrix, group_matrix, &X, &W);

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0], false, "average");

    EMMixtureResult *res =
        EMAlgoritmRowMixture(&X, &W, probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0],
                             log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0], &P,
                             &inputParams, H);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    if (res == NULL)
    {
        releaseMatrix(&X);
        releaseMatrix(&W);
        releaseMatrix(&P);
        Rcpp::stop("EMAlgorithmRowMixture: backend returned NULL.");
    }

    const char *stop_reasons[] = {"Converged", "Maximum time reached", "Maximum iterations reached"};
    const char *stopping_reason =
        (res->finish_id >= 0 && res->finish_id < 3) ? stop_reasons[res->finish_id] : "Unknown";

    const bool do_ll = (compute_ll.size() > 0) ? static_cast<bool>(compute_ll[0]) : true;
    const double out_ll = do_ll ? res->log_likelihood : NA_REAL;

    int C = candidate_matrix.nrow();
    int B = candidate_matrix.ncol();
    int G = group_matrix.ncol();
    int outH = res->mixture_h;
    int assignment_count = res->responsibilities.cols;

    Rcpp::NumericMatrix RfinalProbability(res->probabilities.rows, res->probabilities.cols);
    std::memcpy(RfinalProbability.begin(), res->probabilities.data,
                (size_t)res->probabilities.rows * (size_t)res->probabilities.cols * sizeof(double));

    std::size_t N = std::size_t(B) * G * C;
    Rcpp::NumericVector condProb(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        condProb[i] = res->q[i];
    }
    Rcpp::NumericVector expectedOut(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        expectedOut[i] = res->predicted_votes[i];
    }
    condProb.attr("dim") = Rcpp::IntegerVector::create(G, C, B);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(G, C, B);

    Rcpp::NumericMatrix phi(G, outH);
    for (int g = 0; g < G; ++g)
    {
        for (int h = 0; h < outH; ++h)
        {
            phi(g, h) = res->phi[(size_t)g * (size_t)outH + (size_t)h];
        }
    }

    Rcpp::NumericMatrix responsibilities(B, assignment_count);
    std::memcpy(responsibilities.begin(), res->responsibilities.data,
                (size_t)B * (size_t)assignment_count * sizeof(double));

    Rcpp::NumericVector component_prob((R_xlen_t)G * C * outH);
    for (std::size_t i = 0; i < (std::size_t)G * (std::size_t)C * (std::size_t)outH; ++i)
    {
        component_prob[i] = res->component_prob[i];
    }
    component_prob.attr("dim") = Rcpp::IntegerVector::create(G, C, outH);

    int iter = res->total_iterations;
    double elapsed = res->total_time;
    int finish = res->finish_id;

    cleanupMixtureResult(res);
    releaseMatrix(&X);
    releaseMatrix(&W);
    releaseMatrix(&P);

    return Rcpp::List::create(Rcpp::_["result"] = RfinalProbability, Rcpp::_["log_likelihood"] = out_ll,
                              Rcpp::_["total_iterations"] = iter, Rcpp::_["total_time"] = elapsed,
                              Rcpp::_["stopping_reason"] = stopping_reason, Rcpp::_["finish_id"] = finish,
                              Rcpp::_["q"] = condProb, Rcpp::_["expected_outcome"] = expectedOut,
                              Rcpp::_["phi"] = phi, Rcpp::_["responsibilities"] = responsibilities,
                              Rcpp::_["component_prob"] = component_prob,
                              Rcpp::_["assignment_count"] = assignment_count);
}

// ---- Run Bootstrapping Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                        Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method,
                        Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds,
                        Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold,
                        Rcpp::LogicalVector compute_ll, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size,
                        Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                        Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations, Rcpp::String LP_method,
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
                               monte_method, compute_ll[0], LP_method, project_every[0], false, "average");

    Matrix avgProb = {NULL, 0, 0};
    Matrix sdResult =
        bootstrapA(&XR, &WR, nboot[0], EMAlg.c_str(), probabilityM.c_str(), stopping_threshold[0],
                   log_stopping_threshold[0], maximum_iterations[0], maximum_seconds[0], verbose[0], &P,
                   &inputParams, &avgProb);
    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    // Convert to R's matrix
    Rcpp::NumericMatrix output(sdResult.rows, sdResult.cols);

    std::memcpy(output.begin(), // where to copy
                sdResult.data,  // source
                sdResult.rows * sdResult.cols * sizeof(double));

    Rcpp::NumericMatrix avg_output(avgProb.rows, avgProb.cols);

    std::memcpy(avg_output.begin(), // where to copy
                avgProb.data,       // source
                avgProb.rows * avgProb.cols * sizeof(double));

    freeMatrix(&sdResult);
    freeMatrix(&avgProb);
    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&P);

    return Rcpp::List::create(Rcpp::_["sd"] = output, Rcpp::_["avg_prob"] = avg_output);
}

// ---- Run Parametric EM Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmParametric(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                 Rcpp::NumericMatrix attribute_matrix, Rcpp::NumericMatrix beta,
                                 Rcpp::NumericMatrix alpha, Rcpp::IntegerVector maximum_iterations,
                                 Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector log_stopping_threshold,
                                 Rcpp::IntegerVector maximum_newton, Rcpp::LogicalVector verbose,
                                 Rcpp::String LP_method, Rcpp::LogicalVector project_every,
                                 Rcpp::String em_method, Rcpp::String monte_method,
                                 Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    if (attribute_matrix.nrow() == 0 || attribute_matrix.ncol() == 0)
        Rcpp::stop("Error: V matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    Matrix VR = convertToMatrix(attribute_matrix);
    Matrix BetaR = convertToMatrix(beta);
    Matrix AlphaR = convertToMatrix(alpha);
    std::string adjust_prob_cond_method = LP_method;
    std::string EMAlg = em_method;
    std::string monteMethod = monte_method;
    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, 1000, 3000, monte_iter[0], monte_error[0], 0, monteMethod, true,
                               adjust_prob_cond_method, project_every[0], false, "average");

    double elapsed = 0.0;
    int total_iter = 0;
    double logLikelihood = 0.0;
    Matrix *condProbMat = NULL;
    Matrix *expectedMat = NULL;
    Matrix *finalProb =
        EM_Algorithm_Method(&XR, &WR, &VR, &BetaR, &AlphaR, maximum_iterations[0], maximum_seconds[0],
                            log_stopping_threshold[0], maximum_newton[0], verbose[0], &elapsed, &total_iter,
                            &logLikelihood, &condProbMat, &expectedMat, EMAlg.c_str(), &inputParams,
                            adjust_prob_cond_method.c_str(), project_every[0]);

    // ---- Build the probability array (g x c x b) ---- //
    int S = VR.rows;
    int R = finalProb[0].rows;
    int C = finalProb[0].cols;

    Rcpp::NumericVector probArr(R * C * S);
    probArr.attr("dim") = Rcpp::IntegerVector::create(R, C, S);

    Rcpp::NumericVector condProb(R * C * S);
    condProb.attr("dim") = Rcpp::IntegerVector::create(R, C, S);

    Rcpp::NumericVector expectedOut(R * C * S);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(R, C, S);

    for (int s = 0; s < S; ++s)
    {
        for (int i = 0; i < R; ++i)
        {
            for (int j = 0; j < C; ++j)
            {
                int idx = i + R * (j + C * s);
                probArr[idx] = MATRIX_AT(finalProb[s], i, j);
                if (condProbMat != NULL)
                    condProb[idx] = MATRIX_AT(condProbMat[s], i, j);
                if (expectedMat != NULL)
                    expectedOut[idx] = MATRIX_AT(expectedMat[s], i, j);
            }
        }
        freeMatrix(&finalProb[s]);
        if (condProbMat != NULL)
            freeMatrix(&condProbMat[s]);
        if (expectedMat != NULL)
            freeMatrix(&expectedMat[s]);
    }
    Free(finalProb);
    if (condProbMat != NULL)
        Free(condProbMat);
    if (expectedMat != NULL)
        Free(expectedMat);

    // ---- Copy out final beta and alpha ---- //
    Rcpp::NumericMatrix Rbeta(BetaR.rows, BetaR.cols);
    std::memcpy(Rbeta.begin(), BetaR.data, sizeof(double) * BetaR.rows * BetaR.cols);

    Rcpp::NumericMatrix Ralpha(AlphaR.rows, AlphaR.cols);
    std::memcpy(Ralpha.begin(), AlphaR.data, sizeof(double) * AlphaR.rows * AlphaR.cols);

    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&VR);
    releaseMatrix(&BetaR);
    releaseMatrix(&AlphaR);
    if (inputParams.simulationMethod != nullptr)
        free((void *)inputParams.simulationMethod);
    if (inputParams.prob_cond != nullptr)
        free((void *)inputParams.prob_cond);
    if (inputParams.symmetric_weight_method != nullptr)
        free((void *)inputParams.symmetric_weight_method);

    return Rcpp::List::create(Rcpp::_["prob"] = probArr, Rcpp::_["cond_prob"] = condProb,
                              Rcpp::_["expected_outcome"] = expectedOut, Rcpp::_["beta"] = Rbeta,
                              Rcpp::_["alpha"] = Ralpha, Rcpp::_["time"] = elapsed,
                              Rcpp::_["iter"] = total_iter, Rcpp::_["logLik"] = logLikelihood);
}

// ---- Run Parametric Finite-Mixture EM Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List EMAlgorithmParametricMixture(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                        Rcpp::NumericMatrix attribute_matrix, Rcpp::NumericVector beta,
                                        Rcpp::NumericVector alpha, Rcpp::IntegerVector maximum_iterations,
                                        Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold,
                                        Rcpp::NumericVector log_stopping_threshold,
                                        Rcpp::IntegerVector maximum_newton, Rcpp::LogicalVector compute_ll,
                                        Rcpp::LogicalVector verbose, Rcpp::String em_method,
                                        Rcpp::String monte_method, Rcpp::NumericVector monte_error,
                                        Rcpp::IntegerVector monte_iter, Rcpp::IntegerVector miniterations,
                                        Rcpp::String LP_method, Rcpp::LogicalVector project_every,
                                        Rcpp::IntegerVector mixture_h)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");
    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");
    if (attribute_matrix.nrow() == 0 || attribute_matrix.ncol() == 0)
        Rcpp::stop("Error: V matrix has zero dimensions!");

    int K = (mixture_h.size() > 0) ? mixture_h[0] : 1;
    if (K < 1)
        Rcpp::stop("EMAlgorithmParametricMixture: 'mixture_h' must be a positive integer.");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    Matrix VR = convertToMatrix(attribute_matrix);
    int C = XR.cols;
    int B = XR.rows;
    int G = WR.cols;
    int A = VR.cols;

    Rcpp::IntegerVector betaDims = beta.attr("dim");
    if (betaDims.size() != 2 || betaDims[0] != A || betaDims[1] != K - 1)
        Rcpp::stop("beta must have dimensions (%d x %d) in parametric matrix-mixture mode.", A, K - 1);
    Matrix BetaR;
    BetaR.rows = A;
    BetaR.cols = K - 1;
    BetaR.data = Calloc((size_t)A * (size_t)(K - 1), double);
    std::memcpy(BetaR.data, beta.begin(), (size_t)A * (size_t)(K - 1) * sizeof(double));

    Matrix *ComponentProbR = convertNumericArrayToMatrixArray(alpha, G, C, K, "component_prob");

    std::string EMAlg = em_method;
    std::string monteMethod = monte_method;
    std::string adjust_prob_cond_method = LP_method;
    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, 1000, 3000, monte_iter[0], monte_error[0], miniterations[0], monteMethod,
                               compute_ll[0], adjust_prob_cond_method, project_every[0], false, "average");

    EMParametricMixtureResult *res = EM_Algorithm_Parametric_Mixture(
        &XR, &WR, &VR, ComponentProbR, &BetaR, K, maximum_iterations[0], miniterations[0], maximum_seconds[0],
        stopping_threshold[0], log_stopping_threshold[0], maximum_newton[0], verbose[0], EMAlg.c_str(), &inputParams,
        adjust_prob_cond_method.c_str(), project_every[0]);

    Rcpp::NumericVector probArr((R_xlen_t)G * C * B);
    Rcpp::NumericVector condProb((R_xlen_t)G * C * B);
    Rcpp::NumericVector expectedOut((R_xlen_t)G * C * B);
    probArr.attr("dim") = Rcpp::IntegerVector::create(G, C, B);
    condProb.attr("dim") = Rcpp::IntegerVector::create(G, C, B);
    expectedOut.attr("dim") = Rcpp::IntegerVector::create(G, C, B);

    for (int b = 0; b < B; ++b)
    {
        for (int c = 0; c < C; ++c)
        {
            for (int g = 0; g < G; ++g)
            {
                int idx = g + G * (c + C * b);
                probArr[idx] = MATRIX_AT(res->probabilities[b], g, c);
                condProb[idx] = MATRIX_AT(res->q[b], g, c);
                expectedOut[idx] = MATRIX_AT(res->expected[b], g, c);
            }
        }
    }

    Rcpp::NumericMatrix phi(B, K);
    std::memcpy(phi.begin(), res->phi.data, (size_t)B * (size_t)K * sizeof(double));

    Rcpp::NumericMatrix responsibilities(B, K);
    std::memcpy(responsibilities.begin(), res->responsibilities.data, (size_t)B * (size_t)K * sizeof(double));

    Rcpp::NumericVector componentProb((R_xlen_t)G * C * K);
    std::memcpy(componentProb.begin(), res->component_prob, (size_t)G * (size_t)C * (size_t)K * sizeof(double));
    componentProb.attr("dim") = Rcpp::IntegerVector::create(G, C, K);

    Rcpp::NumericMatrix betaArr(A, K - 1);
    std::memcpy(betaArr.begin(), res->beta.data, (size_t)A * (size_t)(K - 1) * sizeof(double));

    int iter = res->total_iterations;
    double elapsed = res->total_time;
    double logLikelihood = compute_ll[0] ? res->log_likelihood : NA_REAL;
    int finish = res->finish_id;
    const char *stop_reasons[] = {"Converged", "Maximum time reached", "Maximum iterations reached"};
    const char *stopping_reason = (finish >= 0 && finish < 3) ? stop_reasons[finish] : "Unknown";

    cleanupParametricMixtureResult(res);
    releaseMatrix(&BetaR);
    releaseMatrixArray(ComponentProbR, K);
    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&VR);
    if (inputParams.simulationMethod != nullptr)
        free((void *)inputParams.simulationMethod);
    if (inputParams.prob_cond != nullptr)
        free((void *)inputParams.prob_cond);
    if (inputParams.symmetric_weight_method != nullptr)
        free((void *)inputParams.symmetric_weight_method);

    return Rcpp::List::create(Rcpp::_["prob"] = probArr, Rcpp::_["cond_prob"] = condProb,
                              Rcpp::_["expected_outcome"] = expectedOut, Rcpp::_["beta"] = betaArr,
                              Rcpp::_["alpha"] = R_NilValue, Rcpp::_["phi"] = phi,
                              Rcpp::_["responsibilities"] = responsibilities,
                              Rcpp::_["component_prob"] = componentProb,
                              Rcpp::_["time"] = elapsed, Rcpp::_["iter"] = iter,
                              Rcpp::_["logLik"] = logLikelihood,
                              Rcpp::_["stopping_reason"] = stopping_reason,
                              Rcpp::_["finish_id"] = finish);
}

// ---- Run Parametric Bootstrapping Algorithm ---- //
// [[Rcpp::export]]
Rcpp::List bootstrapParametricAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix,
                                  Rcpp::NumericMatrix attribute_matrix, Rcpp::NumericMatrix beta,
                                  Rcpp::NumericMatrix alpha, Rcpp::IntegerVector maximum_iterations,
                                  Rcpp::IntegerVector nboot, Rcpp::NumericVector maximum_seconds,
                                  Rcpp::NumericVector log_stopping_threshold, Rcpp::IntegerVector maximum_newton,
                                  Rcpp::LogicalVector verbose, Rcpp::String LP_method,
                                  Rcpp::LogicalVector project_every)
{
    if (candidate_matrix.nrow() == 0 || candidate_matrix.ncol() == 0)
        Rcpp::stop("Error: X matrix has zero dimensions!");

    if (group_matrix.nrow() == 0 || group_matrix.ncol() == 0)
        Rcpp::stop("Error: W matrix has zero dimensions!");

    if (attribute_matrix.nrow() == 0 || attribute_matrix.ncol() == 0)
        Rcpp::stop("Error: V matrix has zero dimensions!");

    Matrix XR = convertToMatrix(candidate_matrix);
    Matrix WR = convertToMatrix(group_matrix);
    Matrix VR = convertToMatrix(attribute_matrix);
    Matrix BetaR = convertToMatrix(beta);
    Matrix AlphaR = convertToMatrix(alpha);
    std::string adjust_prob_cond_method = LP_method;

    Matrix sdBeta = {NULL, 0, 0};
    Matrix sdAlpha = {NULL, 0, 0};

    bootstrapParametric(&XR, &WR, &VR, &BetaR, &AlphaR, nboot[0], maximum_iterations[0], maximum_seconds[0],
                        log_stopping_threshold[0], maximum_newton[0], verbose[0], &sdBeta, &sdAlpha,
                        adjust_prob_cond_method.c_str(), project_every[0]);

    // ---- Convert outputs to R matrices ---- //
    Rcpp::NumericMatrix RsdBeta(sdBeta.rows, sdBeta.cols);
    std::memcpy(RsdBeta.begin(), sdBeta.data, sizeof(double) * sdBeta.rows * sdBeta.cols);

    Rcpp::NumericMatrix RsdAlpha(sdAlpha.rows, sdAlpha.cols);
    std::memcpy(RsdAlpha.begin(), sdAlpha.data, sizeof(double) * sdAlpha.rows * sdAlpha.cols);

    freeMatrix(&sdBeta);
    freeMatrix(&sdAlpha);
    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&VR);
    releaseMatrix(&BetaR);
    releaseMatrix(&AlphaR);

    return Rcpp::List::create(Rcpp::_["sd_beta"] = RsdBeta, Rcpp::_["sd_alpha"] = RsdAlpha);
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
                               monte_method, compute_ll[0], LP_method, project_every[0], false, "average");

    // We'll hold the boundary indices here
    int G = WR.cols;
    int *cuttingBuffer = Calloc(G, int);
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
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
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
    Free(cuttingBuffer);
    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&P);

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

    // Stores returned partition boundaries. The exhaustive search can return
    // one boundary per original group, so this needs G slots, not G - 1 cuts.
    int G = WR.cols;
    int *boundaries = Calloc(G, int);
    int numCuts = 0;
    Matrix *bestBootstrap = NULL;

    QMethodInput inputParams =
        initializeQMethodInput(EMAlg, samples[0], step_size[0], monte_iter[0], monte_error[0], miniterations[0],
                               monte_method, compute_ll[0], LP_method, project_every[0], false, "average");

    Matrix greedyP =
        aggregateGroupsExhaustive(&XR, &WR, boundaries, &numCuts, set_method.c_str(), nboot[0], sd_threshold[0],
                                  probabilityM.c_str(), EMAlg.c_str(), stopping_threshold[0], log_stopping_threshold[0],
                                  verbose[0], maximum_iterations[0], maximum_seconds[0], &P, &inputParams, &bestLogLL,
                                  &bestQ, &bestExpected, &bestBootstrap, &bestTime, &finishReason, &totalIter);

    if (inputParams.simulationMethod != nullptr)
    {
        free((void *)inputParams.simulationMethod);
    }
    if (inputParams.prob_cond != nullptr)
    {
        free((void *)inputParams.prob_cond);
    }
    if (inputParams.symmetric_weight_method != nullptr)
    {
        free((void *)inputParams.symmetric_weight_method);
    }

    if (numCuts == 0) // Case where there's not any match
    {
        freeMatrix(&greedyP);
        releaseMatrix(&XR);
        releaseMatrix(&WR);
        releaseMatrix(&P);
        Free(boundaries);
        return Rcpp::List::create(Rcpp::_["indices"] = Rcpp::IntegerVector::create(-1));
    }

    // ---- Create human-readable stopping reason ---- //
    const char *stop_reasons[] = {"Converged", "Log-likelihood decrease", "Maximum time reached",
                                  "Maximum iterations reached"};
    const char *stopping_reason = (finishReason >= 0 && finishReason < 4) ? stop_reasons[finishReason] : "Unknown";

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
    Free(bestQ);
    Free(bestExpected);
    freeMatrix(&greedyP);
    freeMatrix(bestBootstrap);
    Free(bestBootstrap);
    releaseMatrix(&XR);
    releaseMatrix(&WR);
    releaseMatrix(&P);

    // Convert to R's integer vector
    Rcpp::IntegerVector result(numCuts);
    for (int i = 0; i < numCuts; i++)
    {
        result[i] = boundaries[i];
    }

    Free(boundaries);

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
