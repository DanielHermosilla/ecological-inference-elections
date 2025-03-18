// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RsetParameters
void RsetParameters(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix);
RcppExport SEXP _fastei_RsetParameters(SEXP candidate_matrixSEXP, SEXP group_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type candidate_matrix(candidate_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type group_matrix(group_matrixSEXP);
    RsetParameters(candidate_matrix, group_matrix);
    return R_NilValue;
END_RCPP
}
// EMAlgorithmFull
Rcpp::List EMAlgorithmFull(Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);
RcppExport SEXP _fastei_EMAlgorithmFull(SEXP em_methodSEXP, SEXP probability_methodSEXP, SEXP maximum_iterationsSEXP, SEXP maximum_secondsSEXP, SEXP stopping_thresholdSEXP, SEXP log_stopping_thresholdSEXP, SEXP verboseSEXP, SEXP step_sizeSEXP, SEXP samplesSEXP, SEXP monte_methodSEXP, SEXP monte_errorSEXP, SEXP monte_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type em_method(em_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type probability_method(probability_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maximum_iterations(maximum_iterationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type maximum_seconds(maximum_secondsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stopping_threshold(stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_stopping_threshold(log_stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type monte_method(monte_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type monte_error(monte_errorSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type monte_iter(monte_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(EMAlgorithmFull(em_method, probability_method, maximum_iterations, maximum_seconds, stopping_threshold, log_stopping_threshold, verbose, step_size, samples, monte_method, monte_error, monte_iter));
    return rcpp_result_gen;
END_RCPP
}
// bootstrapAlg
Rcpp::NumericMatrix bootstrapAlg(Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);
RcppExport SEXP _fastei_bootstrapAlg(SEXP candidate_matrixSEXP, SEXP group_matrixSEXP, SEXP nbootSEXP, SEXP em_methodSEXP, SEXP probability_methodSEXP, SEXP maximum_iterationsSEXP, SEXP maximum_secondsSEXP, SEXP stopping_thresholdSEXP, SEXP log_stopping_thresholdSEXP, SEXP verboseSEXP, SEXP step_sizeSEXP, SEXP samplesSEXP, SEXP monte_methodSEXP, SEXP monte_errorSEXP, SEXP monte_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type candidate_matrix(candidate_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type group_matrix(group_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nboot(nbootSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type em_method(em_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type probability_method(probability_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maximum_iterations(maximum_iterationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type maximum_seconds(maximum_secondsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stopping_threshold(stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_stopping_threshold(log_stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type monte_method(monte_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type monte_error(monte_errorSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type monte_iter(monte_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(bootstrapAlg(candidate_matrix, group_matrix, nboot, em_method, probability_method, maximum_iterations, maximum_seconds, stopping_threshold, log_stopping_threshold, verbose, step_size, samples, monte_method, monte_error, monte_iter));
    return rcpp_result_gen;
END_RCPP
}
// groupAgg
Rcpp::List groupAgg(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold, Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);
RcppExport SEXP _fastei_groupAgg(SEXP sd_statisticSEXP, SEXP sd_thresholdSEXP, SEXP candidate_matrixSEXP, SEXP group_matrixSEXP, SEXP nbootSEXP, SEXP em_methodSEXP, SEXP probability_methodSEXP, SEXP maximum_iterationsSEXP, SEXP maximum_secondsSEXP, SEXP stopping_thresholdSEXP, SEXP log_stopping_thresholdSEXP, SEXP verboseSEXP, SEXP step_sizeSEXP, SEXP samplesSEXP, SEXP monte_methodSEXP, SEXP monte_errorSEXP, SEXP monte_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type sd_statistic(sd_statisticSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sd_threshold(sd_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type candidate_matrix(candidate_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type group_matrix(group_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nboot(nbootSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type em_method(em_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type probability_method(probability_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maximum_iterations(maximum_iterationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type maximum_seconds(maximum_secondsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stopping_threshold(stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_stopping_threshold(log_stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type monte_method(monte_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type monte_error(monte_errorSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type monte_iter(monte_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(groupAgg(sd_statistic, sd_threshold, candidate_matrix, group_matrix, nboot, em_method, probability_method, maximum_iterations, maximum_seconds, stopping_threshold, log_stopping_threshold, verbose, step_size, samples, monte_method, monte_error, monte_iter));
    return rcpp_result_gen;
END_RCPP
}
// groupAggGreedy
Rcpp::List groupAggGreedy(Rcpp::String sd_statistic, Rcpp::NumericVector sd_threshold, Rcpp::NumericMatrix candidate_matrix, Rcpp::NumericMatrix group_matrix, Rcpp::IntegerVector nboot, Rcpp::String em_method, Rcpp::String probability_method, Rcpp::IntegerVector maximum_iterations, Rcpp::NumericVector maximum_seconds, Rcpp::NumericVector stopping_threshold, Rcpp::NumericVector log_stopping_threshold, Rcpp::LogicalVector verbose, Rcpp::IntegerVector step_size, Rcpp::IntegerVector samples, Rcpp::String monte_method, Rcpp::NumericVector monte_error, Rcpp::IntegerVector monte_iter);
RcppExport SEXP _fastei_groupAggGreedy(SEXP sd_statisticSEXP, SEXP sd_thresholdSEXP, SEXP candidate_matrixSEXP, SEXP group_matrixSEXP, SEXP nbootSEXP, SEXP em_methodSEXP, SEXP probability_methodSEXP, SEXP maximum_iterationsSEXP, SEXP maximum_secondsSEXP, SEXP stopping_thresholdSEXP, SEXP log_stopping_thresholdSEXP, SEXP verboseSEXP, SEXP step_sizeSEXP, SEXP samplesSEXP, SEXP monte_methodSEXP, SEXP monte_errorSEXP, SEXP monte_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::String >::type sd_statistic(sd_statisticSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sd_threshold(sd_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type candidate_matrix(candidate_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type group_matrix(group_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nboot(nbootSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type em_method(em_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type probability_method(probability_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type maximum_iterations(maximum_iterationsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type maximum_seconds(maximum_secondsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type stopping_threshold(stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type log_stopping_threshold(log_stopping_thresholdSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type monte_method(monte_methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type monte_error(monte_errorSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type monte_iter(monte_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(groupAggGreedy(sd_statistic, sd_threshold, candidate_matrix, group_matrix, nboot, em_method, probability_method, maximum_iterations, maximum_seconds, stopping_threshold, log_stopping_threshold, verbose, step_size, samples, monte_method, monte_error, monte_iter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastei_RsetParameters", (DL_FUNC) &_fastei_RsetParameters, 2},
    {"_fastei_EMAlgorithmFull", (DL_FUNC) &_fastei_EMAlgorithmFull, 12},
    {"_fastei_bootstrapAlg", (DL_FUNC) &_fastei_bootstrapAlg, 15},
    {"_fastei_groupAgg", (DL_FUNC) &_fastei_groupAgg, 17},
    {"_fastei_groupAggGreedy", (DL_FUNC) &_fastei_groupAggGreedy, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastei(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
