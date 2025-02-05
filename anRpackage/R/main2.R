dyn.load("src/infPackage.so")
library(jsonlite)
library(R6)

EMModel <- R6Class("ecological_inference_model",
    public = list(
        X = NULL, # Candidate matrix
        W = NULL, # Group matrix
        method = NULL, # Method used on the EM Algorithm
        probability = NULL, # Probability matrix
        logLikelihood = NULL, # Numeric array
        total_iterations = NULL, # Log-likelihood iterations
        total_time = NULL, # Time on running the EM-algorithm
        initialize = function(X = NULL, W = NULL, jsonPath = NULL, method = "Multinomial") {
            validMethods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
            if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
                stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
            }
            # Impose a given method
            self$method <- method

            # Verify if the JSON path is valid
            if (!is.null(jsonPath) && nzchar(jsonPath)) {
                data <- fromJSON(jsonPath)
                if (!all(c("X", "W") %in% names(data))) {
                    stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")
                }
                self$X <- as.matrix(data$X)
                self$W <- as.matrix(data$W)
            } else { # Case when there's no JSON path; hand matrices
                if (is.null(X) || is.null(W)) {
                    stop("Either provide X and W or an non-empty JSON path")
                }
                self$X <- as.matrix(X)
                self$W <- as.matrix(X)
            }
            # Send the parameters to C
            RsetParameters(self$X, self$W)
        },
        compute = function(probability_method = "Group proportional",
                           iterations = 1000,
                           stopping_threshold = 0.001,
                           verbose = FALSE, ...) {
            params <- list(...)

            if (self$method != "Hit and Run" || self$method != "MVN CDF") {
                # Run the EM algorithm for the Exact, Multinomial or PDF method.
                resulting_values <- EMAlgorithmAll(
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose
                )
            } else if (self$method == "Hit and Run") {
                # Check for a given step size. If it's not provided, return an error
                if (!"step_size" %in% names(params)) {
                    stop("The 'step_size' (M) wasn't provided for running the Hit and Run method.")
                }
                # Check for a given sample. If it's not provided, return an error
                if (!"samples" %in% names(params)) {
                    stop("The 'samples' (S) wasn't provided for running the Hit and run method.")
                }
                if (!is.integer(step_size) || !is.integer(samples)) {
                    stop("The 'step_size' or 'samples' are an invalid value. They must be integers.")
                }
                # Run the EM algorithm for the Hit and Run method.
                resulting_values <- EMAlgorithmHitAndRun(
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose,
                    step_size,
                    samples
                )
            } else {
                # Check if there's a multivariate method, otherwise, use 'Genz2' as default
                if (!"multivariate_method" %in% names(params)) {
                    multivariate_method <- "Genz2"
                }

                # Check if there's a multivariate error, otherwise, use '0.000001' as default
                if (!"multivariate_error" %in% names(params)) {
                    multivariate_error <- 0.000001
                }

                # Check if there's a multivariate iteration, otherwise, use '5000' as default
                if (!"multivariate_iterations" %in% names(params)) {
                    multivariate_iterations <- 5000
                }
                if (!is.integer(multivariate_iterations) || !is.numeric(multivariate_error) || !is.character(multivariate_method)) {
                    stop("Invalid types values are handed to the EM algorithm method.")
                }
                # Run the EM algorithm for the MVN CDF method.
                resulting_values <- EMAlgorithmCDF(
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose,
                    multivariate_method,
                    multivariate_error,
                    multivariate_iterations
                )
            }
            self$probability <- resulting_values$result
            self$logLikelihood <- resulting_values$log_likelikelihood
            self$total_iterations <- resulting_values$total_iterations
            self$total_time <- resulting_values$total_time
        }
    )
)
# @return Rcpp::List A list with the final probability ("result"), log-likelihood array ("log_likelihood"), total
# iterations that were made ("total_iterations") and the time that was taken ("total_time").

# Rcpp::List EMAlgorithmHitAndRun(Rcpp::String probability_method = "Group proportional",
#                                Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
#                                Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
#                                Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false),
# ¶                                Rcpp::IntegerVector step_size = Rcpp::IntegerVector::create(3000),
#                                Rcpp::IntegerVector samples = Rcpp::IntegerVector(1000));

# Rcpp::List EMAlgorithmCDF(Rcpp::String probability_method = "Group proportional",
#                          Rcpp::IntegerVector maximum_iterations = Rcpp::IntegerVector::create(1000),
#                          Rcpp::NumericVector stopping_threshold = Rcpp::NumericVector::create(0.001),
#                          Rcpp::LogicalVector verbose = Rcpp::LogicalVector::create(false),
#                          Rcpp::String multivariate_method = "Genz2",
## ¶                          Rcpp::NumericVector multivariate_epsilon = Rcpp::NumericVector::create(0.000001),
#                         Rcpp::IntegerVector multivariate_iterations = Rcpp::IntegerVector::create(5000));
