dyn.load("src/infPackage.so")
library(jsonlite)
library(R6)

EMModel <- R6Class("ecological_inference_model",
    public = list(

        #' @field X A (c x b) matrix with the observed results of the candidate votes (c) on a given ballot box (b). Provided manually or loaded from JSON.
        X = NULL,

        #' @field W A (b x g) matrix with the observed results of the demographical group votes (g) on a given ballot box (b). Provided manually or loaded from JSON.
        W = NULL,

        #' @field method A string indicating the EM method. One of: "Multinomial", "Hit and Run", "MVN CDF", "MVN PDF", "Exact".
        method = NULL,

        #' @field probability A (c x g) matrix that would store the final estimated probabilities of having a given group (g) voting for a candidate (c).
        probability = NULL, # Probability matrix

        #' @field logLikelihood A numeric vector that will store the log-likelihood values among the total iterations of the Expected Maximization algorithm.
        logLikelihood = NULL,

        #' @field total_iterations An integer indicating the total iterations that the Expected Maximization algorithm did.
        total_iterations = NULL,

        #' @field total_time The time that the EM algorithm took.
        total_time = NULL,

        #' Constructor for the EMModel class.
        #'
        #' @param X A numeric matrix (c x b) of observed candidate votes (optional, required if `jsonPath` is NULL).
        #' @param W A numeric matrix (b x g) of observed demographic group votes (optional, required if `jsonPath` is NULL).
        #' @param jsonPath A string containing a path to a JSON file with `"X"` and `"W"` matrices. (optional, required if `X` or `W` are NULL)
        #' @param method A string indicating the EM method to use. Defaults to `"Multinomial"`.
        #' @return An initialized EMModel object.
        #' @examples
        #' model <- EMModel$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3), method = "Multinomial")
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

        # Method: compute the Expected Maximization algorithm
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

            invisible(self)
        },

        # Method: store the results to a file
        save_results = function(filename) {
            saveRDS(self, file = filename)
            message("Results saved to ", filename)
            invisible(self)
        }
    )
)
