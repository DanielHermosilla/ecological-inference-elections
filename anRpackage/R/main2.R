# dyn.load("src/infPackage.so")
# library(infPackage)
# library(jsonlite)
# library(R6)


# TODO: ADD A SUMMARY METHOD

#' rxg: An R6 Class for Running an Expectation-Maximization Algorithm
#'
#' This class implements an EM algorithm using different methods for approximating the E-step such as "Multinomial",
#' 		"Hit and Run", "MVN CDF", "MVN PDF", and "Exact".
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
rxg <- R6Class("rxg",
    public = list(

        #' @field X A (b x c) matrix with the observed results of the candidate votes (c) on a given
        #' 		ballot box (b). Provided manually or loaded from JSON.
        X = NULL,

        #' @field W A (b x g) matrix with the observed results of the demographical group votes (g) on
        #' 		a given ballot box (b). Provided manually or loaded from JSON.
        W = NULL,

        #' @field method A string indicating the EM method. One of: "Multinomial", "Hit and Run", "MVN CDF",
        #' 		"MVN PDF", "Exact".
        method = NULL,

        #' @field probability A (g x c) matrix that would store the final estimated probabilities of having
        #' 		a given group (g) voting for a candidate (c).
        probability = NULL,

        #' @field logLikelihood A numeric vector that will store the log-likelihood values among the total
        #' 		iterations of the Expected Maximization algorithm.
        logLikelihood = NULL,

        #' @field total_iterations An integer indicating the total iterations that the Expected Maximization
        #' 		algorithm did.
        total_iterations = NULL,

        #' @field total_time The time that the EM algorithm took.
        total_time = NULL,

        #' Constructor for the EMModel class.
        #'
        #' @param X A numeric matrix (c x b) of observed candidate votes (optional, required if `jsonPath` is NULL).
        #' @param W A numeric matrix (b x g) of observed demographic group votes (optional, required if `jsonPath` is NULL).
        #' @param jsonPath A string containing a path to a JSON file with `"X"` and `"W"` matrices. (optional, required if `X` or `W` are NULL)
        #' @return An initialized EMModel object.
        #' @examples
        #' model <- EMModel$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3), method = "Multinomial")
        initialize = function(X = NULL, W = NULL, jsonPath = NULL) {
            # Verify if the JSON path is valid
            if (!is.null(jsonPath) && nzchar(jsonPath)) {
                data <- jsonlite::fromJSON(jsonPath)
                if (!all(c("X", "W") %in% names(data))) stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")

                self$X <- t(as.matrix(data$X)) # ERASE THE TRANSPOSE LATER
                self$W <- as.matrix(data$W)
            } else { # Case when there's no JSON path; hand matrices
                if (is.null(X) || is.null(W)) stop("Either provide X and W or an non-empty JSON path")

                self$X <- as.matrix(X)
                self$W <- as.matrix(W)
            }

            # ---- C initialization ---- #

            # Convert the parameters to C, must be double in C:
            param_X <- matrix(as.numeric(self$X), nrow(self$X), ncol(self$X))
            param_W <- matrix(as.numeric(self$W), nrow(self$W), ncol(self$W))

            # ---- Check if there's a candidate without votes ---- #
            # Gets an integer vector with the index of rows that sum 0.
            private$empty_candidates <- which(rowSums(self$X == 0) == ncol(self$X))
            if (length(private$empty_candidates) != 0) param_X <- param_X[-private$empty_candidates, drop = FALSE]

            RsetParameters(t(param_X), param_W) # The C code uses the X matrix as (c x b)
        },

        #' @name rxg$compute
        #' Precompute iteration-independent variables that can be reused for optimizing the algorithm.
        #'
        #' The Hit and Run and Exact method offers some calculations that doesn't depend on the current EM iteration, nor
        #' any unseen probability. Hence, calling this method would precompute the reusable values before calling the main
        #' computation. The main computation would run the precomputation if this method haven't been called.
        #'
        #' @param method A string containing the method to precompute. Options: "Exact", "Hit and Run"
        #' @param ... Additional arguments required by specific methods:
        #'   \itemize{
        #'     \item **"Hit and Run" Method:**
        #'       \itemize{
        #'         \item `step_size` (Integer): The step size (`M`) for the Hit and Run algorithm.
        #'         \item `samples` (Integer): The number of samples (`S`) to generate.
        #'       }
        #'   }
        #'
        #' @note The method attribute wont be updated, since the main computation haven't been called and it's possible (but weird) that the
        #' user may want to run another method even if the precomputation was made.
        #'
        #' @return Updates are made on the C internal memory.
        precompute = function(method, ...) {
            params <- list(...)

            if (method == "Hit and Run") {
                # ========= HIT AND RUN ========= #
                # ------ Validation check ------- #
                # Check for the given step size or samples. If it's not provided, return an error

                if (is.null(params$step_size) || params$step_size < 0) stop("precompute():\tA valid 'step_size' wasn't provided")
                if (is.null(params$samples) || params$samples < 0) stop("precompute():\tA valid 'samples' wasn't provided")

                # If the EM parameters differ from the precomputed values => erase the precomputation
                if (private$been_precomputed_hr &&
                    (private$hr_step_size != params$step_size ||
                        private$hr_samples != params$hr_samples)) {
                    clean_hr_precompute()
                }

                private$hr_step_size <- params$step_size
                private$hr_samples <- params$hr_samples
                # ------------- .... ------------ #

                message("Precomputing the Hit and Run method")
                RprecomputeHR(as.integer(private$hr_samples), as.integer(private$hr_step_size))
                # Update output variables
                been_precomputed_hr <- TRUE
                # self$method <- "Hit and Run"
            } else if (method == "Exact") {
                # ========= EXACT ========= #
                message("Precomputing the Exact method")
                RprecomputeExact()
                been_computed_exact <- TRUE
                # Update inner variables
                # self$method <- "Exact"
                # ---------- ... ---------- #
            } else {
                stop("Invalid method for precomputing. Must be either Hit and Run or Exact")
            }
            # NOTE:
            # It won't update the method, since the main computation haven't been called and it's possible (but weird) that the
            # user may want to run another method even if the precomputation was made.

            invisible(self)
        },

        #' Runs the EM algorithm and stores the results.
        #'
        #' @description Executes the Expectation-Maximization (EM) algorithm based on the selected method.
        #' Additional parameters may be required depending on the method.
        #'
        #' @param string main_method The method for estimating the Expectation-Maximization (EM) algorithm. (default: "Multinomial")
        #' @param string probability_method The method for obtaining the initial probability. Options:
        #' 		"Group proportional","Proportional", "Uniform". (default: "Group proportional")
        #' @param integer iterations The maximum amount of iterations to perform on the EM algorithm.
        #' 		(default: 1000)
        #' @param numeric stopping_threshold The minimum difference between consequent probabilities for stopping the
        #' 		iterations. (default: 0.001)
        #' @param verbose bool Boolean indicating wether to print useful messages while iterating. (default: FALSE)
        #' @param ... Additional arguments required by specific methods:
        #'   \itemize{
        #'     \item **"Hit and Run" Method:**
        #'       \itemize{
        #'         \item `step_size` (Integer): The step size (`M`) for the Hit and Run algorithm.
        #'         \item `samples` (Integer): The number of samples (`S`) to generate.
        #'       }
        #'     \item **"MVN CDF" Method:**
        #'       \itemize{
        #'         \item `multivariate_method` (Character): The integration method. Default is `"Genz2"`.
        #'         \item `multivariate_error` (Numeric): The integration error threshold. Default is `1e-6`.
        #'         \item `multivariate_iterations` (Integer): The number of Monte Carlo iterations. Default is `5000`.
        #'       }
        #'   }'
        #'
        #' @return The modified EMModel object (for method chaining).
        #'
        #' @examples
        #' model <- EMModel$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
        #' model$compute()
        #'
        #' # Example for Hit and Run method
        #' model$compute(main_method = "Hit and Run", step_size = 10, samples = 1000)
        #'
        #' # Example for MVN CDF method
        #' model$compute(main_method = "MVN CDF", multivariate_method = "Genz2", multivariate_error = 1e-6, multivariate_iterations = 5000)
        #'
        compute = function(main_method = "Multinomial",
                           probability_method = "Group proportional",
                           iterations = 1000,
                           stopping_threshold = 0.001,
                           verbose = FALSE, ...) {
            params <- list(...)

            # Check if the method provided is valid
            valid_methods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
            if (!is.character(main_method) || length(main_method) != 1 || !(main_method %in% valid_methods)) {
                stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
            }
            # Define the method
            self$method <- main_method

            if (self$method %in% c("Exact", "Multinomial", "MVN PDF")) {
                # === MULTINOMIAL | MVN PDF | EXACT === #
                # Run the EM algorithm for the Exact, Multinomial or PDF method.
                resulting_values <- EMAlgorithmAll(
                    self$method,
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose
                )

                # If the method is exact, the precomputation is going to occur if it wasn't called before
                if (self$method == "Exact") private$been_precomputed_exact <- TRUE
            } else if (self$method == "Hit and Run") {
                # ========= HIT AND RUN ========= #
                # ------ Validation check ------- #
                # Check for the given step size or samples. If it's not provided, return an error

                if (is.null(params$step_size) || params$step_size < 0) stop("compute():\tA valid 'step_size' wasn't provided")
                if (is.null(params$samples) || params$samples < 0) stop("compute():\tA valid 'samples' wasn't provided")

                # If the EM parameters differ from the precomputed values => erase the precomputation
                if (private$been_precomputed_hr &&
                    (private$hr_step_size != params$step_size ||
                        private$hr_samples != params$hr_samples)) {
                    clean_hr_precompute()
                }

                private$hr_step_size <- params$step_size
                private$hr_samples <- params$samples
                # ------------- .... ------------ #

                # Run the EM algorithm for the Hit and Run method.
                resulting_values <- EMAlgorithmHitAndRun(
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose,
                    as.integer(private$hr_step_size),
                    as.integer(private$hr_samples)
                )
                private$been_precomputed_hr <- TRUE
            } else {
                # ========= MVN CDF ========= #
                # ----- Validation check ---- #

                # Check if there's a multivariate method, otherwise, use 'Genz2' as default
                if (is.null(params$multivariate_method) || !params$multivariate_method %in% c("Genz", "Genz2")) {
                    private$mvn_method <- "Genz2"
                } else {
                    private$mvn_method <- params$multivariate_method
                }

                # Check if there's a multivariate error, otherwise, use 1e-5 as default
                if (is.null(params$multivariate_error) || params$multivariate_error < 0) {
                    private$mvn_error <- 1e-5
                } else {
                    private$mvn_error <- as.numeric(params$multivariate_error)
                }

                # Check if there's a multivariate iterations, otherwise, use 5000 as default
                if (is.null(params$multivariate_iterations) || params$mvn_iterations < 0) {
                    private$mvn_iterations <- as.integer(5000)
                } else {
                    private$mvn_iterations <- as.integer(params$mvn_iterations)
                }

                # Check if there's a multivariate error, otherwise, use 0.000001 as default
                if (is.null(params$multivariate_error) || params$multivariate_error < 0) private$mvn_error <- 0.000001
                # -------- ... --------- #

                # Run the EM algorithm for the MVN CDF method.
                resulting_values <- EMAlgorithmCDF(
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose,
                    private$mvn_method,
                    as.numeric(private$mvn_error),
                    as.integer(private$mvn_iterations)
                )
            }

            # ======= LIMIT CASE ======= #
            # --- Handle candidates ---- #
            # If there was a candidate that didn't receive any votes, add them with probability 0.
            if (length(private$empty_candidates) != 0) {
                # Create a matrix full of zeros with the final correct shape (g x c)
                # ncol(self$X) is used just to get the total amount of candidates
                new_probability <- matrix(0, nrow = nrow(resulting_values$result), ncol = ncol(self$X))

                # Identify the columns that have valid probability values
                non_zero_cols <- setdiff(seq_len(ncol(self$X)), private$empty_candidates)

                # Insert the computed probabilities into the correct places
                new_probability[, non_zero_cols] <- resulting_values$result

                self$probability <- new_probability
            } else {
                self$probability <- resulting_values$result
            }
            # ---------- ... ---------- #

            self$logLikelihood <- resulting_values$log_likelikelihood
            self$total_iterations <- resulting_values$total_iterations
            self$total_time <- resulting_values$total_time
            private$been_computed <- TRUE

            invisible(self)
        },

        #' Define the printing method
        #'
        #' @description According to the state of the algorithm (either computed or not), it prints a
        #' 		message with its most relevant parameters
        print = function() {
            cat("RxG ecological inference model\n")
            # Determine if truncation is needed
            truncated_X <- (nrow(self$X) > 5 || ncol(self$X) > 5)
            truncated_W <- (nrow(self$W) > 5 || ncol(self$W) > 5)

            cat("Candidate matrix (X):\n")
            print(self$X[1:min(5, nrow(self$X)), 1:min(5, ncol(self$X))])
            if (truncated_X) cat("...\n")

            cat("Group matrix (W):\n")
            print(self$W[1:min(5, nrow(self$W)), 1:min(5, ncol(self$W))])
            if (truncated_W) cat("...\n")

            if (private$been_computed) {
                cat("Method:\t", self$method, "\n")
                cat("Total Iterations:\t", self$totalIterations, "\n")
                cat("Total Time (s):\t", self$totalTime, "\n")
                cat("Estimated probability\t", self$probability, "\n")
            }
            invisible(self)
        },

        #' @name rxg$summary
        #'
        #' Define the summary method
        #'
        #' @description Shows, in form of a list, a selection of the most important atributes. It'll retrieve
        #' 		the method, amount of candidates, ballots and groups and the principal resuls of the EM algorithm.
        #'
        summary = function() {
            list(
                Method = self$method,
                Candidates = nrow(self$X),
                Ballots = ncol(self$X),
                Groups = ncol(self$W),
                Probabilities = if (!is.null(self$probability)) head(self$probability, 5) else "Not computed yet",
                LogLikelihood = if (!is.null(self$logLikelihood)) tail(self$logLikelihood, 5) else "Not computed yet"
            )
        },

        #' Save Model Results to a File
        #'
        #' @description Saves the current `EMModel` object to a specified file. The results can be saved in:
        #'   \itemize{
        #'     \item **RDS (Binary format)**: Preserves object structure for future use in R.
        #'     \item **JSON**: Saves model data in a human-readable format.
        #'     \item **CSV**: Saves probability matrix in a tabular format.
        #'   }
        #'
        #' @param filename Character. The file name where the results should be saved.
        #'   The file extension determines the format:
        #'   \itemize{
        #'     \item `*.rds` → Saves as **RDS** (default, binary format).
        #'     \item `*.json` → Saves as **JSON** (readable and shareable).
        #'     \item `*.csv` → Saves the **final probability matrix** as CSV.
        #'   }
        #'
        #' @return The modified `EMModel` object (for method chaining).
        #'
        #' @examples
        #' model <- EMModel$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3), method = "Multinomial")
        #' model$computeEM()
        #' model$save_results("results.rds")  # Save as RDS
        #' model$save_results("results.json") # Save as JSON
        #' model$save_results("results.csv")  # Save final probability as CSV
        save_results = function(filename) {
            # Ensure filename is valid
            if (!is.character(filename) || length(filename) != 1) {
                stop("Invalid filename. Please provide a valid file path as a character string.")
            }

            # Get file extension
            file_ext <- tools::file_ext(filename)

            # Save based on the file extension
            if (file_ext == "rds") {
                saveRDS(self, file = filename)
                message("Results saved as RDS: ", filename)
            } else if (file_ext == "json") {
                json_data <- list(
                    method = self$method,
                    probability = as.list(as.data.frame(self$finalProbability)),
                    logLikelihood = self$logLikelihood,
                    totalIterations = self$totalIterations,
                    totalTime = self$totalTime
                )
                jsonlite::write_json(json_data, filename, pretty = TRUE)
                message("Results saved as JSON: ", filename)
            } else if (file_ext == "csv") {
                if (is.null(self$finalProbability)) stop("No probability data to save in CSV format.")
                write.csv(self$finalProbability, filename, row.names = TRUE)
                message("Probability matrix saved as CSV: ", filename)
            } else {
                stop("Unsupported file format. Use '.rds', '.json', or '.csv'.")
            }

            invisible(self) # Enable method chaining
        }
    ),
    active = list(
        #' @field samples Active variable to show the Hit and Run samples if and only if self$method is "Hit and Run"
        samples = function() {
            if (self$method == "Hit and Run") {
                # It could be the case that the Hit and Run method was ran before, but later another
                # method is called, hence; only show this for the last method used.
                return(private$hr_samples)
            } else {
                return(NULL)
            }
        },

        #' @field step_size Active variable to show the Hit and Run step size if and only if self$method is "Hit and Run"
        step_size = function() {
            if (self$method == "Hit and Run") {
                # It could be the case that the Hit and Run method was ran before, but later another
                # method is called, hence; only show this for the last method used.
                return(private$hr_step_size)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_method Active variable to show the method used to estimate the Multivariate
        #' Normal CDF, if and only if self$method is "MVN CDF"
        multivariate_method = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_method)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_error Active variable to show the error threshold for the Montecarlo simulation of the
        #' Multivariate Normal CDF
        multivariate_error = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_error)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_error Active variable to show the error threshold for the Montecarlo simulation of the
        #' Multivariate Normal CDF
        multivariate_iterations = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_iterations)
            } else {
                return(NULL)
            }
        }
    ),
    private = list(
        #' @keywords internal been_called Boolean that determines if the object has been called before. It's mainly used for triggering
        #' the C cleanup and updating its global variables
        been_called = FALSE,

        #' @keywords internal been_precomputed_hr Boolean that determines if the object has been precomputed for the Hit and Run method.
        been_precomputed_hr = FALSE,

        #' @keywords internal been_precomputed_exact Boolean that determines if the object has been precomputed for the Exact method.
        been_precomputed_exact = FALSE,

        #' @keywords internal empty_candidates Integer vector with the index of candidates that didn't receive any vote.
        #' It's used for handling border cases and optimizing.
        empty_candidates = integer(0),

        #' @keywords internal been_computed Boolean that determines if the EM-algorithm have been computed.
        been_computed = FALSE,

        #' @keywords internal hr_samples The samples used on the computation of the EM-algorithm under the Hit and Run method.
        hr_samples = 0,

        #' @keywords internal hr_step_size The step size used on the computation of the EM-algorithm under the Hit and Run method.
        hr_step_size = 0,

        #' @keywords internal mvn_method Save the method used to calculate the multivariate normal cdf.
        mvn_method = "Genz2",

        #' @keywords internal mvn_error Saves the error used for estimating the Montecarlo simulation of the multivariate normal cdf.
        mvn_error = 0.000001,

        #' @keywords internal mvn_iterations Saves the maximum iterations made to estimate the Montecarlo simulation of the
        #' multivariate normal cdf.
        mvn_iterations = 5000,

        #' Destructor to clean allocated memory
        #'
        #' When the object gets removed or garbage collected, it cleans the allocated memory on C that could had been
        #' reused for the object instance. It's used for avoiding memory leaks due to the C behavior.
        #'
        #' @note For having the finalizer on the private access, R6 must be >= 2.4.0.
        finalize = function() {
            if (private$been_precomputed_exact) {
                clean_exact_precompute()
            }
            if (private$been_precomputed_hr) {
                clean_hr_precompute()
            }
            if (private$been_called) {
                clean_everything()
            }
        }
    )
)

#' Summarize the object main attributes.
#'
#' This method is an S3 styled wrapper around `summary()`, allowing to get, as a list, the
#' main details of the object. Additional parameters are shown depending on the state of
#' the object (computed/uncomputed).
#'
#' @inheritParams rxg$compute  # Inherit params from compute()
#' @export
#' @method summary rxg
summary.rxg <- function(object, ...) {
    object$summary()
}

#' Predict Probability using an EM Model
#'
#' This method is a wrapper around `compute()`, allowing prediction of probabilities
#' based on specified parameters. The results are stored in `object$probability`. Refer
#' to object$compute method for more information about its function and parameters. However
#' it will esentially trigger the EM algorithm.
#'
#' @inheritParams rxg$compute  # Inherit params from compute()
#' @return A matrix of estimated probabilities.
#' @export
#' @method predict rxg
predict.rxg <- function(object, ...) {
    params <- list(...)
    do.call(object$compute, params) # Calls compute() with the right arguments
    return(object$probability)
}

#' Convert EMModel Object to a Probability Matrix
#'
#' Extracts the probability matrix from the model, making it able to manipulate it as
#' a matrix
#'
#' @inheritParams EMModel$compute  # Inherit docs from compute()
#' @return A matrix containing the estimated probabilities.
#' @export
#' @method as.matrix rxg
as.matrix.rxg <- function(object, ...) {
    if (is.null(object$probability)) {
        stop(paste0("Probability matrix not available. Run compute() or ", object, "$compute() first."))
    }
    return(object$probability)
}

#' Update an existing RxG model with a new EM algorithm computation
#'
#' This function updates an object with a new Expected Maximization computation
#' with other parameters.
#'
#' @param object An `rxg` object.
#' @param ... New parameters to pass to `compute()`.
#' @return The updated object.
#' @export
#' @method update rxg
update.rxg <- function(object, ...) {
    params <- list(...)

    # Ensure the model has been computed before updating
    if (is.null(object$probability)) {
        stop(paste0("Model must be computed before updating. Run compute() or ", object, "$compute() first."))
    }

    # Recompute with new parameters
    do.call(object$compute, params)

    return(object)
}


#' Returns the current amount of candidates.
#'
#' Given a initialized object, it returns the amount of candidates that it has. It's equivalent
#' of running ncol(object$X).
#'
#' @param object An `rxg` object.
#' @return The amount of candidates.
#' @export
#' @method candidates rxg
candidates.rxg <- function(object) {
    # Ensure the model has been computed before updating
    if (is.null(object$X)) stop("The object must be initialized.")

    return(ncol(object$X))
}


#' Returns the current amount of groups.
#'
#' Given a initialized object, it returns the amount of candidates that it has. It's equivalent
#' of running ncol(object$W).
#'
#' @param object An `rxg` object.
#' @return The amount of groups.
#' @export
#' @method grouops rxg
groups.rxg <- function(object) {
    # Ensure the model has been computed before updating
    if (is.null(object$W)) stop("The object must be initialized.")

    return(ncol(object$W))
}

#' Returns the current amount of ballots boxes.
#'
#' Given a initialized object, it returns the amount of ballot boxes it has. It's equivalent
#' of running nrow(object$W) or nrow(object$X).
#'
#' @param object An `rxg` object.
#' @return The amount of ballots.
#' @export
#' @method ballots rxg
ballots.rxg <- function(object) {
    # Ensure the model has been computed before updating
    if (is.null(object$W)) stop("The object must be initialized.")

    return(nrow(object$W))
}

#' Returns the total amount of voters in the system.
#'
#' Given a initialized object, it returns the total amount of voters. It's equivalent
#' of running sum(object$W) or sum(object$X)
#'
#' @param object An `rxg` object.
#' @return The amount of voters.
#' @export
#' @method sum rxg
sum.rxg <- function(object) {
    # Ensure the model has been computed before updating
    if (is.null(object$W)) stop("The object must be initialized.")

    return(sum(object$W))
}
