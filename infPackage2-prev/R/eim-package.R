
#' An R6 Class for estimating an eim matrix from an ecological inference problem.
#'
#' @description This class implements an EM algorithm using different methods for approximating the E-step such as "Multinomial",
#' "Hit and Run", "MVN CDF", "MVN PDF", and "Exact".
#'
#' @importFrom R6 R6Class
#' @export
eim <- R6::R6Class("eim",
    public = list(

        #' @field X (matrix) A (b x c) matrix with the observed results of the candidate votes (c) on a given
        #' ballot box (b). Provided manually or loaded from JSON.
		X = NULL,

        #' @field W (matrix) A (b x g) matrix with the observed results of the demographical group votes (g) on
        #' ballot box (b). Provided manually or loaded from JSON.
        W = NULL,

        #' @field method (character) A string indicating the EM method. One of: "Multinomial", "Hit and Run", "MVN CDF",
        #' "MVN PDF", "Exact".
        method = NULL,

        #' @field probability (matrix) A (g x c) matrix that would store the final estimated probabilities of having
        #' a given group (g) voting for a candidate (c).
        probability = NULL,

        #' @field logLikelihood (numeric) A numeric vector that will store the log-likelihood values among the total
        #' iterations of the Expected Maximization algorithm.
        logLikelihood = NULL,

        #' @field total_iterations (integer(1)) An integer indicating the total iterations that the Expected Maximization
        #' algorithm did.
        total_iterations = NULL,

        #' @field total_time (numeric(1)) The time that the EM algorithm took.
        total_time = NULL,

        #' @description Creates the object by defining the X and W matrix attributes.
		#'
        #' @param X (matrix) A matrix (c x b) of observed candidates (c) votes per ballot boxes (b) (optional, required if json_path is NULL).
		#'
        #' @param W (matrix) A matrix (b x g) of observed demographic (g) group votes per ballot boxes (b) (optional, required if json_path is NULL).
		#'
        #' @param json_path (character) A string containing a path to a JSON file with "X" and "W" matrices. (optional, required if X or W are NULL)
		#'
        #' @return An initialized eim object.
		#'
        #' @examples
        #' # Example 1: Create a eim object from a matrix
        #' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
        #'
        #' # Example 2: Create a eim object from a JSON file
        #' \dontrun{
        #' model2 <- eim$new(json_path = "a/file/path/to/a/file.json")
        #' }
        initialize = function(X = NULL, W = NULL, json_path = NULL) {

			# ======= JSON PATH READING ======= #
            # Verify if the JSON path is valid
            if (!is.null(json_path) && nzchar(json_path)) {
				# Reads the JSON
                data <- jsonlite::fromJSON(json_path)
				# Validations
				private$validate_JSON_data(data)
				private$validate_X_W(data$X, data$W)
				# Update parameters
                self$X <- as.matrix(data$X)
                self$W <- as.matrix(data$W)

			# ===== DIRECT MATRIX READING ==== #
            } else { # Case when there's no JSON path; hand matrices
				# Validations
                if (is.null(X) || is.null(W))
				stop("Either provide X and W or a valid JSON path")
				private$validate_X_W(X, W)

				# Update parameters
                self$X <- as.matrix(X)
                self$W <- as.matrix(W)
            }

            # ====== C initialization ====== #

            # Convert the parameters to C, must be double in C:
            param_X <- matrix(as.numeric(self$X), nrow(self$X), ncol(self$X))
			param_X <- t(param_X)
            param_W <- matrix(as.numeric(self$W), nrow(self$W), ncol(self$W))

            # ---- Check if there's a candidate without votes ---- #
            # Gets an integer vector with the index of rows that sum 0.
            private$empty_candidates <- which(rowSums(self$X == 0) == ncol(self$X))
            if (length(private$empty_candidates) != 0) param_X <- param_X[-private$empty_candidates, drop = FALSE]

            RsetParameters(param_X, param_W) # The C code uses the X matrix as (c x b)
        },

        #' @description Calculates the EM-independent variables for the Hit and Run and Exact method.
		#'
        #' @param method The method for precomputing. Options: "Hit and Run", "Exact"
		#'
        #' @param ... Additional arguments required by specific methods:
        #'   \itemize{
        #'     \item \strong{"Hit and Run" Method:}
        #'       \itemize{
        #'         \item \code{step_size} (integer(1)) The step size (M) for the Hit and Run algorithm. Must be a positive integer.
        #'         \item \code{samples} (integer(1)) The number of samples (S) to generate. Must be a integer
        #'       }
        #'   }
		#'
        #' @return The modified eim object (for method chaining). Updates are made on the C internal memory.
		#'
        #' @examples
        #' # Example 1: Precompute the Hit and Run method
        #' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
        #' model$precompute("Hit and Run",
        #' 					step_size = 1000,
        #' 					samples = 5000)
        #' # Changes are made on C internals API
        #'
        #' # Example 2: Precompute the Exact method
        #' model$precompute("Exact")
        precompute = function(method, ...) {
            params <- list(...)
			# Check if the method provided is valid
			private$validate_methods(method, params)

			# ========= HIT AND RUN ========= #
            if (method == "Hit and Run") {
				# Define the attributes
                private$hr_step_size <- params$step_size
                private$hr_samples <- params$samples

                message("Precomputing the Hit and Run method")
                RprecomputeHR(as.integer(private$hr_samples), as.integer(private$hr_step_size))
                # Update output variables
                private$been_precomputed_hr <- TRUE

            } 
			else if (method == "Exact") {
			    # ========= EXACT ========= #
				if (private$been_precomputed_exact) clean_exact_precompute()
                message("Precomputing the Exact method")
                RprecomputeExact()
                private$been_precomputed_exact <- TRUE
                # ---------- ... ---------- #
            } else {
                stop("Invalid method for precomputing. Must be either `Hit and Run` or `Exact`")
            }
            # NOTE:
            # It won't update the method, since the main computation haven't been called and it's possible (but weird) that the
            # user may want to run another method even if the precomputation was made.

            invisible(self)},

        #' @description Executes the Expectation-Maximization (EM) algorithm based on the selected method. Additional parameters may be required depending on the method.
		#'
        #' @param main_method (character) The method for estimating the Expectation-Maximization (EM) algorithm. Options:
        #' "Multinomial", "MVN CDF", "MVN PDF", "Hit and Run" and "Exact" (default: "Multinomial").
		#'
        #' @param probability_method (character) The method for obtaining the initial probability. Options: "Group proportional",
        #' "Proportional", "Uniform". (default: "Group proportional").
		#'
        #' @param iterations (integer(1)) The maximum amount of iterations to perform on the EM algorithm. (default: 1000).
		#'
        #' @param stopping_threshold (numeric(1)) The minimum difference between consequent probabilities for stopping the iterations.
        #' (default: 0.001).
		#'
        #' @param verbose (boolean(1)) Boolean indicating whether to print useful messages while iterating. (default: FALSE).
		#'
        #' @param ... Additional arguments required by specific methods:
        #'   \itemize{
        #'     \item \strong{"Hit and Run" Method:}
        #'       \itemize{
        #'         \item step_size (integer(1)): The step size (M) for the Hit and Run algorithm. \cr
        #'         \item samples (integer(1)): The number of samples (S) to generate.
        #'       }
        #'     \item \strong{"MVN CDF" Method:}
        #'       \itemize{
        #'         \item multivariate_method (character): The integration method. Default is "Genz2". \cr
        #'         \item multivariate_error (numeric(1)): The integration error threshold. Default is 1e-6. \cr
        #'         \item multivariate_iterations (integer(1)): The number of Monte Carlo iterations. Default is 5000.
        #'       }
        #'   }
		#'
        #' @return The modified eim object (for method chaining).
        #'
        #' @examples
        #' # Example 1: Compute the Expected Maximization with default values
        #' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
        #' model$compute() # Retrieves the object with updated attributes
        #'
        #' # Example 2: Compute the Expected Maximization for the Hit and Run method
        #' \donttest{
		#' 	model$compute(main_method = "Hit and Run", 
		#' 					step_size = 3000, 
		#' 					samples = 1000)
        #' 	}
        #'
        #' # Example 3: Omit arguments to the Hit and Run method
        #' \dontrun{
        #' 	model$compute(main_method = "Hit and Run", 
		#' 					step_size = 3000) 
		#' 	# Error; must hand in samples parameter too
        #' 	}
        #'
        #' # Example 4: Run the MVN CDF with default values
        #' model$compute(main_method = "MVN CDF")
        #'
        #' # Example 5: Run a Exact estimation with user defined parameters
        #' model$compute(main_method = "Exact", 
		#' 				probability_method = "Uniform", 
		#'  			iterations = 5, 
		#' 				stopping_threshold = 1e-3) 
		#' # Verbose was omitted
        compute = function(main_method = "Multinomial",
                           probability_method = "Group proportional",
                           iterations = 1000,
                           stopping_threshold = 0.001,
                           verbose = FALSE, ...) {
            params <- list(...)

			# Check if the method provided is valid
			private$validate_methods(main_method, params)
			# Define the method
            self$method <- main_method

			# ========= HIT AND RUN ========= #
			if (self$method == "Hit and Run") {
			
				# Define the attributes
                private$hr_step_size <- params$step_size
                private$hr_samples <- params$samples

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

			}
			# ========= MVN CDF ========= #
			else if (self$method == "MVN CDF") {

                # Check if there's a multivariate method, otherwise, use 'Genz2' as default
                if (is.null(params$multivariate_method)) {
                    private$mvn_method <- "Genz2"
                } else {
                    private$mvn_method <- params$multivariate_method
                }

                # Check if there's a multivariate error, otherwise, use 1e-5 as default
                if (is.null(params$multivariate_error)) {
                    private$mvn_error <- 1e-5
                } else {
                    private$mvn_error <- as.numeric(params$multivariate_error)
                }

                # Check if there's a multivariate iterations, otherwise, use 5000 as default
                if (is.null(params$multivariate_iterations)) {
                    private$mvn_iterations <- as.integer(5000)
                } else {
                    private$mvn_iterations <- as.integer(params$mvn_iterations)
                }

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
			# ==== MULTINOMIAL | MVN PDF | EXACT ==== #
			else {
				# Run the EM algorithm
				resulting_values <- EMAlgorithmAll(
                    self$method,
                    probability_method,
                    iterations,
                    stopping_threshold,
                    verbose
                )

                # If the method is exact, the precomputation is going to occur if it wasn't called before
                if (self$method == "Exact") private$been_precomputed_exact <- TRUE
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

                self$probability <- as.matrix(new_probability)
            } else {
                self$probability <- as.matrix(resulting_values$result)
            }
            # ---------- ... ---------- #

            self$logLikelihood <- as.numeric(resulting_values$log_likelihood)
            self$total_iterations <- as.integer(resulting_values$total_iterations)
            self$total_time <- resulting_values$total_time
            private$been_computed <- TRUE

            invisible(self)
        },

        #' @description According to the state of the algorithm (either computed or not), it prints a message with its most relevant parameters
		#'
        #' @return The own object for method chaining.
		#'
        #' @examples
        #' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
        #' print(model) # Will print the X and W matrix.
        #'
        #' model$compute()
        #' print(model) # Will print the Xand W matrix among the EM results.
        print = function() {
            cat("eim ecological inference model\n")
            # Determine if truncation is needed
            truncated_X <- (nrow(self$X) > 5 || ncol(self$X) > 5)
            truncated_W <- (nrow(self$W) > 5 || ncol(self$W) > 5)

            cat("Candidate matrix (X) [b x c]:\n")
            print(self$X[1:min(5, nrow(self$X)), 1:min(5, ncol(self$X))])
            if (truncated_X) cat("...\n")

            cat("Group matrix (W) [b x g]:\n")
            print(self$W[1:min(5, nrow(self$W)), 1:min(5, ncol(self$W))])
            if (truncated_W) cat("...\n")

            if (private$been_computed) {
                cat("Estimated probability [g x c]:\n")
                print(self$probability[1:min(5, nrow(self$probability)), 1:min(5, ncol(self$probability))])
                cat("Method:\t", self$method, "\n")
                cat("Total Iterations:\t", self$total_iterations, "\n")
                cat("Total Time (s):\t", self$total_time, "\n")
				cat("Log-likelihood:\t", tail(self$logLikelihood, 1))
            }
            invisible(self)
        },

        #' @description Shows, in form of a list, a selection of the most important attributes. It'll retrieve the method, amount of candidates, ballots and groups and the principal results of the EM algorithm.
		#'
        #' @return (list) A list with the method, candidates, ballots, group, probabilities and log-likelihood.
		#'
        #' @examples
        #' model <- eim$new(X = matrix(1:15, 5, 3), W = matrix(1:10, 2, 5))
        #' a_list <- model$summary()
        #' a_list$method # Not computed yet
        #' a_list$groups # 2
        #' a_list$ballots # 5
        #' names(a_list) 
		#' # "candidates" "groups" "ballots" "method" "probabilities" "logLikelihood"
        summary = function() {
            list(
                candidates = nrow(self$X),
                groups = ncol(self$W),
                ballots = ncol(self$X),
                method = if (!is.null(self$method)) self$method else "Not computed yet",
                probabilities = if (!is.null(self$probability)) self$probability else "Not computed yet",
                logLikelihood = if (!is.null(self$logLikelihood)) self$logLikelihood else "Not computed yet"
            )
        },

        #' @description Saves the current eim object to a specified file. The results can be saved in:
        #'   \itemize{
        #'     \item \strong{RDS (Binary format):} Preserves object structure for future use in R. \cr
        #'     \item \strong{JSON}: Saves model data in a human-readable format. \cr
        #'     \item \strong{CSV}: Saves probability matrix in a tabular format.
        #'   }
		#'
        #' @param filename (character) The file name where the results should be saved with its extension.
        #'   The file extension determines the format:
        #'   \itemize{
        #'     \item \var{*.rds} → Saves as \file{RDS} (default, binary format). \cr
        #'     \item \var{*.json} → Saves as \file{JSON} (readable and shareable). \cr
        #'     \item \var{*.csv} → Saves the \emph{final probability matrix} as \file{CSV}.
        #'   }
		#'
        #' @return The modified eim object (for method chaining).
		#'
        #' @examples
        #' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
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
                    probability = as.list(as.data.frame(self$probability)),
                    logLikelihood = self$logLikelihood,
                    totalIterations = self$total_iterations,
                    totalTime = self$total_time
                )
                jsonlite::write_json(json_data, filename, pretty = TRUE)
                message("Results saved as JSON: ", filename)
            } else if (file_ext == "csv") {
                if (is.null(self$probability)) stop("No probability data to save in CSV format.")
                write.csv(self$probability, filename, row.names = TRUE)
                message("Probability matrix saved as CSV: ", filename)
            } else {
                stop("Unsupported file format. Use '.rds', '.json', or '.csv'.")
            }

            invisible(self) # Enable method chaining
        }
    ),
    active = list(
        #' @field samples (integer(1)) Active variable to show the Hit and Run samples if and only if the method is "Hit and Run"
        samples = function() {
            if (self$method == "Hit and Run") {
                return(private$hr_samples)
            } else {
                return(NULL)
            }
        },

        #' @field step_size (integer(1)) Active variable to show the Hit and Run step size if and only if the method is "Hit and Run"
        step_size = function() {
            if (self$method == "Hit and Run") {
                return(private$hr_step_size)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_method (character) Active variable to show the method used to estimate the Multivariate
        #' Normal CDF, if and only if self$method is "MVN CDF"
        multivariate_method = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_method)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_error (numeric(1)) Active variable to show the error threshold for the Montecarlo simulation of the Multivariate Normal CDF
        multivariate_error = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_error)
            } else {
                return(NULL)
            }
        },

        #' @field multivariate_iterations (numeric(1)) Active variable to show the number of iterations for the Montecarlo simulation of the Multivariate Normal CDF
        multivariate_iterations = function() {
            if (self$method == "MVN CDF") {
                return(private$mvn_iterations)
            } else {
                return(NULL)
            }
        }
    ),
    private = list(
        # PRIVATE FIELD: been_called Boolean that determines if the object has been called before (used for C cleanup).
        been_called = FALSE,
        # PRIVATE FIELD been_precomputed_hr Boolean that indicates if precomputation for the Hit and Run method has been performed.
        been_precomputed_hr = FALSE,
        # PRIVATE FIELD been_precomputed_exact Boolean that indicates if precomputation for the Exact method has been performed.
        been_precomputed_exact = FALSE,
        # PRIVATE FIELD empty_candidates Integer vector with indices of candidates that received no votes.
        empty_candidates = integer(0),
        # PRIVATE FIELD been_computed Boolean that indicates if the EM algorithm has been computed.
        been_computed = FALSE,
        # PRIVATE FIELD hr_samples The number of samples used in the Hit and Run method.
        hr_samples = 0,
        # PRIVATE FIELD hr_step_size The step size used in the Hit and Run method.
        hr_step_size = 0,
        # PRIVATE FIELD mvn_method The integration method used for the MVN CDF.
        mvn_method = "Genz2",
        # PRIVATE FIELD mvn_error The error threshold used for the MVN CDF Monte Carlo simulation.
        mvn_error = 0.000001,
        # PRIVATE FIELD mvn_iterations The number of iterations for the MVN CDF Monte Carlo simulation.
        mvn_iterations = 5000,
		# PRIVATE METHOD validate_JSON_data Validates an input JSON file data.
		validate_JSON_data = function(data) {
			# Check that the keys "X" and "W" exist
			if (!all(c("X", "W") %in% names(data))) {
				stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")
			}
			# Check that they are not NULL
			if (is.null(data$X) || is.null(data$W)) {
				stop("'X' and 'W' cannot be NULL")
			}
			# Check that they are not empty
			if (length(data$X) == 0 || length(data$W) == 0) {
				stop("'X' and 'W' must not be empty")
			}
			# Check for missing values
			if (any(is.na(data$X)) || any(is.na(data$W))) {
				stop("Matrices 'X' and 'W' cannot contain missing values")
			}
			# Check that both matrices have the same number of rows (i.e., ballot boxes)
			if (nrow(data$X) != nrow(data$W)) {
				stop("The ballot boxes dimensions are different between 'X' and 'W'")
			}
		},
		# PRIVATE METHOD validate_X_W Validates an X or W matrix.
		validate_X_W = function(X, W) {
			# Convert to matrices if they are not already
			X <- as.matrix(X)
			W <- as.matrix(W)

			# Check that both matrices have the same number of rows
			if (nrow(X) != nrow(W)) {
				stop("The ballot boxes dimensions are different between 'X' and 'W'")
			}
			if (ncol(X) <= 1) {
				stop("The candidate matrix dimension must be 2 or greater")
			}
			if (ncol(W) <= 1) {
				stop("The group matrix dimension must be 2 or greater")
			}
		},
		# PRIVATE METHOD validate_methods Checks the EM algorithm method input. Params would be a list of optional inputs.
		validate_methods = function(method, params) {
		    
		    # ====== METHOD INPUT ====== #
		    valid_methods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
		    
		    if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
		        stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
		    }

			# ====== HIT AND RUN ======= #
			if (method == "Hit and Run") {
                tolerance <- .Machine$double.eps^0.5
				# Check for the step size
				if (is.null(params$step_size) 
					|| params$step_size < 0
					|| abs(params$step_size %% 1) > tolerance) # Last condition: check if it has decimal part
				stop("Hit and Run:\tA valid 'step_size' wasn't provided")

				# Check for samples
                if (is.null(params$samples) 
					|| params$samples < 0
				    || abs(params$step_size %% 1) > tolerance) # Last condition: check if it has decimal part
				stop("Hit and Run:\tA valid 'samples' wasn't provided")

				# If the EM parameters differ from the precomputed values => erase the precomputation
                if (private$been_precomputed_hr &&
                    (private$hr_step_size != params$step_size ||
                        private$hr_samples != params$hr_samples)) {
					message("Warning:\tRunning the Hit and Run method with different values than the precomputed ones")
                    clean_hr_precompute()
                }

			}
			# ========= MVN CDF ========= #
			else if (method == "MVN CDF") {

				valid_cdf_methods <- c("Genz, Genz2")
				# If there's a parameter for MVN CDF; it has to be valid (could be omitted) 
				if (!is.null(params$multivariate_method) && !(params$multivariate_method %in% valid_cdf_methods))
					stop("MVN CDF\t: The method for estimating the CDF isn't supported. Must be one of: ", 
						paste(valid_cdf_methods, collapse = ", "))

				# If there's a multivariate error, it has to be valid.
				if (!is.null(params$multivariate_error) && params$multivariate_error < 0)
					stop("MVN CDF\t: The error threshold for the Montecarlo simulation isn't valid")
				
				# Check if there's a multivariate iterations, otherwise, use 5000 as default
				tolerance <- .Machine$double.eps^0.5
				if (!is.null(params$multivariate_iterations) && 
						(params$mvn_iterations < 0 || abs(params$mvn_iterations %% 1) > tolerance))
					stop("MVN CDF\t: An invalid iteration parameter was handed")
			
			}

		},
		# PRIVATE METHOD finalize Destructor to clean allocated memory via C routines.
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
#' Return a list with the principal attributes
#'
#' @description This method is an S3 styled wrapper around summary(), allowing to get, as a list,
#' the main details of the object. Additional parameters are shown depending on the state of
#' the object (computed/uncomputed).
#'
#' @param A eim object
#'
#' @return (list) A list with the principal attributes
#'
#' @examples
#' model <- eim$new(X = matrix(1:15, 5, 3), W = matrix(1:10, 2, 5))
#' a_list <- summary(model)
#' a_list$method # Not computed yet
#' a_list$groups # 2
#' a_list$ballots # 5
#' names(a_list) # "candidates" "groups" "ballots" "method" "probabilities" "logLikelihood"
#'
#' @export
summary.eim <- function(object) {
    object$summary()
}

#' Predict the R x C probability using an EM Model
#'
#' Given the parameters of an Expected Maximization, computes and updates the results of the
#' input algorithm.
#'
#' @description This method is a wrapper around compute(), allowing prediction of probabilities
#' based on specified parameters. The results are stored in object$probability. Refer
#' to object$compute method for more information about its function and parameters. However
#' it will essentially trigger the EM algorithm.
#'
#' @param object An eim object
#'
#' @param ... Additional arguments required by specific methods:
#'   \itemize{
#' 		\item \var{main_method} (character) The method for estimating the Expectation-Maximization (EM) algorithm. Options:
#' 			"Multinomial", "MVN CDF", "MVN PDF", "Hit and Run" and "Exact" (default: "Multinomial"). \cr
#'
#' 		\item \var{probability_method} (character) The method for obtaining the initial probability. Options:
#' 			"Group proportional", "Proportional", "Uniform". (default: "Group proportional"). \cr
#'
#' 		\item \var{iterations} (integer(1)) The maximum amount of iterations to perform on the EM algorithm. (default: 1000). \cr
#'
#' 		\item \var{stopping_threshold} (numeric(1)) The minimum difference between consequent probabilities for stopping the
#' 		 iterations (default: 0.001). \cr
#'
#' 		\item \var{verbose} (boolean(1)) Boolean indicating whether to print useful messages while iterating. (default: FALSE). \cr
#'
#'     \item \strong{"Hit and Run" Method:}
#'       \itemize{
#'         \item \verb{step_size} (integer(1)): The step size (M) for the Hit and Run algorithm. \cr
#'         \item \verb{samples} (integer(1)): The number of samples (S) to generate. \cr
#'       }
#'     \item \strong{"MVN CDF" Method:}
#'       \itemize{
#'         \item \verb{multivariate_method} (character): The integration method. Default is "Genz2". \cr
#'         \item \verb{multivariate_error} (numeric(1)): The integration error threshold. Default is 1e-6. \cr
#'         \item \verb{multivariate_iterations} (integer(1)): The number of Monte Carlo iterations. Default is 5000. \cr
#'       }
#'   }
#'
#'
#' @return (matrix) A matrix of estimated probabilities.
#'
#' @examples
#' # Example 1: Predict the Hit and Run method
#' model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#' predict(model, "Hit and Run", 
#' 				step_size = 1000, 
#' 				samples = 5000) 
#' Returns the matrix of probabilities
#' nrow(model) # 3
#'
#' # Example 2:  Predict the Exact method
#' predict(model, "Exact")
#'
#' @export
predict.eim <- function(object, ...) {
    params <- list(...)
    do.call(object$compute, params) # Calls compute() with the right arguments
    return(object$probability)
}

#' Returns the eim's model probability matrix.
#'
#' Extracts the probability matrix from the model, making it able to manipulate it as
#' a matrix
#'
#' @param object
#'
#' @return (matrix) A matrix containing the estimated probabilities.
#'
#' @examples
#' model <- eim$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' predict(model)
#' ncol(as.matrix(model)) # Returns the column of probabilities
#'
#' @export
as.matrix.eim <- function(object) {
    if (is.null(object$probability)) {
        stop(paste0("Probability matrix not available. Run compute() or ", object, "$compute() first."))
    }
    return(object$probability)
}

#' Update an existing eim model with a new EM algorithm computation
#'
#' This function updates an object with a new Expected Maximization computation
#' with other parameters.
#'
#' @param object An eim object.
#'
#' @param ... New parameters to pass to compute().
#'
#' @return Changes made on the updated object.
#'
#' @examples
#' model <- eim$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' predict(model)
#' update(model, "MVN PDF")
#' model$method # Returns 'MVN PDF'
#'
#' @export
update.eim <- function(object, ...) {
    params <- list(...)

    # Ensure the model has been computed before updating
    if (is.null(object$probability)) {
        stop(paste0("Model must be computed before updating. Run compute() or ", object, "$compute() first."))
    }

    # Recompute with new parameters
    do.call(object$compute, params)

    return(object)
}

#' Returns the total amount of voters in the system.
#'
#' Given a initialized object, it returns the total amount of voters. It's equivalent
#' of running sum(object$W) or sum(object$X)
#'
#' @param object An eim object.
#'
#' @return The amount of voters.
#'
#' @examples
#' model <- eim$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' sum(model) # 1395
#'
#' @export
sum.eim <- function(object) {
    # Ensure the model has been computed before updating
    if (is.null(object$W)) stop("The object must be initialized.")

    return(sum(object$W))
}
