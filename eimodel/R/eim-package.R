#' An R6 Class for estimating an probability matrix from an R x C ecological inference problem.
#'
#' @name eim
#'
#' @author Charles Thraves, Pablo Ubilla
#'
#' @references Thraves, C. and Ubilla, P.: 'Fast Ecological Inference Algorithm for the R×C Case'
#'
#' @description This class implements an EM algorithm using different methods for approximating the E-step such as \code{"Multinomial",
#' "Hit and Run", "MVN CDF", "MVN PDF"}, and \code{"Exact"}.
#'
#' @importFrom R6 R6Class
#' @export
eim <- R6::R6Class("eim",
  public = list(

    #' @field X \emph{(matrix)} An (b x c) matrix with the observed results of the candidate votes (c) on a given
    #' ballot box (b). Provided manually or loaded from JSON.
    X = NULL,

    #' @field W \emph{(matrix)} An (b x g) matrix with the observed results of the demographical group votes (g) on
    #' a given ballot box (b). Provided manually or loaded from JSON.
    W = NULL,

    #' @field method \emph{(character)} A string indicating the EM method. One of: \code{"Multinomial", "Hit and Run", "MVN CDF",
    #' "MVN PDF"} or \code{"Exact"}.
    method = NULL,

    #' @field probability \emph{(matrix)} An (g x c) matrix that would store the final estimated probabilities of having
    #' a given group (g) voting for a candidate (c).
    probability = NULL,

    #' @field logLikelihood \emph{(numeric)} An numeric vector that will store the log-likelihood values among the total
    #' iterations of the Expected Maximization algorithm.
    logLikelihood = NULL,

    #' @field total_iterations \emph{(integer(1))} An integer indicating the total iterations that the Expected Maximization
    #' algorithm did.
    total_iterations = NULL,

    #' @field total_time \emph{(numeric(1))} The time that the EM algorithm took.
    total_time = NULL,

    #' @field finish_state \emph{Character} The reason of the algorithm's finishing. It could either be "Maximum iterations", "Log-likelihood decrease", "Convergence" or "Early exit".
    finish_state = NULL,

    #' @description Creates the object by defining the X and W matrix attributes.
    #'
    #' @param X \emph{(matrix)} A matrix (c x b) of observed candidates (c) votes per ballot boxes (b) (optional, required if json_path is NULL).
    #'
    #' @param W \emph{(matrix)} A matrix (b x g) of observed demographic (g) group votes per ballot boxes (b) (optional, required if json_path is NULL).
    #'
    #' @param json_path \emph{(character)} A string containing a path to a JSON file with "X" and "W" matrices. (optional, required if X or W are NULL)
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
        if (is.null(X) || is.null(W)) {
          stop("Either provide X and W or a valid JSON path")
        }
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

    #' @description Calculates the EM-independent variables for the \code{Hit and Run} and \code{Exact} method.
    #'
    #' @param method \emph{(character)} The method for precomputing. Options: \code{"Hit and Run"} or \code{"Exact"}.
    #'
    #' @param ... Additional arguments required by specific methods:
    #'   \itemize{
    #'     \item \strong{"Hit and Run" Method:}
    #'       \itemize{
    #'         \item \code{step_size} \emph{(integer(1))} The step size (M) for the Hit and Run algorithm. Must be a positive integer.
    #'         \item \code{samples} \emph{(integer(1))} The number of samples (S) to generate. Must be a integer
    #'       }
    #'   }
    #'
    #' @note Precomputing would eventually accelerate the Expected Maximization algorithm computation time.
    #'
    #' @return The modified eim object (for method chaining). Updates are made on the C internal memory.
    #'
    #' @examples
    #' # Example 1: Precompute the Hit and Run method
    #' simulations <- simulate_elections(num_ballots = 20,
    #' 				num_candidates = 5,
    #' 				num_groups = 3,
    #' 				ballot_voters = rep(100, 20))
    #'
    #' model <- eim$new(simulations$X, simulations$W)
    #'
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
      } else if (method == "Exact") {
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

      invisible(self)
    },

    #' @description Executes the Expectation-Maximization (EM) algorithm based on the selected method. Additional parameters may be required depending on the method.
    #'
    #' @note This method can also be called as an S3 method with \code{predict()}.
    #'
    #' @param method \emph{(character)} The method for estimating the Expectation-Maximization (EM) algorithm. Options:
    #' \code{"Multinomial", "MVN CDF", "MVN PDF", "Hit and Run"} and \code{"Exact"} (default: \code{"Multinomial"}).
    #'
    #' @param probability_method \emph{(character)} The method for obtaining the initial probability. Options: \code{"Group proportional"},
    #' \code{"Proportional"} or \code{"Uniform"}. (default: \code{"Group proportional"}).
    #'
    #' @param iterations \emph{(integer(1))} The maximum amount of iterations to perform on the EM algorithm. (default: \code{1000}).
    #'
    #' @param maximum_time \emph{(integer(1))} The maximum time (in minutes) that the algorithm will iterate. (default: \code{1440}).
    #'
    #' @param stopping_threshold \emph{(numeric(1))} The minimum difference between consequent probabilities for stopping the iterations.
    #' (default: \code{0.001}).
    #'
    #' @param verbose \emph{(boolean(1))} Boolean indicating whether to print useful messages while iterating. (default: \code{FALSE}).
    #'
    #' @param ... Additional arguments required by specific methods:
    #'   \itemize{
    #'     \item \strong{"Hit and Run" Method:}
    #'       \itemize{
    #'         \item \code{step_size} \emph{(integer(1))}: The step size (M) for the Hit and Run algorithm. \cr
    #'         \item \code{samples} \emph{(integer(1))}: The number of samples (S) to generate.
    #'       }
    #'     \item \strong{"MVN CDF" Method:}
    #'       \itemize{
    #'         \item \code{multivariate_method} \emph{(character)}: The integration method. Default is \code{"Genz2"}. \cr
    #'         \item \code{multivariate_error} \emph{(numeric(1))}: The integration error threshold. Default is \code{1e-6}. \cr
    #'         \item \code{multivariate_iterations} \emph{(integer(1))}: The number of Monte Carlo iterations. Default is \code{5000}.
    #'       }
    #'   }
    #'
    #' @references Check \emph{Thraves, C. and Ubilla, P.: 'Fast Ecological Inference Algorithm for the R×C Case'} for benchmarking results.
    #'
    #' @return The modified eim object (for method chaining).
    #'
    #' @examples
    #' # Example 1: Compute the Expected Maximization with default values
    #'
    #' simulations <- simulate_elections(num_ballots = 20,
    #' 				num_candidates = 5,
    #' 				num_groups = 3,
    #' 				ballot_voters = rep(100, 20))
    #'
    #' model <- eim$new(simulations$X, simulations$W)
    #' model$compute() # Retrieves the object with updated attributes
    #'
    #' # Example 2: Compute the Expected Maximization for the Hit and Run method
    #' \donttest{
    #' 	model$compute(method = "Hit and Run",
    #' 					step_size = 3000,
    #' 					samples = 1000)
    #' 	}
    #'
    #' # Example 3: Omit arguments to the Hit and Run method
    #' \dontrun{
    #' 	model$compute(method = "Hit and Run",
    #' 					step_size = 3000)
    #' 	# Error; must hand in samples parameter too
    #' 	}
    #'
    #' # Example 4: Run the MVN CDF with default values
    #' model$compute(method = "MVN CDF")
    #'
    #' # Example 5: Run a Exact estimation with user defined parameters
    #' model$compute(method = "Exact",
    #' 				probability_method = "Uniform",
    #'  			iterations = 5,
    #' 				stopping_threshold = 1e-3)
    #' # Verbose was omitted
    compute = function(method = "Multinomial",
                       probability_method = "Group proportional",
                       iterations = 1000,
                       maximum_time = 1440,
                       stopping_threshold = 0.001,
                       verbose = FALSE, ...) {
      params <- list(...)

      # Check if the method provided is valid
      private$validate_methods(method, params)
      # Define the method
      self$method <- method

      # ========= HIT AND RUN ========= #
      if (self$method == "Hit and Run") {
        # Define the attributes
        private$hr_step_size <- params$step_size
        private$hr_samples <- params$samples

        # Run the EM algorithm for the Hit and Run method.
        resulting_values <- EMAlgorithmHitAndRun(
          probability_method,
          iterations,
          maximum_time,
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
          maximum_time,
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
          maximum_time,
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
        # ncol(self$X) is used just to get the total number of candidates
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
      self$finish_state <- resulting_values$stopping_reason
      private$been_computed <- TRUE

      invisible(self)
    },

    #' @description Runs a bootstrap algorithm to estimate the \strong{standard deviation} of probabilities
    #'
    #' @note This method can also be called as an S3 method with \code{std()}.
    #'
    #' Computes the EM algorithm an "\code{bootstrap_iterations}" times by sampling without replacement an "\code{ballot_boxes}" amount of tables. New samples are made
    #' per iteration. It'll obtain the standard deviation as a result from every iteration.
    #'
    #' @param bootstrap_iterations \emph{(integer(1))} Amount of EM computations.
    #'
    #' @param ballot_boxes \emph{(integer(1))} Amount of ballot boxes to use as a sample. It must be strictly smaller than the current amount of ballot boxes.
    #'
    #' @param method \emph{(character)} The method for estimating the Expectation-Maximization (EM) algorithm. Options:
    #' \code{"Multinomial", "MVN CDF", "MVN PDF", "Hit and Run"} and \code{"Exact"} (default: \code{"Multinomial"}).
    #'
    #' @param probability_method \emph{(character)} The method for obtaining the initial probability. Options: \code{"Group proportional"},
    #' \code{"Proportional"} or \code{"Uniform"}. (default: \code{"Group proportional"}).
    #'
    #' @param iterations \emph{(integer(1))} The maximum amount of iterations to perform on the EM algorithm. (default: \code{1000}).
    #'
    #' @param maximum_time \emph{(integer(1))} The maximum time (in minutes) that the EM algorithm will iterate. (default: \code{1440}).
    #'
    #' @param stopping_threshold \emph{(numeric(1))} The minimum difference between consequent probabilities for stopping the iterations.
    #' (default: \code{0.001}).
    #'
    #' @param verbose \emph{(boolean(1))} Boolean indicating whether to print useful messages while iterating. (default: \code{FALSE}).
    #'
    #' @param ... Additional arguments required by specific methods:
    #'   \itemize{
    #'     \item \strong{"Hit and Run" Method:}
    #'       \itemize{
    #'         \item \code{step_size} \emph{(integer(1))}: The step size (M) for the Hit and Run algorithm. \cr
    #'         \item \code{samples} \emph{(integer(1))}: The number of samples (S) to generate.
    #'       }
    #'     \item \strong{"MVN CDF" Method:}
    #'       \itemize{
    #'         \item \code{multivariate_method} \emph{(character)}: The integration method. Default is \code{"Genz2"}. \cr
    #'         \item \code{multivariate_error} \emph{(numeric(1))}: The integration error threshold. Default is \code{1e-6}. \cr
    #'         \item \code{multivariate_iterations} \emph{(integer(1))}: The number of Monte Carlo iterations. Default is \code{5000}.
    #'       }
    #'   }
    #'
    #' @return The modified eim object, under the field of \code{$standard_deviation}.
    #'
    #' @examples
    #'
    #' simulations <- simulate_elections(num_ballots = 20,
    #' 				num_candidates = 5,
    #' 				num_groups = 3,
    #' 				ballot_voters = rep(100, 20))
    #'
    #' model <- eim$new(X = simulations$X, W = simulations$W)
    #'
    #' model$bootstrap(30, 10)
    #'
    #' model$std # An estimate of the probabilities sd.
    bootstrap = function(bootstrap_iterations, ballot_boxes,
                         method = "Multinomial",
                         probability_method = "Group proportional",
                         iterations = 1000,
                         stopping_threshold = 0.001,
                         verbose = FALSE, ...) {
      params <- list(...)
      # By default we won't verbose the computing intermediate results, since they got missing values.
      compute_params <- c(
        list(
          method = method,
          probability_method = probability_method,
          iterations = iterations,
          stopping_threshold = stopping_threshold,
          verbose = FALSE
        ),
        params
      )

      # Parameter checking
      if (ballot_boxes >= nrow(self$X)) {
        stop("Bootstrap error: The sampling size must be smaller than the amount of ballot boxes")
      }
      if (method == "Exact" && bootstrap_iterations >= 5) {
        message("Warning: The bootstraping method may take a while.")
      }

      # Array for storing intermediate results
      results_array <- vector("list", length = bootstrap_iterations)

      # Main loop
      for (i in seq_len(bootstrap_iterations)) {
        if (verbose) {
          message("Running iteration:\t", i)
        }
        # Sample according the ballot boxes
        sample <- sample(1:nrow(self$X), size = ballot_boxes, replace = FALSE)
        iteration_X <- self$X[sample, , drop = FALSE]
        iteration_W <- self$W[sample, , drop = FALSE]
        # Create a temporary eim object and compute
        sample_object <- eim$new(iteration_X, iteration_W)
        do.call(sample_object$compute, compute_params)

        results_array[[i]] <- sample_object$probability
        rm(sample_object)
      }

      # Convert results into a 3D tensor
      final_array <- array(unlist(results_array), dim = c(ncol(self$W), ncol(self$X), bootstrap_iterations))
      sd_matrix <- apply(final_array, MARGIN = c(1, 2), FUN = sd)

      # Final results
      private$bootstrap_result <- sd_matrix
      private$bootstrap_called <- TRUE

      invisible(self)
    },

    #' @description According to the state of the algorithm (either computed or not), it prints a message with its most relevant parameters
    #'
    #' @note This method can also be called as an S3 method with \code{print()}.
    #'
    #' @return The own object for method chaining.
    #'
    #' @examples
    #'
    #' simulations <- simulate_elections(num_ballots = 20,
    #' 				num_candidates = 5,
    #' 				num_groups = 3,
    #' 				ballot_voters = rep(100, 20))
    #'
    #' model <- eim$new(simulations$X, simulations$W)
    #' print(model) # Will print the X and W matrix.
    #'
    #' model$compute()
    #' print(model) # Will print the X and W matrix among the EM results.
    print = function() {
      cat("eim ecological inference model\n")
      # Determine if truncation is needed
      truncated_X <- (nrow(self$X) > 5)
      truncated_W <- (nrow(self$W) > 5)

      cat("Candidate matrix (X) [b x c]:\n")
      print(self$X[1:min(5, nrow(self$X)), ])
      if (truncated_X) cat("...\n")

      cat("Group matrix (W) [b x g]:\n")
      print(self$W[1:min(5, nrow(self$W)), ])
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

    #' @description Shows, in form of a list, a selection of the most important attributes. It'll retrieve the method, number of candidates, ballots and groups and the principal results of the EM algorithm.
    #'
    #' @note This method can also be called as an S3 method with \code{summary()}.
    #'
    #' @return \emph{(list)} A list with the method, candidates, ballots, group, probabilities and log-likelihood.
    #'
    #' @examples
    #'
    #' simulations <- simulate_elections(num_ballots = 5,
    #' 				num_candidates = 3,
    #' 				num_groups = 2,
    #' 				ballot_voters = rep(100, 5))
    #'
    #' model <- eim$new(simulations$X, simulations$W)
    #' a_list <- model$summary()
    #' a_list$method # Not computed yet
    #' a_list$groups # 2
    #' a_list$ballots # 5
    #' names(a_list)
    #' # "candidates" "groups" "ballots" "method" "probabilities" "logLikelihood"
    summary = function() {
      list(
        candidates = ncol(self$X),
        groups = ncol(self$W),
        ballots = ncol(self$X),
        method = if (!is.null(self$method)) self$method else "Not computed yet",
        probabilities = if (!is.null(self$probability)) self$probability else "Not computed yet",
        logLikelihood = if (!is.null(self$logLikelihood)) self$logLikelihood else "Not computed yet",
        std = if (!is.null(self$std)) self$std else "Not computed yet"
      )
    },

    #' @description Saves the current eim object to a specified file. The results can be saved in:
    #'   \itemize{
    #'     \item \strong{RDS (Binary format):} Preserves object structure for future use in R. \cr
    #'     \item \strong{JSON}: Saves model data in a human-readable format. \cr
    #'     \item \strong{CSV}: Saves probability matrix in a tabular format.
    #'   }
    #'
    #' @param filename \emph{(character)} The file name where the results should be saved with its extension.
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
    #'
    #' simulations <- simulate_elections(num_ballots = 20,
    #' 				num_candidates = 5,
    #' 				num_groups = 3,
    #' 				ballot_voters = rep(100, 20))
    #'
    #' model <- eim$new(simulations$X, simulations$W)
    #' model$compute()
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

        # Add the standard deviation only if it exists
        if (private$bootstrap_called) {
          json_data$std <- self$std
        }

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
    #' @field samples \emph{(integer(1))} Active variable to show the \code{Hit and Run} samples if and only if the method is \code{"Hit and Run"}
    samples = function() {
      if (self$method == "Hit and Run") {
        return(private$hr_samples)
      } else {
        return(NULL)
      }
    },

    #' @field step_size \emph{(integer(1))} Active variable to show the \code{Hit and Run} step size if and only if the method is \code{"Hit and Run"}
    step_size = function() {
      if (self$method == "Hit and Run") {
        return(private$hr_step_size)
      } else {
        return(NULL)
      }
    },

    #' @field multivariate_method \emph{(character)} Active variable to show the method used to estimate the Multivariate
    #' Normal CDF, if and only if \code{self$method} is \code{"MVN CDF"}
    multivariate_method = function() {
      if (self$method == "MVN CDF") {
        return(private$mvn_method)
      } else {
        return(NULL)
      }
    },

    #' @field multivariate_error \emph{(numeric(1))} Active variable to show the error threshold for the Montecarlo simulation of the Multivariate Normal CDF
    multivariate_error = function() {
      if (self$method == "MVN CDF") {
        return(private$mvn_error)
      } else {
        return(NULL)
      }
    },

    #' @field multivariate_iterations \emph{(numeric(1))} Active variable to show the number of iterations for the Montecarlo simulation of the Multivariate Normal CDF
    multivariate_iterations = function() {
      if (self$method == "MVN CDF") {
        return(private$mvn_iterations)
      } else {
        return(NULL)
      }
    },
    #' @field std \emph{(matrix)} Active variable to show an estimate of the probability standard deviation if the bootstrapping method has been called
    std = function() {
      if (private$bootstrap_called) {
        return(private$bootstrap_result)
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
    # PRIVATE_FIELD bootstrap_called Used to see if the bootstrapping method has been called.
    bootstrap_called = FALSE,
    # PRIVATE_FIELD boostrap_result The result of the bootstrapping method.
    bootstrap_result = 0,
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
        if (is.null(params$step_size) ||
          params$step_size < 0 ||
          abs(params$step_size %% 1) > tolerance) { # Last condition: check if it has decimal part
          stop("Hit and Run:\tA valid 'step_size' wasn't provided")
        }

        # Check for samples
        if (is.null(params$samples) ||
          params$samples < 0 ||
          abs(params$step_size %% 1) > tolerance) { # Last condition: check if it has decimal part
          stop("Hit and Run:\tA valid 'samples' wasn't provided")
        }

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
        if (!is.null(params$multivariate_method) && !(params$multivariate_method %in% valid_cdf_methods)) {
          stop(
            "MVN CDF\t: The method for estimating the CDF isn't supported. Must be one of: ",
            paste(valid_cdf_methods, collapse = ", ")
          )
        }

        # If there's a multivariate error, it has to be valid.
        if (!is.null(params$multivariate_error) && params$multivariate_error < 0) {
          stop("MVN CDF\t: The error threshold for the Montecarlo simulation isn't valid")
        }

        # Check if there's a multivariate iterations, otherwise, use 5000 as default
        tolerance <- .Machine$double.eps^0.5
        if (!is.null(params$multivariate_iterations) &&
          (params$mvn_iterations < 0 || abs(params$mvn_iterations %% 1) > tolerance)) {
          stop("MVN CDF\t: An invalid iteration parameter was handed")
        }
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

#' @rdname eim
#' @export
summary.eim <- function(object) {
  object$summary()
}

#' @rdname eim
#' @export
predict.eim <- function(object, ...) {
  params <- list(...)
  do.call(object$compute, params) # Calls compute() with the right arguments
  return(object$probability)
}

#' @rdname eim
#' @export
as.matrix.eim <- function(object) {
  if (is.null(object$probability)) {
    stop(paste0("Probability matrix not available. Run compute() or ", object, "$compute() first."))
  }
  return(object$probability)
}

#' @rdname eim
#' @export
std.eim <- function(object, ...) {
  params <- list(...)
  do.call(object$bootstrap, params)
  return(object$std)
}

#' Update an existing eim model with a new EM algorithm computation
#'
#' This function updates an object with a new Expected Maximization computation
#' with other parameters.
#'
#' @param \emph{object} An eim object.
#'
#' @param ... New parameters to passed to \code{compute()}.
#'
#' @return Changes made on the updated object.
#'
#' @examples
#' simulations <- simulate_elections(
#'   num_ballots = 20,
#'   num_candidates = 5,
#'   num_groups = 3,
#'   ballot_voters = rep(100, 20)
#' )
#'
#' model <- eim$new(simulations$X, simulations$W)
#'
#' predict(model)
#'
#' update(model, "MVN PDF")
#'
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
#' of running \code{sum(object$W)} or \code{sum(object$X)}.
#'
#' @param object An eim object.
#'
#' @return The amount of voters.
#'
#' @examples
#' simulations <- simulate_elections(
#'   num_ballots = 20,
#'   num_candidates = 5,
#'   num_groups = 3,
#'   ballot_voters = rep(100, 20)
#' )
#'
#' model <- eim$new(simulations$X, simulations$W)
#' sum(model) # 2000
#'
#' @export
sum.eim <- function(object) {
  # Ensure the model has been computed before updating
  if (is.null(object$W)) stop("The object must be initialized.")

  return(sum(object$W))
}
