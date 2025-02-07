#' @title Estimate RxC Matrix with EM Algorithm
#' @name rxc
#' @description This R6 class implements an EM algorithm for estimating an RxC matrix from an ecological inference problem. It supports multiple methods for approximating the E-step including "Multinomial", "Hit and Run", "MVN CDF", "MVN PDF", and "Exact". The class can be initialized with candidate and group matrices or by reading a JSON file, and it provides methods for precomputation, running the algorithm, and saving results.
#' @export
rxc <- R6::R6Class("rxc",
                   public = list(
                       ## Public Members
                       #' @field X A (b x c) matrix with the observed candidate votes. Provided manually or loaded from JSON.
                       X = NULL,
                       #' @field W A (b x g) matrix with the observed group votes. Provided manually or loaded from JSON.
                       W = NULL,
                       #' @field method A character string indicating the selected EM method. Options: "Multinomial", "Hit and Run", "MVN CDF", "MVN PDF", "Exact".
                       method = NULL,
                       #' @field probability A (g x c) matrix storing the final estimated probabilities.
                       probability = NULL,
                       #' @field logLikelihood A numeric vector storing the log-likelihood values during the EM iterations.
                       logLikelihood = NULL,
                       #' @field total_iterations An integer representing the total number of iterations performed by the EM algorithm.
                       total_iterations = NULL,
                       #' @field total_time A numeric value representing the total time taken by the EM algorithm.
                       total_time = NULL,
                       
                       ## Class Constructor:
                       #' @section Class Constructor:
                       #' \describe{
                       #'   \item{\code{new(X = NULL, W = NULL, json_path = NULL)}}{
                       #'     \itemize{
                       #'       \item \strong{X}: A (b x c) matrix of candidate votes. Optional if \code{json_path} is provided.
                       #'       \item \strong{W}: A (b x g) matrix of group votes. Optional if \code{json_path} is provided.
                       #'       \item \strong{json_path}: A character string specifying the path to a JSON file containing matrices \code{X} and \code{W}.
                       #'     }
                       #'   }
                       #' }
                       #' @export
                       initialize = function(X = NULL, W = NULL, json_path = NULL) {
                           if (!is.null(json_path) && nzchar(json_path)) {
                               data <- jsonlite::fromJSON(json_path)
                               if (!all(c("X", "W") %in% names(data))) stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")
                               self$X <- t(as.matrix(data$X)) # ERASE THE TRANSPOSE LATER
                               self$W <- as.matrix(data$W)
                           } else {
                               if (is.null(X) || is.null(W)) stop("Either provide X and W or a non-empty JSON path")
                               self$X <- as.matrix(X)
                               self$W <- as.matrix(W)
                           }
                           # ---- C initialization ---- #
                           param_X <- matrix(as.numeric(self$X), nrow(self$X), ncol(self$X))
                           param_W <- matrix(as.numeric(self$W), nrow(self$W), ncol(self$W))
                           private$empty_candidates <- which(rowSums(self$X == 0) == ncol(self$X))
                           if (length(private$empty_candidates) != 0) param_X <- param_X[-private$empty_candidates, drop = FALSE]
                           RsetParameters(t(param_X), param_W)
                       },
                       
                       ## Public Methods:
                       #' @section Public Methods:
                       #' \describe{
                       #'   \item{\code{precompute(method, ...)}}{
                       #'     \itemize{
                       #'       \item \strong{Args:}
                       #'         \itemize{
                       #'           \item \code{method}: A character string specifying the method to precompute. Options: "Exact" or "Hit and Run".
                       #'           \item \code{...}: Additional arguments. For "Hit and Run", supply \code{step_size} (a positive integer) and \code{samples} (an integer).
                       #'         }
                       #'       \item \strong{Returns}: The modified \code{rxc} object (invisibly) for method chaining.
                       #'     }
                       #'   }
                       #'   
                       #'   \item{\code{compute(main_method, probability_method, iterations, stopping_threshold, verbose, ...)}}{
                       #'     \itemize{
                       #'       \item \strong{Args:}
                       #'         \itemize{
                       #'           \item \code{main_method}: A character string indicating the EM method. Options: "Multinomial", "MVN CDF", "MVN PDF", "Hit and Run", "Exact".
                       #'           \item \code{probability_method}: A character string for obtaining the initial probability. Options: "Group proportional", "Proportional", "Uniform".
                       #'           \item \code{iterations}: An integer specifying the maximum number of iterations (default: 1000).
                       #'           \item \code{stopping_threshold}: A numeric value for the stopping criterion (default: 0.001).
                       #'           \item \code{verbose}: A boolean indicating whether to print messages during computation (default: FALSE).
                       #'           \item \code{...}: Additional arguments required by specific methods. For example, for "Hit and Run", supply \code{step_size} and \code{samples}; for "MVN CDF", supply \code{multivariate_method}, \code{multivariate_error}, and \code{multivariate_iterations}.
                       #'         }
                       #'       \item \strong{Returns}: The modified \code{rxc} object (invisibly) for method chaining.
                       #'     }
                       #'   }
                       #'   
                       #'   \item{\code{print()}}{
                       #'     \itemize{
                       #'       \item \strong{Description}: Prints the model's principal attributes, including the candidate and group matrices and, if computed, the EM results.
                       #'       \item \strong{Returns}: The \code{rxc} object (invisibly) for method chaining.
                       #'     }
                       #'   }
                       #'   
                       #'   \item{\code{summary()}}{
                       #'     \itemize{
                       #'       \item \strong{Description}: Returns a list with key attributes such as the number of candidates, groups, ballots, selected method, a preview of the probabilities, and log-likelihood values.
                       #'       \item \strong{Returns}: A list with the principal attributes.
                       #'     }
                       #'   }
                       #'   
                       #'   \item{\code{save_results(filename)}}{
                       #'     \itemize{
                       #'       \item \strong{Args:}
                       #'         \itemize{
                       #'           \item \code{filename}: A character string specifying the file name (with extension) to save the results. Supported extensions are ".rds", ".json", or ".csv".
                       #'         }
                       #'       \item \strong{Returns}: The modified \code{rxc} object (invisibly) for method chaining.
                       #'     }
                       #'   }
                       #' }
                       #' @export
                       precompute = function(method, ...) {
                           params <- list(...)
                           if (method == "Hit and Run") {
                               if (is.null(params$step_size) || params$step_size < 0) stop("precompute(): A valid 'step_size' wasn't provided")
                               if (is.null(params$samples) || params$samples < 0) stop("precompute(): A valid 'samples' wasn't provided")
                               if (private$been_precomputed_hr &&
                                   (private$hr_step_size != params$step_size ||
                                    private$hr_samples != params$hr_samples)) {
                                   clean_hr_precompute()
                               }
                               private$hr_step_size <- params$step_size
                               private$hr_samples <- params$samples
                               message("Precomputing the Hit and Run method")
                               RprecomputeHR(as.integer(private$hr_samples), as.integer(private$hr_step_size))
                               private$been_precomputed_hr <- TRUE
                           } else if (method == "Exact") {
                               message("Precomputing the Exact method")
                               RprecomputeExact()
                               private$been_precomputed_exact <- TRUE
                           } else {
                               stop("Invalid method for precomputing. Must be either Hit and Run or Exact")
                           }
                           invisible(self)
                       },
                       
                       compute = function(main_method = "Multinomial",
                                          probability_method = "Group proportional",
                                          iterations = 1000,
                                          stopping_threshold = 0.001,
                                          verbose = FALSE, ...) {
                           params <- list(...)
                           valid_methods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
                           if (!is.character(main_method) || length(main_method) != 1 || !(main_method %in% valid_methods)) {
                               stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
                           }
                           self$method <- main_method
                           if (self$method %in% c("Exact", "Multinomial", "MVN PDF")) {
                               resulting_values <- EMAlgorithmAll(
                                   self$method,
                                   probability_method,
                                   iterations,
                                   stopping_threshold,
                                   verbose
                               )
                               if (self$method == "Exact") private$been_precomputed_exact <- TRUE
                           } else if (self$method == "Hit and Run") {
                               if (is.null(params$step_size) || params$step_size < 0) stop("compute(): A valid 'step_size' wasn't provided")
                               if (is.null(params$samples) || params$samples < 0) stop("compute(): A valid 'samples' wasn't provided")
                               if (private$been_precomputed_hr &&
                                   (private$hr_step_size != params$step_size ||
                                    private$hr_samples != params$hr_samples)) {
                                   clean_hr_precompute()
                               }
                               private$hr_step_size <- params$step_size
                               private$hr_samples <- params$samples
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
                               if (is.null(params$multivariate_method) || !params$multivariate_method %in% c("Genz", "Genz2")) {
                                   private$mvn_method <- "Genz2"
                               } else {
                                   private$mvn_method <- params$multivariate_method
                               }
                               if (is.null(params$multivariate_error) || params$multivariate_error < 0) {
                                   private$mvn_error <- 1e-5
                               } else {
                                   private$mvn_error <- as.numeric(params$multivariate_error)
                               }
                               if (is.null(params$multivariate_iterations) || params$mvn_iterations < 0) {
                                   private$mvn_iterations <- as.integer(5000)
                               } else {
                                   private$mvn_iterations <- as.integer(params$mvn_iterations)
                               }
                               if (is.null(params$multivariate_error) || params$multivariate_error < 0) private$mvn_error <- 0.000001
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
                           if (length(private$empty_candidates) != 0) {
                               new_probability <- matrix(0, nrow = nrow(resulting_values$result), ncol = ncol(self$X))
                               non_zero_cols <- setdiff(seq_len(ncol(self$X)), private$empty_candidates)
                               new_probability[, non_zero_cols] <- resulting_values$result
                               self$probability <- new_probability
                           } else {
                               self$probability <- resulting_values$result
                           }
                           self$logLikelihood <- resulting_values$log_likelikelihood
                           self$total_iterations <- resulting_values$total_iterations
                           self$total_time <- resulting_values$total_time
                           private$been_computed <- TRUE
                           invisible(self)
                       },
                       
                       print = function() {
                           cat("rxc ecological inference model\n")
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
                               cat("Total Iterations:\t", self$total_iterations, "\n")
                               cat("Total Time (s):\t", self$total_time, "\n")
                               cat("Estimated probability:\n")
                               print(self$probability)
                           }
                           invisible(self)
                       },
                       
                       summary = function() {
                           list(
                               candidates = nrow(self$X),
                               groups = ncol(self$W),
                               ballots = ncol(self$X),
                               method = if (!is.null(self$method)) self$method else "Not computed yet",
                               probabilities = if (!is.null(self$probability)) head(self$probability, 5) else "Not computed yet",
                               logLikelihood = if (!is.null(self$logLikelihood)) tail(self$logLikelihood, 5) else "Not computed yet"
                           )
                       },
                       
                       save_results = function(filename) {
                           if (!is.character(filename) || length(filename) != 1) {
                               stop("Invalid filename. Please provide a valid file path as a character string.")
                           }
                           file_ext <- tools::file_ext(filename)
                           if (file_ext == "rds") {
                               saveRDS(self, file = filename)
                               message("Results saved as RDS: ", filename)
                           } else if (file_ext == "json") {
                               json_data <- list(
                                   method = self$method,
                                   probability = as.list(as.data.frame(self$finalProbability)),
                                   logLikelihood = self$logLikelihood,
                                   totalIterations = self$total_iterations,
                                   totalTime = self$total_time
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
                           invisible(self)
                       }
                   ),
                   active = list(
                       #' @field samples Active binding showing the Hit and Run samples if the selected method is "Hit and Run".
                       samples = function() {
                           if (self$method == "Hit and Run") {
                               return(private$hr_samples)
                           } else {
                               return(NULL)
                           }
                       },
                       #' @field step_size Active binding showing the Hit and Run step size if the selected method is "Hit and Run".
                       step_size = function() {
                           if (self$method == "Hit and Run") {
                               return(private$hr_step_size)
                           } else {
                               return(NULL)
                           }
                       },
                       #' @field multivariate_method Active binding for the MVN CDF integration method if the selected method is "MVN CDF".
                       multivariate_method = function() {
                           if (self$method == "MVN CDF") {
                               return(private$mvn_method)
                           } else {
                               return(NULL)
                           }
                       },
                       #' @field multivariate_error Active binding for the MVN CDF error threshold if the selected method is "MVN CDF".
                       multivariate_error = function() {
                           if (self$method == "MVN CDF") {
                               return(private$mvn_error)
                           } else {
                               return(NULL)
                           }
                       },
                       #' @field multivariate_iterations Active binding for the MVN CDF number of Monte Carlo iterations if the selected method is "MVN CDF".
                       multivariate_iterations = function() {
                           if (self$method == "MVN CDF") {
                               return(private$mvn_iterations)
                           } else {
                               return(NULL)
                           }
                       }
                   ),
                   private = list(
                       #' @keywords internal been_called Boolean that determines if the object has been called before (used for C cleanup).
                       been_called = FALSE,
                       #' @keywords internal been_precomputed_hr Boolean that indicates if precomputation for the Hit and Run method has been performed.
                       been_precomputed_hr = FALSE,
                       #' @keywords internal been_precomputed_exact Boolean that indicates if precomputation for the Exact method has been performed.
                       been_precomputed_exact = FALSE,
                       #' @keywords internal empty_candidates Integer vector with indices of candidates that received no votes.
                       empty_candidates = integer(0),
                       #' @keywords internal been_computed Boolean that indicates if the EM algorithm has been computed.
                       been_computed = FALSE,
                       #' @keywords internal hr_samples The number of samples used in the Hit and Run method.
                       hr_samples = 0,
                       #' @keywords internal hr_step_size The step size used in the Hit and Run method.
                       hr_step_size = 0,
                       #' @keywords internal mvn_method The integration method used for the MVN CDF.
                       mvn_method = "Genz2",
                       #' @keywords internal mvn_error The error threshold used for the MVN CDF Monte Carlo simulation.
                       mvn_error = 0.000001,
                       #' @keywords internal mvn_iterations The number of iterations for the MVN CDF Monte Carlo simulation.
                       mvn_iterations = 5000,
                       #' @keywords internal finalize Destructor to clean allocated memory via C routines.
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

## S3 Methods

#' @title Summarize rxc Object
#' @name summary.rxc
#' @description Returns a list with the principal attributes of an rxc object.
#' @param object An rxc object.
#' @return A list with the principal attributes.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:10, 2, 5))
#' summary(model)
#' @export
summary.rxc <- function(object) {
    object$summary()
}

#' @title Predict RxC Probabilities
#' @name predict.rxc
#' @description Predicts the RxC probability using the EM algorithm of an rxc object. This function is a wrapper around \code{compute()}.
#' @param object An rxc object.
#' @param ... Additional arguments to pass to \code{compute()}.
#' @return A matrix of estimated probabilities.
#' @examples
#' model <- rxc$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#' predict(model, "Hit and Run", step_size = 1000, samples = 5000)
#' @export
predict.rxc <- function(object, ...) {
    params <- list(...)
    do.call(object$compute, params)
    return(object$probability)
}

#' @title Extract Probability Matrix from rxc Object
#' @name as.matrix.rxc
#' @description Extracts the probability matrix from an rxc object.
#' @param object An rxc object.
#' @return A matrix containing the estimated probabilities.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' as.matrix(model)
#' @export
as.matrix.rxc <- function(object) {
    if (is.null(object$probability)) {
        stop("Probability matrix not available. Run compute() or object$compute() first.")
    }
    return(object$probability)
}

#' @title Update an rxc Object
#' @name update.rxc
#' @description Updates an existing rxc model with a new EM algorithm computation.
#' @param object An rxc object.
#' @param ... New parameters to pass to \code{compute()}.
#' @return The updated rxc object.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' update(model, "MVN PDF")
#' @export
update.rxc <- function(object, ...) {
    params <- list(...)
    if (is.null(object$probability)) {
        stop("Model must be computed before updating. Run compute() or object$compute() first.")
    }
    do.call(object$compute, params)
    return(object)
}

#' @title Get Number of Candidates
#' @name candidates.rxc
#' @description Returns the current number of candidates in an rxc object.
#' @param object An rxc object.
#' @return An integer representing the number of candidates.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' candidates(model)
#' @export
candidates.rxc <- function(object) {
    if (is.null(object$X)) stop("The object must be initialized.")
    return(ncol(object$X))
}

#' @title Get Number of Groups
#' @name groups.rxc
#' @description Returns the current number of groups in an rxc object.
#' @param object An rxc object.
#' @return An integer representing the number of groups.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' groups(model)
#' @export
groups.rxc <- function(object) {
    if (is.null(object$W)) stop("The object must be initialized.")
    return(ncol(object$W))
}

#' @title Get Number of Ballot Boxes
#' @name ballots.rxc
#' @description Returns the number of ballot boxes in an rxc object.
#' @param object An rxc object.
#' @return An integer representing the number of ballots.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' ballots(model)
#' @export
ballots.rxc <- function(object) {
    if (is.null(object$W)) stop("The object must be initialized.")
    return(nrow(object$W))
}

#' @title Get Total Number of Voters
#' @name sum.rxc
#' @description Returns the total number of voters in the system from an rxc object.
#' @param object An rxc object.
#' @return An integer representing the total number of voters.
#' @examples
#' model <- rxc$new(X = matrix(1:15, 5, 3), W = matrix(1:50, 5, 10))
#' sum(model)
#' @export
sum.rxc <- function(object) {
    if (is.null(object$W)) stop("The object must be initialized.")
    return(sum(object$W))
}
