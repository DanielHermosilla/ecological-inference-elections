library(jsonlite)

#' S3 Object for the Expectation-Maximization Algorithm
#'
#' This constructor creates an `"eim"` S3 object, either by using matrices
#' `X` and `W` directly or by reading them from a JSON file. Each
#' `eim` object encapsulates the data (votes for candidates and demographic
#' groups) required by the underlying Expectation-Maximization algorithm.
#'
#' @param X A `(b x c)` matrix of candidate votes per ballot box.
#'   This is optional if a valid `json_path` is provided.
#' @param W A `(b x g)` matrix of group votes per ballot box.
#'   This is optional if a valid `json_path` is provided.
#' @param json_path A character string specifying the path to a JSON file that
#'   contains the `"X"` and `"W"` matrices. If both `X` and `W`
#'   are not supplied, they will be read from this file.
#'
#' @details
#' If `X` and `W` are directly supplied, they must match the correct
#' dimensions of ballot boxes (b). Alternatively, if `json_path` is provided, the function expects
#' the JSON file to contain elements named `"X"` and `"W"` under the
#' top-level object.
#'
#' Internally, this function also initializes the corresponding instance within
#' the low-level (C-based) API, ensuring the data is correctly registered for
#' further processing by the EM algorithm.
#'
#' @return A list of class `"eim"` containing:
#' \describe{
#'   \item{\code{X}}{The candidate votes matrix \eqn{(b \times c)}.}
#'   \item{\code{W}}{The group votes matrix \eqn{(b \times g)}.}
#' }
#'
#' @note
#' For an alternative way to generate `X` and `W`, you can use the
#' [simulate_elections] function, which provides several parameters
#' (such as `group_proportions`, `candidate_proportions`, and
#' `probability`) for controlling the underlying distributions of votes.
#' Afterwards, you may pass the resulting matrices to `eim()`.
#'
#' @section Methods:
#' In addition to this constructor, the `"eim"` class provides several
#' S3 methods for common operations. Some of these methods are fully documented,
#' while others are ommited due to its straightfoward implementantion. The available methods are:
#' - [compute] -- Runs the Expected-Maximization algorithm.
#' - [sd] -- Estimates the standard deviation of probabilities
#' - [save] -- Save the results to a specified file.
#' - `var` -- Estimates the variance of probabilities.
#' - `print` -- Print useful information about the object.
#' - `summary` -- Shows, in form of a list, the most important attributes.
#' - `as.matrix` -- Returns the probability matrix.
#' - `write.csv` -- Writes the probability matrix in a `.csv` file.
#' - `dput` -- Writes the object in a `.rda` file.
#' - `logLik` -- Returns the log-likelihood from the last iteration.
#'
#' @examples
#' # Example 1: Create an eim object from existing matrices
#' X_mat <- matrix(1:9, nrow = 3, ncol = 3)
#' W_mat <- matrix(1:9, nrow = 3, ncol = 3)
#' model <- eim(X = X_mat, W = W_mat)
#'
#' # Example 2: Create an eim object from a JSON file
#' \dontrun{
#' model2 <- eim(json_path = "path/to/file.json")
#' }
#'
#' # Example 3: Use simulate_elections with optional parameters, then create an eim object
#' \dontrun{
#' # Simulate data for 2 ballot boxes, 3 candidates, 2 groups,
#' # and custom group/candidate proportions
#' sim_result <- simulate_elections(
#'     num_ballots = 2,
#'     num_candidates = 3,
#'     num_groups = 2,
#'     group_proportions = c(0.4, 0.6),
#'     candidate_proportions = c(0.2, 0.3, 0.5)
#' )
#'
#' # Construct an eim object with the simulated data
#' model3 <- eim(X = sim_result$X, W = sim_result$W)
#' }
#'
#' @export
eim <- function(X = NULL, W = NULL, json_path = NULL) {
    # Load data from JSON if a path is provided
    if (!is.null(json_path) && nzchar(json_path)) {
        matrices <- .validate_json_eim(json_path) # nolint
        X <- as.matrix(matrices$X)
        W <- as.matrix(matrices$W)
    }
    # Perform matricial validation
    .validate_eim(X, W) # nolint

    # Convert matrices for C API
    param_X <- matrix(as.numeric(X), nrow(X), ncol(X))
    param_X <- t(param_X)
    param_W <- matrix(as.numeric(W), nrow(W), ncol(W))
    RsetParameters(param_X, param_W)

    # Create the S3 object
    obj <- list(
        X = X,
        W = W
    )

    class(obj) <- "eim"
    obj
}
#' @title Compute the Expected-Maximization Algorithm
#'
#' @description
#' Executes the Expectation-Maximization (EM) algorithm using a specified method for estimating
#' the E-step. Certain methods may require additional arguments, which can be passed through `...`.
#'
#' @param object An object of class `eim`, which can be created using the [eim] function.
#'   This parameter shouldn't be called if you prefer to supply the `X` and `W` matrices or a `json_path`.
#'
#' @param X A `(b x c)` matrix representing candidate votes per ballot box. This is required if neither an
#'   `eim` object or a JSON path is provided.
#'
#' @param W A `(b x g)` matrix representing group votes per ballot box. This is required if neither an
#'   `eim` object or a JSON path is provided.
#'
#' @param json_path A path to a JSON file containing `X` and `W`. This is required if neither an
#'   `eim` object or matrices are provided.
#'
#' @param method An optional string specifying the method used for estimating the E-step. Valid
#'   options are `mult`, `mvn_cdf`, `mvn_pdf`, `hnr`, or `exact`.
#'   The default is `mult`.
#'
#' @param initial_prob An optional string specifying the method used to obtain the initial
#'   probability. Accepted values are `group_proportional`, `proportional`, or `uniform`.
#'   The default is `group_proportional`.
#'
#' @param maxiter An optional integer indicating the maximum number of Expected-Maximization algorithm iterations.
#'   The default value is `1000`.
#'
#' @param maxtime An optional integer specifying the maximum running time (in minutes) for the
#'   algorithm. The default value is `1440`, which corresponds to one day.
#'
#' @param stop_threshold An optional numeric value indicating the minimum difference between
#'   consecutive probability estimates required to stop iterating. The default value is `0.001`.
#'
#' @param verbose An optional boolean indicating whether to print informational messages during
#'   iteration. The default value is `FALSE`.
#'
#' @param step_size An optional integer specifying the step size (`M`) for the **Hit and Run**
#'   algorithm. This parameter is only applicable when using the **Hit and Run** method and will
#'   be ignored otherwise. The default value is `3000`.
#'
#' @param samples An optional integer indicating the number of samples (`S`) to generate for the
#'   **Hit and Run** method. This parameter is only relevant when the **Hit and Run** method is selected.
#'   The default value is `1000`.
#'
#' @param monte_method An optional string specifying the method used to estimate the **Multivariate Normal CDF (MVN CDF)**
#'   via a Monte Carlo simulation. Accepted values are "`genz`" and "`genz2`", with "`genz2`"
#'   set as the default.
#'
#' @param monte_error An optional numeric value defining the error threshold for the Monte Carlo
#'   simulation when estimating the **MVN CDF**. The default value is `1e-6`.
#'
#' @param monte_sample An optional integer specifying the number of Monte Carlo
#'   samples for the **MVN CDF** method. The default value is `5000`.
#'
#' @references
#' [Thraves, C. and Ubilla, P.: *"Fast Ecological Inference Algorithm for the R×C Case"*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834). Aditionally, the MVN CDF is computed by the methods introduced in [Genz, A. (2000). Numerical computation of multivariate normal probabilities. *Journal of Computational and Graphical Statistics*](https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities)
#'
#' @note
#' This function can be executed using one of three mutually exclusive approaches:
#'
#' 1. By providing an existing `eim` object.
#' 2. By supplying both input matrices (`X` and `W`) directly.
#' 3. By specifying a JSON file (`json_path`) containing the matrices.
#'
#' These input methods are **mutually exclusive**, meaning that you must provide **exactly one** of
#' these options. Attempting to provide more than one or none of these inputs will result in an error.
#'
#' When called with an `eim` object, the function updates the object with the computed results.
#' If an `eim` object is not provided, the function will create one internally using either the
#' supplied matrices or the data from the JSON file before executing the algorithm.
#'
#' @seealso The [eim] object implementation.
#'
#' @return
#' The function updates (or creates) the input object by adding the following attributes:
#' \describe{
#'   \item{probability}{The estimated probability matrix.}
#'   \item{logLik}{An array containing the log-likelihood values from each iteration.}
#'   \item{iterations}{The total number of iterations performed by the EM algorithm.}
#'   \item{time}{The total execution time of the algorithm.}
#'   \item{status}{The final status ID of the algorithm upon completion.}
#'   \item{message}{The finishing status displayed as a message.}
#' }
#' Aditionally, it will update the object with additional parameters in case the method provided is `hnr` or `mvn_cdf`.
#'
#' @usage
#' compute(
#' 		object = NULL,
#'      X = NULL,
#'      W = NULL,
#'      json_path = NULL,
#'      method = "Multinomial",
#'      initial_prob = "Group proportional",
#'      maxiter = 1000,
#'      maxtime = 1440,
#'      stop_threshold = 0.001,
#'      verbose = FALSE,
#'      ...
#' )
#'
#' @examples
#' # Example 1: Compute the Expected-Maximization with default settings
#' simulations <- simulate_elections(
#'     num_ballots = 20,
#'     num_candidates = 5,
#'     num_groups = 3,
#'     ballot_voters = rep(100, 20)
#' )
#' model <- eim(simulations$X, simulations$W)
#' compute(model) # Retrieves the object with updated attributes
#'
#' # Example 2: Compute the Expected-Maximization using the Hit and Run method
#' \donttest{
#' compute(
#'     object = model,
#'     method = "hnr",
#'     step_size = 3000,
#'     samples = 1000
#' )
#' }
#'
#' # Example 3: Run the MVN CDF method with default settings
#' compute(object = model, method = "mvn_cdf")
#'
#' # Example 4: Perform an Exact estimation using user-defined parameters
#' \dontrun{
#' compute(
#'     json_path = "a/json/file.json",
#'     method = "exact",
#'     initial_prob = "uniform",
#'     maxiter = 5,
#'     stop_threshold = 1e-3
#' )
#' }
#'
#' @name compute
#' @export
compute.eim <- function(object = NULL,
                        X = NULL,
                        W = NULL,
                        json_path = NULL,
                        method = "mult",
                        initial_prob = "group proportional",
                        maxiter = 1000,
                        maxtime = 1440,
                        stop_threshold = 0.001,
                        verbose = FALSE, ...) {
    params <- list(...)

    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }
    # Check if the method provided is valid
    valid_methods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")

    if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
        stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    object$method <- method

    if (method == "Hit and Run") {
        # ========= HIT AND RUN ========= #
        # Validate step size and samples
        tolerance <- .Machine$double.eps^0.5
        # Check for the step size
        if (is.null(params$step_size) ||
            params$step_size < 0 ||
            abs(params$step_size %% 1) > tolerance) { # Last condition: check if it has decimal part
            stop("Hit and Run:\tA valid 'step_size' wasn't provided")
        } else {
            # Use 3000 as a default step size
            params$step_size <- 3000
        }

        # Check for samples
        if (is.null(params$samples) ||
            params$samples < 0 ||
            abs(params$samples %% 1) > tolerance) { # Last condition: check if it has decimal part
            stop("Hit and Run:\tA valid 'samples' wasn't provided")
        } else {
            # Use 1000 samples as default
            params$samples <- 1000
        }

        # Run the EM algorithm for the Hit and Run method.
        resulting_values <- EMAlgorithmHitAndRun(
            initial_prob,
            iterations,
            time,
            stop_threshold,
            verbose,
            as.integer(params$step_size),
            as.integer(params$samples)
        )

        object$samples <- params$samples
        object$step_size <- params$step_size
        # Deallocates memory from C
        clean_hr_precompute()
    } else if (method == "MVN CDF") {
        # ========= MVN CDF ========= #
        # Check if there's a montecarlo method, otherwise, use 'Genz2' as default
        valid_cdf_methods <- c("Genz, Genz2")
        # If there's a parameter for MVN CDF; it has to be valid (could be omitted)
        if (!is.null(params$montecarlo_method) && !(params$montecarlo_method %in% valid_cdf_methods)) {
            stop(
                "MVN CDF\t: The method for estimating the CDF isn't supported. Must be one of: ",
                paste(valid_cdf_methods, collapse = ", ")
            )
        } else {
            # Use Genz2 as default
            params$montecarlo_method <- "Genz2"
        }

        # If there's a montecarlo error, it has to be valid.
        if (!is.null(params$montecarlo_error) && params$montecarlo_error < 0) {
            stop("MVN CDF\t: The error threshold for the Montecarlo simulation isn't valid")
        } else {
            # Use 1e-6 as default montecarlo error
            params$montecarlo_error <- 1e-6
        }

        # Check if there's montecarlo iterations, otherwise, use 5000 as default
        tolerance <- .Machine$double.eps^0.5
        if (!is.null(params$montecarlo_iterations) &&
            (params$montecarlo_iterations < 0 || abs(params$montecarlo_iterations %% 1) > tolerance)) {
            stop("MVN CDF\t: An invalid iteration parameter was handed")
        } else {
            params$montecarlo_iterations <- 5000
        }

        # Run the EM algorithm for the MVN CDF method.
        resulting_values <- EMAlgorithmCDF(
            initial_prob,
            iterations,
            time,
            stop_threshold,
            verbose,
            params$montecarlo_method,
            as.numeric(params$montecarlo_error),
            as.integer(params$montecarlo_iterations)
        )
        object$montecarlo_method <- params$montecarlo_method
        object$montecarlo_error <- params$montecarlo_error
        object$montecarlo_iterations <- params$montecarlo_iterations
    } else {
        # ==== MULTINOMIAL | MVN PDF | EXACT ==== #
        # Run the EM algorithm
        resulting_values <- EMAlgorithmAll(
            method,
            initial_prob,
            iterations,
            time,
            stop_threshold,
            verbose
        )

        # If the method is exact, clean the exact set
        if (method == "Exact") clean_exact_precompute()
    }
    # ---------- ... ---------- #
    # Deallocates memory from C
    clean_everything()


    object$probability <- as.matrix(resulting_values$result)
    object$logLik <- as.numeric(resulting_values$log_likelihood)
    object$iterations <- as.numeric(resulting_values$total_iterations)
    object$time <- resulting_values$total_time
    object$status <- resulting_values$stopping_reason

    invisible(object) # Updates the object.
}

#' Generic method for establishing the `compute` call.
#'
#' @param object The object to call
#' @param ... Aditional arguments for methods
#' @return Returns a list according the received object
#' @examples
#' eimObject <- eim(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#' compute(eimObject)
#' @seealso [compute.eim]
#' @noRd
#' @name computeGeneric
#' @export
compute <- function(object, ...) {
    UseMethod("compute", object)
}

#' Runs a Bootstrap to Estimate the **Standard Deviation** of Predicted Probabilities
#'
#' @description
#' This function computes the Expected-Maximization (EM) algorithm "`boot_iter`" times by sampling,
#' with replacement, all of the ballot boxes out of the original data for each
#' iteration. A new sample is drawn for every iteration, and from all of these iterations,
#' the **standard deviation** of the estimated probabilities is computed.
#'
#' @param boot_iter Integer specifying how many times to run the
#'   EM algorithm.
#'
#' @param get_p Boolean that determines if the probability is estimated aswell. Defaults to `TRUE`.
#'
#' @inheritParams compute
#'
#' @param ... Additional arguments passed to the [compute] function that will execute the EM algorithm.
#'
#' @inherit compute note
#'
#' @usage
#' bootstrap(
#'  boot_iter,
#'  object = NULL,
#'  X = NULL,
#'  W = NULL,
#'  json_path = NULL,
#' 	get_p = TRUE,
#'  ...
#' )
#'
#' @seealso The [eim] object and [compute] implementation.
#'
#' @return
#' Updates the `eim` object with a new `sd` field containing the estimated
#' standard deviations of the probabilities.
#'
#' @examples
#' # Simulate data for 20 ballot boxes, 5 candidates, 3 groups
#'
#' simulations <- simulate_elections(
#'     num_ballots = 20,
#'     num_candidates = 5,
#'     num_groups = 3,
#'     ballot_voters = rep(100, 20)
#' )
#'
#' # Create an eim object
#' model <- eim(X = simulations$X, W = simulations$W)
#'
#' # Estimate the standard deviation using 30 EM iterations with
#' # MVN PDF method
#'
#' model <- bootstrap(model, boot_iter = 30, method = "mvn_cdf")
#'
#' # The eim object now contains the standard deviations in model$sd
#' print(model$sd)
#'
#' # Probability is estimated too
#' print(model$probability)
#'
#' @name bootstrap
#' @export
bootstrap <- function(boot_iter,
                      object = NULL,
                      X = NULL,
                      W = NULL,
                      json_path = NULL,
                      get_p = TRUE,
                      ...) {
    # Check if parameters were provided with ...
    # Retrieve the default values from compute() as a list
    compute_defaults <- formals(compute.eim)
    # Update the values with user-handed parameters
    compute_args <- modifyList(as.list(compute_defaults), list(...))

    # Check if each source is provided
    object_provided <- !is.null(object)
    xw_provided <- !is.null(X) || !is.null(W)
    json_provided <- !is.null(json_path)

    # Ensure only one source is provided
    if (sum(object_provided, xw_provided, json_provided) != 1) {
        stop(
            "You must provide exactly one of the following:\n",
            "(1)\tan `eim` object (initialized with `eim`)\n",
            "(2)\t`X` and `W`\n",
            "(3)\ta `json_path`"
        )
    }

    if (!object_provided) {
        # Create the object
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    if (compute_args$method == "Exact" && boot_iter >= 5 && ncol(object$X) > 3) {
        message("Warning: The bootstraping method may take a while.")
    }

    # Array for storing intermediate results
    results_array <- vector("list", length = boot_iter)

    # Main loop
    for (i in seq_len(boot_iter)) {
        # Sample according the ballot boxes
        sample_indices <- sample(1:nrow(object$X), size = nrow(object$X), replace = TRUE)
        if (compute_args$verbose) {
            message("Running iteration:\t", i, "\nThe sampled ballot boxes indices are:\t", sample_indices)
        }
        iteration_X <- object$X[sample_indices, , drop = FALSE]
        iteration_W <- object$W[sample_indices, , drop = FALSE]

        # Create a temporary eim object and compute
        sample_object <- eim(X = iteration_X, W = iteration_W)
        compute_args$object <- sample_object
        result <- do.call(compute.eim, compute_args)
        results_array[[i]] <- result$probability

        # Deallocate memory
        rm(sample_object)
        rm(result)
    }

    # Obtain the probabilities too
    if (get_p) {
        if (compute_args$verbose) {
            message("Bootstraping finished, now estimating the probabilities.")
        }
        compute_args$object <- object
        object <- do.call(compute.eim, compute_args)
    }

    # Convert results into a 3D tensor
    final_array <- array(unlist(results_array), dim = c(ncol(object$W), ncol(object$X), boot_iter))
    sd_matrix <- apply(final_array, MARGIN = c(1, 2), FUN = sd)

    if (compute_args$verbose) {
        message("Bootstraping finished!")
        message("Estimated standard deviation matrix:")
        print(sd_matrix)
        if (get_p) {
            message("Estimated probability matrix")
            object$probability
        }
    }

    object$sd <- sd_matrix
    invisible(object) # Updates the object
}

#' @description According to the state of the algorithm (either computed or not), it prints a message with its most relevant parameters
#'
#' @return Doesn't return anything. Yields messages on the console.
#'
#' @examples
#' simulations <- simulate_elections(
#'     num_ballots = 20,
#'     num_candidates = 5,
#'     num_groups = 3,
#'     ballot_voters = rep(100, 20)
#' )
#'
#' model <- eim(simulations$X, simulations$W)
#' print(model) # Will print the X and W matrix.
#' compute(model)
#' print(model) # Will print the X and W matrix among the EM results.
#' @noRd
#' @export
print.eim <- function(object, ...) {
    cat("eim ecological inference model\n")
    # Determine if truncation is needed
    truncated_X <- (nrow(object$X) > 5)
    truncated_W <- (nrow(object$W) > 5)

    cat("Candidate matrix (X) [b x c]:\n")
    print(object$X[1:min(5, nrow(object$X)), ]) # nolint
    if (truncated_X) cat(".\n.\n.\n")

    cat("Group matrix (W) [b x g]:\n")
    print(object$W[1:min(5, nrow(object$W)), ]) # nolint
    if (truncated_W) cat(".\n.\n.\n")

    if (!is.null(object$method)) {
        cat("Estimated probability [g x c]:\n")
        truncated_P <- (nrow(object$probability))
        print(object$probability[1:min(5, nrow(object$probability)), ]) # nolint
        if (truncated_P) cat(".\n.\n.\n")
        cat("Method:\t", object$method, "\n")
        if (object$method == "Hit and Run") {
            cat("Step size (M):\t", object$step_size)
            cat("Samples (S):\t", object$samples)
        } else if (object$method == "MVN PDF") {
            cat("Montecarlo method:\t", object$montecarlo_method)
            cat("Montecarlo iterations:\t", object$montecarlo_iterations)
            cat("Montecarlo error:\t", object$montecarlo_error)
        }
        cat("Total Iterations:\t", object$total_iterations, "\n")
        cat("Total Time (s):\t", object$total_time, "\n")
        cat("Log-likelihood:\t", tail(object$logLiklihood, 1))
    }
}

#' @description Shows, in form of a list, a selection of the most important attributes. It'll retrieve the method, number of candidates, ballots, groups and total votes as well as the principal results of the EM algorithm.
#'
#' @param object An `"eim"` object.
#' @param ... Additional arguments that are ignored.
#' @return A list with the chosen attributes
#'
#' @examples
#'
#' simulations <- simulate_elections(
#'     num_ballots = 5,
#'     num_candidates = 3,
#'     num_groups = 2,
#'     ballot_voters = rep(100, 5)
#' )
#'
#' model <- eim(simulations$X, simulations$W)
#' summarised <- summary(model)
#' names(summarised)
#' # "candidates" "groups" "ballots" "votes"
#' @noRd
#' @export
summary.eim <- function(object, ...) {
    # Generates the list with the chore attribute.
    object_core_attr <- list(
        candidates = ncol(object$X),
        groups = ncol(object$W),
        ballots = ncol(object$X),
        votes = sum(object$X)
    )

    # A list with attributes to display if the EM is computed.
    if (!is.null(object$method)) {
        object_compute_attr <- list(
            method = object$method,
            probabilities = object$probability,
            logLik = object$logLik,
            status = object$status
        )
    } else {
        object_compute_attr <- list(
            status = "Not computed yet"
        )
    }

    # Display sd if the bootstrapping method has been called.
    if (!is.null(object$sd)) {
        object_compute_attr$sd <- object$sd
    }

    final_list <- c(object_core_attr, object_compute_attr)
    final_list
}

#' Returns the object estimated probability
#'
#' @param object An `"eim"` object.
#' @param ... Additional arguments that are ignored.
#' @return The probability matrix
#' @noRd
#' @export
as.matrix.eim <- function(object, ...) {
    if (is.null(object$probability)) {
        stop(paste0(
            "Probability matrix not available. Run compute().",
        ))
    }
    return(object$probability)
}

#' Save an `eim` object to a file
#'
#' This function saves an `eim` object to a specified file format. Supported formats are
#' **RDS**, **JSON**, and **CSV**. The function dynamically extracts and saves all available
#' attributes when exporting to JSON. If the `probability` field exists, it is saved when using CSV;
#' otherwise, it yields an error.
#'
#' @param object An `eim` object.
#' @param filename A character string specifying the file path, including the desired file extension (`.rds`, `.json`, or `.csv`).
#' @param ... Additional arguments (currently unused but included for compatibility).
#'
#' @usage save(object, filename, ...)
#'
#' @details
#' - If the file extension is **RDS**, the entire object is saved using `saveRDS()`.
#' - If the file extension is **JSON**, all available attributes of the object are stored in JSON format.
#' - If the file extension is **CSV**:
#'   - If the object contains a `probability` field, only that field is saved as a CSV.
#'   - Otherwise, returns an error.
#'
#' @return The function does not return anything explicitly but saves the object to the specified file.
#'
#' @seealso The [eim] object implementation.
#'
#' @examples
#'
#' model <- eim(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#'
#' compute(model)
#'
#' # Save as RDS
#' save(model, "model_results.rds")
#'
#' # Save as JSON
#' save(model, "model_results.json")
#'
#' # Save as CSV
#' save(model, "model_results.csv")
#'
#' @name save
#' @export
eim.save <- function(object, filename, ...) {
    # Ensure filename is a valid string
    if (!is.character(filename) || length(filename) != 1) {
        stop("Invalid filename. Please provide a valid file path as a character string.")
    }

    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    # Get file extension
    file_ext <- tools::file_ext(filename)

    # Save as RDS
    if (file_ext == "rds") {
        saveRDS(object, file = filename)
        message("Results saved as RDS: ", filename)

        # Save as JSON (with all available attributes)
    } else if (file_ext == "json") {
        json_data <- list()

        # Dynamically extract all attributes and store them
        for (name in names(object)) {
            json_data[[name]] <- object[[name]]
        }

        jsonlite::write_json(json_data, filename, pretty = TRUE, auto_unbox = TRUE)
        message("Results saved as JSON: ", filename)

        # Save as CSV
    } else if (file_ext == "csv") {
        if (!is.null(object$probability)) {
            write.csv(object$probability, filename, row.names = TRUE)
            message("Probability matrix saved as CSV: ", filename)
        } else {
            stop("The `compute()` method must be called for saving a '.csv' file.")
        }
    } else {
        stop("Unsupported file format. Use '.rds', '.json', or '.csv'.")
    }
}

#' @noRd
#' @export
eim.write.csv <- function(object, filename, ...) {
    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }
    if (!is.character(filename) || length(filename) != 1) {
        stop("Invalid filename. Please provide a valid file path as a character string.")
    }

    # Get file extension
    file_ext <- tools::file_ext(filename)

    if (file_ext != "csv") {
        stop("The filepath provided must end with '.csv'")
    }

    if (!is.null(object$probability)) {
        write.csv(object$probability, filename, row.names = TRUE)
        message("Probability matrix saved as CSV: ", filename)
    } else {
        stop("The `compute()` method must be called for saving a '.csv' file.")
    }
}

#' @noRd
#' @export
eim.dput <- function(object, filename, ...) {
    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }
    if (!is.character(filename) || length(filename) != 1) {
        stop("Invalid filename. Please provide a valid file path as a character string.")
    }

    # Get file extension
    file_ext <- tools::file_ext(filename)

    if (file_ext != "rda") {
        stop("The filepath provided must end with '.rda'")
    }
    saveRDS(object, file = filename)
    message("Results saved as RDS: ", filename)
}

#' Returns the log-likelihood of the last iteration
#'
#' @param object An `eim` object
#' @param ... Additional parameters that won't do nothing
#'
#' @return A numeric value with the log-likelihood
#' @noRd
#' @export
eim.logLik <- function(object, ...) {
    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    if (is.null(object$logLik)) {
        stop("The `compute()` method must be called for getting the log-likelihood.")
    }
    tail(object$logLik, 1)
}
