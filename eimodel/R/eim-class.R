library(jsonlite)

#' S3 Object for the Expectation-Maximization Algorithm
#'
#' This constructor creates an `eim` S3 object, either by using matrices
#' `X` and `W` directly or by reading them from a JSON file. Each
#' `eim` object encapsulates the data (votes for candidates and demographic
#' groups) required by the underlying Expectation-Maximization algorithm.
#'
#' @param X A `(b x c)` matrix representing candidate votes per ballot box.
#'
#' @param W A `(b x g)` matrix representing group votes per ballot box.
#'
#' @param json_path A path to a JSON file containing `X` and `W` fields, stored as nested arrays.
#'
#' @param mismatch Boolean to allow a mismatch between group and candidate votes. By default it is `FALSE`.
#'
#' @details
#' If `X` and `W` are directly supplied, they must match the correct
#' dimensions of ballot boxes (b). Alternatively, if `json_path` is provided, the function expects
#' the JSON file to contain elements named `"X"` and `"W"` under the
#' top-level object. This two approaches are **mutually exclusable**, yielding an error otherwise.
#'
#' Internally, this function also initializes the corresponding instance within
#' the low-level (C-based) API, ensuring the data is correctly registered for
#' further processing by the EM algorithm.
#'
#' @return A list of class `eim` containing:
#' \describe{
#'   \item{\code{X}}{The candidate votes matrix \code{(b x c)}.}
#'   \item{\code{W}}{The group votes matrix \code{(b x g)}.}
#' }
#'
#' @note
#' For an alternative way to generate `X` and `W`, you can use the
#' [simulate_election] function, which provides several parameters
#' (such as `group_proportions`, `lambda`, and
#' `probability`) for controlling the underlying distributions of votes.
#' Afterwards, you may pass the resulting matrices to `eim()`.
#'
#' @section Methods:
#' In addition to this constructor, the "eim" class provides several
#' S3 methods for common operations. Some of these methods are fully documented,
#' while others are ommited due to its straightfoward implementantion. The available methods are:
#' - [compute] -- Runs the EM algorithm.
#' - [bootstrap] -- Estimates the standard deviation of probabilities
#' - [save] -- Save the results to a specified file.
#' - `print` -- Print useful information about the object.
#' - `summary` -- Shows, in form of a list, the most important attributes.
#' - `as.matrix` -- Returns the probability matrix.
#' - `write.csv` -- Writes the probability matrix in a `.csv` file.
#' - `dput` -- Writes the object in a `.rda` file.
#' - `logLik` -- Returns the log-likelihood from the last iteration.
#'
#' @examples
#'
#' # Example 1: Create an eim object from a JSON file
#' \dontrun{
#' model1 <- eim(json_path = "path/to/file.json")
#' }
#'
#' # Example 2: Use simulate_election with optional parameters, then create an eim object
#' # from matrices
#'
#' # Simulate data for 50 ballot boxes, 4 candidates and 5 groups
#' sim_result <- simulate_election(
#'     num_ballots = 50,
#'     num_candidates = 4,
#'     num_groups = 5,
#'     group_proportions = c(0.2, 0.2, 0.4, 0.1, 0.1),
#' )
#'
#' model2 <- eim(X = sim_result$X, W = sim_result$W)
#'
#' # Example 3: Create an object from a user defined matrix
#'
#' x_mat <- matrix(c(
#'     10, 20, 30, 40,
#'     15, 25, 35, 45,
#'     12, 22, 32, 42
#' ), nrow = 3, ncol = 4, byrow = TRUE)
#'
#' w_mat <- matrix(c(
#'     100, 150, 200, 250, 300,
#'     110, 160, 210, 260, 310,
#'     120, 170, 220, 270, 320
#' ), nrow = 3, ncol = 5, byrow = TRUE)
#'
#' model3 <- eim(X = x_mat, W = w_mat, mismatch = TRUE)
#' @export
#' @aliases NULL
eim <- function(X = NULL, W = NULL, json_path = NULL, mismatch = FALSE) {
    # Load data from JSON if a path is provided
    if (!is.null(json_path) && nzchar(json_path)) {
        matrices <- .validate_json_eim(json_path) # nolint
        X <- as.matrix(matrices$X)
        W <- as.matrix(matrices$W)
    } else if (is.null(X) || is.null(W)) {
        stop("Providing either 'X' or 'W' and a 'json_path' is not allowed.")
    }
    # Perform matricial validation
    .validate_eim(X, W, mismatch) # nolint

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
#' @inheritParams eim
#'
#' @param object An object of class `eim`, which can be created using the [eim] function.
#'   This parameter shouldn't be called if you prefer to supply the `X` and `W` matrices or a `json_path`.
#'
#' @param method An optional string specifying the method used for estimating the E-step. Valid
#'   options are `mult`, `mvn_cdf`, `mvn_pdf`, `hnr`, or `exact`.
#'   The default is `mult`.
#'
#' @param initial_prob An optional string specifying the method used to obtain the initial
#'   probability. Accepted values are:
#' - `uniform`: Assigns equal probability to every candidate within each group.
#' - `proportional`: Assigns probabilities to each group based on the proportion of candidates votes.
#' - `group_proportional`: Computes the probability matrix by taking into account both group and candidate proportions. This is the default method.
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
#'   via a Monte Carlo simulation. Accepted values are `genz` and `genz2`, with `genz2`
#'   set as the default.
#'
#' @param monte_error An optional numeric value defining the error threshold for the Monte Carlo
#'   simulation when estimating the **MVN CDF**. The default value is `1e-6`.
#'
#' @param monte_samples An optional integer specifying the number of Monte Carlo
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
#'   \item{status}{
#'     The final status ID of the algorithm upon completion:
#'     \itemize{
#'       \item \code{0}: Convergence achieved.
#'       \item \code{1}: Log-likelihood decrease.
#'       \item \code{2}: Maximum time reached.
#'       \item \code{3}: Maximum iterations reached.
#'     }
#'   }
#'   \item{message}{The finishing status displayed as a message.}
#' }
#' Aditionally, it will create additional parameters if the method provided is `hnr` or `mvn_cdf`.
#'
#' @usage
#' compute(
#' 	object = NULL,
#'  X = NULL,
#'  W = NULL,
#'  json_path = NULL,
#'  mismatch = FALSE,
#'  method = "mult",
#'  initial_prob = "group proportional",
#'  maxiter = 1000,
#'  maxtime = 1440,
#'  stop_threshold = 0.001,
#'  verbose = FALSE,
#'  ...
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
compute <- function(object = NULL,
                    X = NULL,
                    W = NULL,
                    json_path = NULL,
                    mismatch = FALSE,
                    method = "mult",
                    initial_prob = "group proportional",
                    maxiter = 1000,
                    maxtime = 1440,
                    stop_threshold = 0.001,
                    verbose = FALSE, ...) {
    params <- list(...)

    if (is.null(object)) {
        object <- eim(X, W, json_path, mismatch)
    } else if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    # Check if the method provided is valid
    valid_methods <- c("hnr", "exact", "mvn_cdf", "mvn_pdf", "mult")

    if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
        stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    # Check for invalid mismatch methods
    if (mismatch && method %in% c("hnr", "exact")) {
        stop("Mismatch isn't allowed for computing with method '", paste0(method), "'")
    }

    object$method <- method

    if (method == "hnr") {
        # ========= HIT AND RUN ========= #
        # Validate step size -> if not provided use 3000 as default
        object$step_size <- as.integer(.validate_arg_compute(params, "step_size", 3000, int = TRUE)) # nolint
        # Validate samples -> if not provided use 1000 as default
        object$samples <- as.integer(.validate_arg_compute(params, "samples", 1000, int = TRUE)) # nolint
        # Run the EM algorithm for the Hit and Run method.
        resulting_values <- EMAlgorithmHitAndRun(
            initial_prob,
            maxiter,
            maxtime,
            stop_threshold,
            verbose,
            object$step_size,
            object$samples
        )
        # Deallocates memory from C
        clean_hr_precompute()
    } else if (method == "mvn_cdf") {
        # ========= MVN CDF ========= #
        # Check if there's a montecarlo method, otherwise, use 'genz' as default
        valid_cdf_methods <- c("genz, genz2")
        # Validate the method
        if ("monte_method" %in% names(params) && !params$monte_method %in% valid_cdf_methods) {
            stop(
                "MVN CDF:\tThe supplied method for estimating the CDF isn't supported. Must be one of: ",
                paste(valid_cdf_methods, collapse = ", ")
            )
        } else {
            # Use 'genz2' as default
            params$monte_method <- "genz2"
        }
        object$monte_method <- params$monte_method
        # Validate error -> if not provided use 1e-6 as default.
        object$monte_error <- as.integer(.validate_arg_compute(params, "monte_error", 1e-6, int = FALSE)) # nolint
        # Validate samples -> if not provided use 1e-6 as default.
        object$monte_samples <- as.numeric(.validate_arg_compute(params, "monte_samples", 5000, int = TRUE))

        # Run the EM algorithm for the MVN CDF method.
        resulting_values <- EMAlgorithmCDF(
            initial_prob,
            maxiter,
            maxtime,
            stop_threshold,
            verbose,
            object$monte_method,
            object$monte_error,
            object$monte_samples
        )
    } else {
        # ==== MULTINOMIAL | MVN PDF | EXACT ==== #
        # Run the EM algorithm
        resulting_values <- EMAlgorithmAll(
            method,
            initial_prob,
            maxiter,
            maxtime,
            stop_threshold,
            verbose
        )

        # If the method is exact, clean the exact set
        if (method == "Exact") clean_exact_precompute()
    }
    # ---------- ... ---------- #
    # Deallocates memory from C
    clean_everything()


    object$prob <- as.matrix(resulting_values$result)
    object$logLik <- as.numeric(resulting_values$log_likelihood)
    object$iterations <- as.numeric(resulting_values$total_iterations)
    object$time <- resulting_values$total_time
    object$message <- resulting_values$stopping_reason
    object$status <- as.integer(resulting_values$finish_id)

    invisible(object) # Updates the object.
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
#'  mismatch = FALSE,
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
#' # Example 1: Using an 'eim' object directly
#' simulations <- simulate_election(
#'     num_ballots = 20,
#'     num_candidates = 5,
#'     num_groups = 3,
#'     ballot_voters = rep(100, 20)
#' )
#'
#' # Initialize the 'eim' object
#' model <- eim(X = simulations$X, W = simulations$W)
#'
#' # Run bootstrap with an existing 'eim' object
#' model <- bootstrap(
#'     boot_iter = 30,
#'     object = model,
#'     method = "mvn_cdf",
#'     maxiter = 500,
#'     verbose = TRUE,
#'     monte_samples = 5000
#' )
#'
#' # Access standard deviation and probability results
#' print(model$sd)
#' print(model$probability)
#'
#'
#' # Example 2: Providing 'X' and 'W' matrices directly
#' model <- bootstrap(
#'     boot_iter = 50,
#'     get_p = FALSE,
#'     X = simulations$X,
#'     W = simulations$W,
#'     method = "exact",
#'     maxiter = 100,
#'     maxtime = 15,
#'     stop_threshold = 0.01
#' )
#'
#' print(model$sd)
#'
#' # Example 3: Using a JSON file as input
#' \dontrun{
#' model <- bootstrap(
#'     boot_iter = 20,
#'     json_path = "path/to/election_data.json",
#'     method = "mult",
#'     get_p = FALSE
#' )
#'
#' print(model$sd)
#' }
#'
#' @name bootstrap
#' @export
bootstrap <- function(boot_iter,
                      object = NULL,
                      X = NULL,
                      W = NULL,
                      json_path = NULL,
                      mismatch = FALSE,
                      get_p = TRUE,
                      ...) {
    # Check if parameters were provided with ...
    # Retrieve the default values from compute() as a list
    compute_defaults <- formals(compute)
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
        result <- do.call(compute, compute_args)
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
        object <- do.call(compute, compute_args)
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
            object$prob
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
        truncated_P <- (nrow(object$prob))
        print(object$prob[1:min(5, nrow(object$prob)), ]) # nolint
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
            probabilities = object$prob,
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
    if (is.null(object$prob)) {
        stop(paste0(
            "Probability matrix not available. Run compute().",
        ))
    }
    return(object$prob)
}

#' Save an `eim` object to a file
#'
#' This function saves an `eim` object to a specified file format. Supported formats are
#' **RDS**, **JSON**, and **CSV**. The function dynamically extracts and saves all available
#' attributes when exporting to JSON. If the `prob` field exists, it is saved when using CSV;
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
#'   - If the object contains a `prob` field, only that field is saved as a CSV.
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
#' @aliases NULL
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
        if (!is.null(object$prob)) {
            write.csv(object$prob, filename, row.names = TRUE)
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

    if (!is.null(object$prob)) {
        write.csv(object$prob, filename, row.names = TRUE)
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
