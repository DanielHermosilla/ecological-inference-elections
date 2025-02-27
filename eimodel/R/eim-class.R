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
#' @details
#' If `X` and `W` are directly supplied, they must match the
#' dimensions of ballot boxes `(b)`. Alternatively, if `json_path` is provided, the function expects
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
#' A way to generate synthetic data for `X` and `W` is by using the [simulate_election] function. See Example 2 below.
#'
#' @section Methods:
#' In addition to this constructor, the "eim" class provides several
#' S3 methods for common operations. Some of these methods are fully documented,
#' while others are ommited due to its straightfoward implementantion. The available methods are:
#' - [run_em] -- Runs the EM algorithm.
#' - [bootstrap] -- Estimates the standard deviation of the probabilities.
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
#' # Simulate data for 500 ballot boxes, 4 candidates and 5 groups
#' sim_result <- simulate_election(
#'     num_ballots = 500,
#'     num_candidates = 3,
#'     num_groups = 5,
#'     group_proportions = c(0.2, 0.2, 0.4, 0.1, 0.1),
#' )
#'
#' model2 <- eim(X = sim_result$X, W = sim_result$W)
#'
#' # Example 3: Create an object from a user defined matrix with 8 ballot boxes,
#' # 2 candidates and 7 groups.
#'
#' x_mat <- matrix(c(
#'     57, 90,
#'     60, 84,
#'     43, 102,
#'     72, 71,
#'     63, 94,
#'     52, 80,
#'     60, 72,
#'     54, 77,
#' ), nrow = 8, ncol = 2, byrow = TRUE)
#'
#' w_mat <- matrix(c(
#'     10, 15, 25, 21, 10, 40, 26,
#'     11, 21, 37, 32, 8, 23, 12,
#'     17, 12, 43, 27, 12, 19, 15,
#'     20, 18, 25, 15, 22, 17, 26,
#'     21, 19, 27, 16, 23, 22, 29,
#'     18, 16, 20, 14, 19, 22, 23,
#'     10, 15, 21, 18, 20, 16, 32,
#'     12, 17, 19, 22, 15, 18, 28
#' ), nrow = 8, ncol = 7, byrow = TRUE)
#'
#' model3 <- eim(X = x_mat, W = w_mat)
#'
#' @export
#' @aliases NULL
eim <- function(X = NULL, W = NULL, json_path = NULL) {
    x_provided <- !is.null(X)
    w_provided <- !is.null(W)
    xw_provided <- x_provided || w_provided
    json_provided <- !is.null(json_path)

    if (sum(x_provided, w_provided) == 1) {
        stop("eim: If providing a matrix, 'X' and 'W' must be provided.")
    }

    if (sum(xw_provided, json_provided) != 1) {
        stop(
            "eim: You must provide exactly one of the following:\n",
            "(1)\tan `eim` object (initialized with `eim`)\n",
            "(2)\t`X` and `W`\n"
        )
    }

    # Load data from JSON if a path is provided
    if (json_provided) {
        matrices <- .validate_json_eim(json_path) # nolint
        X <- as.matrix(matrices$X)
        W <- as.matrix(matrices$W)
    }

    # Perform matricial validation
    .validate_eim(X, W) # nolint

    # Convert matrices for C API. TODO: Remove this, or use R_PreserveObject on C.
    param_X <- matrix(as.numeric(X), nrow(X), ncol(X))
    param_X <- t(param_X)
    param_W <- matrix(as.numeric(W), nrow(W), ncol(W))
    # RsetParameters(param_X, param_W)

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
#' Executes the Expectation-Maximization (EM) algorithm indicating the approximation method to use in the *E*-step.
#' Certain methods may require additional arguments, which can be passed through `...` (see [fastei-package] for more details).
#'
#' @param object An object of class `eim`, which can be created using the [eim] function. This parameter should not be used if either (i) `X` and `W` matrices or (ii) `json_path` is supplied. See **Note**.
#'
#' @inheritParams eim
#'
#' @param method An optional string specifying the method used for estimating the *E*-step. Valid
#'   options are:
#' - `mult`: The default method, using a single sum of Multinomial distributions.
#' - `mvn_cdf`: Uses a Multivariate Normal CDF distribution to approximate the conditional probability.
#' - `mvn_pdf`: Uses a Multivariate Normal PDF distribution to approximate the conditional probability.
#' - `hnr`: Uses MCMC to sample vote outcomes. This is used to estimate the conditional probability of the *E*-step.
#' - `exact`: Solves the *E*-step using the Total Probability Law.
#'
#' For a detailed description of each method, see [fastei-package] and **References**.
#'
#' @param initial_prob An optional string specifying the method used to obtain the initial
#'   probability. Accepted values are:
#' - `uniform`: Assigns equal probability to every candidate within each group.
#' - `proportional`: Assigns probabilities to each group based on the proportion of candidates votes.
#' - `group_proportional`: Computes the probability matrix by taking into account both group and candidate proportions. This is the default method.
#'
#' @param mismatch Boolean, if `TRUE` allows a mismatch between the voters and votes for each ballot-box, only works if method is `mvn_cdf`, `mvn_pdf`, and `mult`. If `FALSE`, throws an error if there is a mismatch. By default it is `FALSE`.
#'
#' @param maxiter An optional integer indicating the maximum number of EM iterations.
#'   The default value is `1000`.
#'
#' @param maxtime An optional numeric specifying the maximum running time (in seconds) for the
#'   algorithm. This is checked at every iteration of the EM algorithm. The default value is `3600`, which corresponds to an hour.
#'
#' @param stop_threshold An optional numeric value indicating the minimum difference between
#'   consecutive log-likelihood values required to stop iterating. The default value is `0.001`.
#'
#' @param verbose An optional boolean indicating whether to print informational messages during the EM
#'   iterations. The default value is `FALSE`.
#'
#' @param step_size An optional integer specifying the step size for the `hnr`
#'   algorithm. This parameter is only applicable when `method = hnr` and will
#'   be ignored otherwise. The default value is `3000`.
#'
#' @param samples An optional integer indicating the number of samples to generate for the
#'   **Hit and Run** method. This parameter is only relevant when `method = hnr`.
#'   The default value is `1000`.
#'
#' @param mc_method An optional string specifying the method used to estimate the `mvn_cdf` method
#'   via a Monte Carlo simulation. Accepted values are `genz` and `genz2`, with `genz2`
#'   set as the default. This parameter is only applicable when `method = mvn_cdf`. See **References** for more details.
#'
#' @param mc_error An optional numeric value defining the error threshold for the Monte Carlo
#'   simulation when estimating the `mvn_cdf` method. The default value is `1e-6`. This parameter is only relevant
#' when `method = mvn_cdf`.
#'
#' @param mc_samples An optional integer specifying the number of Monte Carlo
#'   samples for the `mvn_cdf` method. The default value is `5000`. This argument is only applicable when `method = mvn_cdf`.
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
#' The function returns an 'eim' object with the following attributes:
#' \describe{
#'   \item{prob}{The estimated probability matrix `(g x c)`.}
#' 	 \item{cond_prob}{A `(b x g x c)` 3d-array with the probability that a at each ballot-box a voter of each group voted for each candidate, given the observed outcome at the particular ballot-box.}
#'   \item{logLik}{The log-likelihood value from the last iteration.}
#'   \item{iterations}{The total number of iterations performed by the EM algorithm.}
#'   \item{time}{The total execution time of the algorithm in seconds.}
#'   \item{status}{
#'     The final status ID of the algorithm upon completion:
#'     \itemize{
#'       \item \code{0}: Convergence achieved.
#'       \item \code{1}: Maximum time reached.
#'       \item \code{2}: Maximum iterations reached.
#'     }
#'   }
#'   \item{message}{The finishing status displayed as a message.}
#'   \item{method}{The method for estimating the conditional probability in the E-step.}
#' }
#' Aditionally, it will create `samples` and `step_size` parameters if the specified method is `hnr`, or `mc_method`, `mc_error` and `mc_samples` if the method is `mvn_cdf`.
#'
#' @usage
#' run_em(
#' 	object = NULL,
#'  X = NULL,
#'  W = NULL,
#'  json_path = NULL,
#'  method = "mult",
#'  initial_prob = "group_proportional",
#'  mismatch = FALSE,
#'  maxiter = 1000,
#'  maxtime = 3600,
#'  stop_threshold = 0.001,
#'  verbose = FALSE,
#'  ...
#' )
#'
#' @examples
#' # Example 1: Compute the Expected-Maximization with default settings
#' simulations <- simulate_elections(
#'     num_ballots = 300,
#'     num_candidates = 5,
#'     num_groups = 3,
#' )
#' model <- eim(simulations$X, simulations$W)
#' model <- run_em(model) # Returns the object with updated attributes
#'
#' # Example 2: Compute the Expected-Maximization using the Hit and Run method
#' model <- run_em(
#'     object = model,
#'     method = "hnr",
#'     step_size = 1500,
#'     samples = 800
#' )
#'
#' # Example 3: Run the MVN CDF method with default settings
#' model <- run_em(object = model, method = "mvn_cdf")
#'
#' # Example 4: Perform an Exact estimation using user-defined parameters
#' \dontrun{
#' run_em(
#'     json_path = "a/json/file.json",
#'     method = "exact",
#'     initial_prob = "uniform",
#'     maxiter = 10,
#'     maxtime = 600,
#'     stop_threshold = 1e-3,
#'     verbose = TRUE
#' )
#' }
#'
#' @name run_em
#' @export
run_em <- function(object = NULL,
                   X = NULL,
                   W = NULL,
                   json_path = NULL,
                   method = "mult",
                   initial_prob = "group_proportional",
                   mismatch = FALSE,
                   maxiter = 1000,
                   maxtime = 3600,
                   stop_threshold = 0.001,
                   verbose = FALSE, ...) {
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint

    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("run_em: The object must be initialized with the eim() function.")
    }

    # Note: Mismatch restricted methods are checked inside .validate_compute
    if (!mismatch && (sum(X) != sum(W))) {
        stop(
            "run_em: Mismatch in the number of votes: 'X' has ", sum(X),
            " votes, but 'W' has ", sum(W), " rows. To allow a mismatch, set the argument to TRUE."
        )
    }

    object$method <- method

    # Default values
    if (method == "hnr") {
        # Step size
        object$step_size <- as.integer(if ("step_size" %in% names(all_params)) all_params$step_size else 3000)
        # Samples
        object$samples <- as.integer(if ("samples" %in% names(all_params)) all_params$samples else 1000)
    } else if (method == "mvn_cdf") {
        # Montecarlo method
        object$mc_method <- if ("mc_method" %in% names(all_params)) all_params$mc_method else "genz2"
        # Montecarlo samples
        object$mc_samples <- if ("mc_samples" %in% names(all_params)) all_params$mc_samples else 5000
        # Montecarlo error
        object$mc_error <- if ("mc_error" %in% names(all_params)) all_params$mc_error else 1e-6
    }

    RsetParameters(t(object$X), object$W)

    resulting_values <- EMAlgorithmFull(
        method,
        initial_prob,
        maxiter,
        maxtime,
        stop_threshold,
        verbose,
        as.integer(if (!is.null(object$samples)) object$samples else 3000),
        as.integer(if (!is.null(object$step_size)) object$step_size else 1000),
        if (!is.null(object$mc_method)) object$mc_method else "genz2",
        as.numeric(if (!is.null(object$mc_samples)) object$mc_samples else 5000),
        as.numeric(if (!is.null(object$mc_error)) object$mc_error else 1e-6)
    )
    # ---------- ... ---------- #

    object$prob <- as.matrix(resulting_values$result)
    object$logLik <- as.numeric(resulting_values$log_likelihood)
    object$iterations <- as.numeric(resulting_values$total_iterations)
    object$time <- resulting_values$total_time
    object$message <- resulting_values$stopping_reason
    object$status <- as.integer(resulting_values$finish_id)
    object$q <- resulting_values$q

    invisible(object) # Updates the object.
}

#' Runs a Bootstrap to Estimate the **Standard Deviation** of Predicted Probabilities
#'
#' @description
#' This function computes the Expected-Maximization (EM) algorithm "`nboot`" times. It then computes the standard deviation from the `nboot` estimated probability matrices on each component.
#'
#' @param nboot Integer specifying how many times to run the
#'   EM algorithm.
#'
#' @inheritParams run_em
#'
#' @inheritParams simulate_election
#'
#' @param ... Additional arguments passed to the [run_em] function that will execute the EM algorithm.
#'
#' @inherit run_em note
#'
#' @usage
#' bootstrap(
#'  object = NULL,
#'  X = NULL,
#'  W = NULL,
#'  json_path = NULL,
#'  nboot = 50,
#'  seed = NULL
#'  ...
#' )
#'
#' @seealso The [eim] object and [run_em] implementation.
#'
#' @return
#' Returns an `eim` object with the `sd` field containing the estimated standard deviations of the probabilities. If an `eim` object is provided, its attributes (see [run_em]) are retained in the returned object.
#'
#' @examples
#' # Example 1: Using an 'eim' object directly
#' simulations <- simulate_election(
#'     num_ballots = 200,
#'     num_candidates = 5,
#'     num_groups = 3,
#' )
#'
#' model <- eim(X = simulations$X, W = simulations$W)
#'
#' model <- bootstrap(
#'     object = model,
#'     nboot = 30,
#'     method = "mvn_cdf",
#'     maxiter = 500,
#'     verbose = TRUE,
#'     mc_samples = 5000
#' )
#'
#' # Access standard deviation throughout 'model'
#' print(model$sd)
#'
#'
#' # Example 2: Providing 'X' and 'W' matrices directly
#' model <- bootstrap(
#'     X = simulations$X,
#'     W = simulations$W,
#'     nboot = 70,
#'     method = "mvn_pdf",
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
#'     json_path = "path/to/election_data.json",
#'     nboot = 70,
#'     method = "mult",
#' )
#'
#' print(model$sd)
#' }
#'
#' @name bootstrap
#' @export
bootstrap <- function(object = NULL,
                      X = NULL,
                      W = NULL,
                      json_path = NULL,
                      nboot = 50,
                      seed = NULL,
                      ...) {
    # Retrieve the default values from run_em() as a list
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint # It would validate nboot too.

    # Retrieve default values from run_em() and update with user parameters
    run_em_defaults <- formals(run_em)
    run_em_args <- modifyList(as.list(run_em_defaults), all_params)
    run_em_args <- run_em_args[names(run_em_args) != "..."] # Remove ellipsis

    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)

    # Initialize eim object if needed
    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("Bootstrap: The object must be initialized with the `eim()` function.")
    }

    # Extract parameters with defaults if missing
    method <- all_params$method %||% "mult"
    initial_prob <- all_params$initial_prob %||% "group_proportional"
    maxiter <- as.integer(all_params$maxiter %||% 1000)
    maxtime <- as.numeric(all_params$maxtime %||% 3600)
    stop_threshold <- as.numeric(all_params$stop_threshold %||% 0.01)
    verbose <- as.logical(all_params$verbose %||% FALSE)

    # Handle method-specific defaults
    if (method == "hnr") {
        step_size <- as.integer(all_params$step_size %||% 3000)
        samples <- as.integer(all_params$samples %||% 1000)
    } else if (method == "mvn_cdf") {
        mc_method <- all_params$mc_method %||% "genz2"
        mc_samples <- as.integer(all_params$mc_samples %||% 5000)
        mc_error <- as.numeric(all_params$mc_error %||% 1e-6)
    } else {
        step_size <- 0L
        samples <- 0L
        mc_method <- ""
        mc_samples <- 0L
        mc_error <- 0.0
    }

    # Call C++ bootstrap function
    result <- bootstrapAlg(
        t(object$X), object$W, as.integer(nboot),
        method, initial_prob, as.integer(maxiter),
        maxtime, stop_threshold, verbose,
        as.integer(step_size), as.integer(samples),
        mc_method, mc_error, as.integer(mc_samples)
    )

    # Compute standard deviation from the bootstrapped probability matrices
    reshaped_result <- array(result, dim = c(nrow(object$X), ncol(object$W), nboot))
    object$sd <- apply(reshaped_result, c(1, 2), sd)

    return(object)
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
#' run_em(model)
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
        truncated_P <- (nrow(object$prob) > 5)
        print(round(object$prob[1:min(5, nrow(object$prob)), ], 2)) # nolint
        if (truncated_P) cat(".\n.\n.\n")
        cat("Method:\t", object$method, "\n")
        if (object$method == "Hit and Run") {
            cat("Step size (M):", object$step_size)
            cat("Samples (S):", object$samples)
        } else if (object$method == "MVN PDF") {
            cat("Montecarlo method:", object$mc_method)
            cat("Montecarlo iterations:", object$mc_samples)
            cat("Montecarlo error:", object$mc_error)
        }
        cat("Total Iterations:", object$iterations, "\n")
        cat("Total Time (s):", object$time, "\n")
        cat("Log-likelihood:", tail(object$logLik, 1))
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
        ballots = nrow(object$X),
        votes = sum(object$X)
    )

    # A list with attributes to display if the EM is computed.
    if (!is.null(object$method)) {
        object_run_em_attr <- list(
            method = object$method,
            prob = object$prob,
            logLik = object$logLik,
            status = object$status
        )
    } else {
        object_run_em_attr <- list(
            status = -1 # Not computed yet
        )
    }

    # Display sd if the bootstrapping method has been called.
    if (!is.null(object$sd)) {
        object_run_em_attr$sd <- object$sd
    }

    final_list <- c(object_core_attr, object_run_em_attr)
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
            "Probability matrix not available. Run run_em()."
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
#' @usage save_eim(object, filename, ...)
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
#' run_em(model)
#'
#' # Save as RDS
#' save_eim(model, "model_results.rds")
#'
#' # Save as JSON
#' save_eim(model, "model_results.json")
#'
#' # Save as CSV
#' save_eim(model, "model_results.csv")
#'
#' @name save_eim
#' @aliases NULL
#' @export
save_eim <- function(object, filename, ...) {
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
            write.csv(as.matrix(object$prob), filename, row.names = TRUE)
            message("Probability matrix saved as CSV: ", filename)
        } else {
            stop("The `run_em()` method must be called for saving a '.csv' file.")
        }
    } else {
        stop("Unsupported file format. Use '.rds', '.json', or '.csv'.")
    }
}

#' @noRd
#' @export
write.csv.eim <- function(object, filename, ...) {
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
        write.csv(as.matrix(object$prob), filename, row.names = TRUE)
        message("Probability matrix saved as CSV: ", filename)
    } else {
        stop("The `run_em()` method must be called for saving a '.csv' file.")
    }
}

#' @noRd
#' @export
dput.eim <- function(object, filename, ...) {
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
logLik.eim <- function(object, ...) {
    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    if (is.null(object$logLik)) {
        stop("The `run_em()` method must be called for getting the log-likelihood.")
    }
    tail(object$logLik, 1)
}
