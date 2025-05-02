library(jsonlite)

#' S3 Object for the Expectation-Maximization Algorithm
#'
#' This constructor creates an `eim` S3 object, either by using matrices
#' `X` and `W` directly or by reading them from a JSON file. Each
#' `eim` object encapsulates the data (votes for candidates and demographic groups) required by the underlying Expectation-Maximization algorithm.
#'
#' @param X A `(b x c)` matrix representing candidate votes per ballot box.
#'
#' @param W A `(b x g)` matrix representing group votes per ballot box.
#'
#' @param json_path A path to a JSON file containing `X` and `W` fields, stored as nested arrays. It may contain additional fields with other attributes, which will be added to the returned object.
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
#'
#' \itemize{
#'   \item \code{\link{run_em}} - Runs the EM algorithm.
#'   \item \code{\link{bootstrap}} - Estimates the standard deviation.
#'   \item \code{\link{save_eim}} - Saves the object to a file.
#'   \item \code{\link{get_agg_proxy}} - Estimates an ideal group aggregation given their standard deviations.
#'   \item \code{\link{get_agg_opt}} - Estimates an ideal group aggregation among all combinations, given the log-likelihood.
#'   \item \code{print.eim} - Print info about the object.
#'   \item \code{summary.eim} - Summarize the object.
#'   \item \code{as.matrix.eim} - Returns the probability matrix.
#'   \item \code{logLik.eim} - Returns the final log-likelihood.
#' }
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
#'     group_proportions = c(0.2, 0.2, 0.4, 0.1, 0.1)
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
#'     54, 77
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
#' @aliases eim()
eim <- function(X = NULL, W = NULL, json_path = NULL) {
    x_provided <- !is.null(X)
    w_provided <- !is.null(W)
    xw_provided <- x_provided || w_provided
    json_provided <- !is.null(json_path)

    if (sum(x_provided, w_provided) == 1) {
        stop("eim: If providing a matrix, 'X' and 'W' must be provided.")
    }

    if (sum(xw_provided, json_provided) != 1) {
        stop(paste(
            "eim: You must provide exactly one of the following:\n",
            "(1)\tAn `eim` object (initialized with `eim`)\n",
            "(2)\t`X` and `W`"
        ))
    }
    extra_params <- list()

    # Load data from JSON if a path is provided
    if (json_provided) {
        matrices <- .validate_json_eim(json_path) # nolint
        X <- as.matrix(matrices$X)
        W <- as.matrix(matrices$W)
        allowed_params <- c(
            "prob",
            "initial_prob",
            "iterations",
            "nboot",
            "logLik",
            "convergence",
            "maxiter",
            "maxtime",
            "ll_threshold",
            "cond_prob",
            "sd",
            "group_agg",
            "message",
            "status",
            "time",
            "method",
            "W_agg",
            "samples",
            "step_size",
            "mc_method",
            "mc_samples",
            "mc_error"
        )
        extra_params <- matrices[names(matrices) %in% allowed_params] # TODO: Validate them
    }

    # Perform matricial validation
    .validate_eim(X, W) # nolint

    # Create the S3 object
    obj <- list(
        X = X,
        W = W
    )

    # Add optional parameters if they exist
    if (length(extra_params) > 0) {
        obj <- c(obj, extra_params)
    }

    class(obj) <- "eim"
    obj
}
#' @title Compute the Expected-Maximization Algorithm
#'
#' @description
#' Executes the Expectation-Maximization (EM) algorithm indicating the approximation method to use in the E-step.
#' Certain methods may require additional arguments, which can be passed through `...` (see [fastei-package] for more details).
#'
#' @param object An object of class `eim`, which can be created using the [eim] function. This parameter should not be used if either (i) `X` and `W` matrices or (ii) `json_path` is supplied. See **Note**.
#'
#' @inheritParams eim
#'
#' @param method An optional string specifying the method used for estimating the E-step. Valid
#'   options are:
#' - `mult`: The default method, using a single sum of Multinomial distributions.
#' - `mvn_cdf`: Uses a Multivariate Normal CDF distribution to approximate the conditional probability.
#' - `mvn_pdf`: Uses a Multivariate Normal PDF distribution to approximate the conditional probability.
#' - `mcmc`: Uses MCMC to sample vote outcomes. This is used to estimate the conditional probability of the E-step.
#' - `exact`: Solves the E-step using the Total Probability Law.
#'
#' For a detailed description of each method, see [fastei-package] and **References**.
#'
#' @param initial_prob An optional string specifying the method used to obtain the initial
#'   probability. Accepted values are:
#' - `uniform`: Assigns equal probability to every candidate within each group.
#' - `proportional`: Assigns probabilities to each group based on the proportion of candidates votes.
#' - `group_proportional`: Computes the probability matrix by taking into account both group and candidate proportions. This is the default method.
#' - `random`: Use randomized values to fill the probability matrix.
#'
#' @param allow_mismatch Boolean, if `TRUE`, allows a mismatch between the voters and votes for each ballot-box, only works if `method` is `"mvn_cdf"`, `"mvn_pdf"`, `"mult"` and `"mcmc"`. If `FALSE`, throws an error if there is a mismatch. By default it is `TRUE`.
#'
#' @param maxiter An optional integer indicating the maximum number of EM iterations.
#'   The default value is `1000`.
#'
#' @param maxtime An optional numeric specifying the maximum running time (in seconds) for the
#'   algorithm. This is checked at every iteration of the EM algorithm. The default value is `3600`, which corresponds to an hour.
#'
#' @param param_threshold An optional numeric value indicating the minimum difference between
#'   consecutive probability values required to stop iterating. The default value is `0.001`. Note that the algorithm will stop if either `ll_threshold` **or** `param_threshold` is accomplished.
#'
#' @param ll_threshold An optional numeric value indicating the minimum difference between consecutive log-likelihood values to stop iterating. The default value is `inf`, essentially deactivating
#' the threshold. Note that the algorithm will stop if either `ll_threshold` **or** `param_threshold` is accomplished.
#'
#' @param verbose An optional boolean indicating whether to print informational messages during the EM
#'   iterations. The default value is `FALSE`.
#'
#' @param seed An optional integer indicating the random seed for the randomized algorithms. This argument is only applicable if `initial_prob = "random"` or `method` is either `"mcmc"` or `"mvn_cdf"`.
#'
#' @param step_size An optional integer specifying the step size for the `mcmc`
#'   algorithm. This parameter is only applicable when `method = "mcmc"` and will
#'   be ignored otherwise. The default value is `3000`.
#'
#' @param samples An optional integer indicating the number of samples to generate for the
#'   **MCMC** method. This parameter is only relevant when `method = "mcmc"`.
#'   The default value is `1000`.
#'
#' @param mc_method An optional string specifying the method used to estimate the `mvn_cdf` method
#'   via a Monte Carlo simulation. Accepted values are `genz` and `genz2`, with `genz2`
#'   set as the default. This parameter is only applicable when `method = "mvn_cdf"`. See **References** for more details.
#'
#' @param mc_error An optional numeric value defining the error threshold for the Monte Carlo
#'   simulation when estimating the `mvn_cdf` method. The default value is `1e-6`. This parameter is only relevant
#' when `method = "mvn_cdf"`.
#'
#' @param mc_samples An optional integer specifying the number of Monte Carlo
#'   samples for the `mvn_cdf` method. The default value is `5000`. This argument is only applicable when `method = "mvn_cdf"`.
#'
#' @param ... Added for compability
#'
#' @references
#' [Thraves, C. and Ubilla, P.: *"Fast Ecological Inference Algorithm for the RxC Case"*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834). Aditionally, the MVN CDF is computed by the methods introduced in [Genz, A. (2000). Numerical computation of multivariate normal probabilities. *Journal of Computational and Graphical Statistics*](https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities)
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
#' The function returns an `eim` object with the function arguments and the following attributes:
#' \describe{
#'   \item{prob}{The estimated probability matrix `(g x c)`.}
#' 	 \item{cond_prob}{A `(b x g x c)` 3d-array with the probability that a at each ballot-box a voter of each group voted for each candidate, given the observed outcome at the particular ballot-box.}
#'   \item{logLik}{The log-likelihood value from the last iteration.}
#'   \item{iterations}{The total number of iterations performed by the EM algorithm.}
#'   \item{time}{The total execution time of the algorithm in seconds.}
#'   \item{status}{
#'     The final status ID of the algorithm upon completion:
#'     \itemize{
#'       \item \code{0}: Converged
#'       \item \code{1}: Maximum time reached.
#'       \item \code{2}: Maximum iterations reached.
#'     }
#'   }
#'   \item{message}{The finishing status displayed as a message, matching the status ID value.}
#'   \item{method}{The method for estimating the conditional probability in the E-step.}
#' }
#' Aditionally, it will create `samples` and `step_size` parameters if the specified `method = "mcmc"`, or `mc_method`, `mc_error` and `mc_samples` if `method = "mvn_cdf"`.
#'
#' Also, if the eim object supplied is created with the function [simulate_election], it also returns the real probability with the name `real_prob`. See [simulate_election].
#'
#' @examples
#' # Example 1: Compute the Expected-Maximization with default settings
#' simulations <- simulate_election(
#'     num_ballots = 300,
#'     num_candidates = 5,
#'     num_groups = 3,
#' )
#' model <- eim(simulations$X, simulations$W)
#' model <- run_em(model) # Returns the object with updated attributes
#'
#' # Example 2: Compute the Expected-Maximization using the MCMC method
#' model <- run_em(
#'     object = model,
#'     method = "mcmc",
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
#'     param_threshold = 1e-3,
#'     ll_threshold = 1e-5,
#'     verbose = TRUE
#' )
#' }
#'
#' @name run_em
#' @aliases run_em()
#' @export
run_em <- function(object = NULL,
                   X = NULL,
                   W = NULL,
                   json_path = NULL,
                   method = "mult",
                   initial_prob = "group_proportional",
                   allow_mismatch = TRUE,
                   maxiter = 1000,
                   maxtime = 3600,
                   param_threshold = 0.001,
                   ll_threshold = as.double(-Inf),
                   seed = NULL,
                   verbose = FALSE,
                   samples = 1000,
                   step_size = 3000,
                   mc_method = "genz2",
                   mc_error = 1e-5,
                   mc_samples = 5000,
                   ...) {
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("run_em: The object must be initialized with the eim() function.")
    }

    # Note: Mismatch restricted methods are checked inside .validate_compute
    mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))
    if (!allow_mismatch && length(mismatch_rows) > 0) {
        stop(
            "run_em: Row-wise mismatch in vote totals detected.\n",
            "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
            "To allow mismatches, set `allow_mismatch = TRUE`."
        )
    } else if (method == "exact" && length(mismatch_rows) > 0) {
        stop("run_em: Exact method isn't supported with mismatch")
    }

    object$method <- method

    # Default values
    if (method == "mcmc") {
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
        param_threshold,
        ll_threshold,
        verbose,
        as.integer(if (!is.null(object$step_size)) object$step_size else 1000),
        as.integer(if (!is.null(object$samples)) object$samples else 3000),
        if (!is.null(object$mc_method)) object$mc_method else "genz2",
        as.numeric(if (!is.null(object$mc_samples)) object$mc_samples else 5000),
        as.numeric(if (!is.null(object$mc_error)) object$mc_error else 1e-6)
    )
    # ---------- ... ---------- #

    object$prob <- as.matrix(resulting_values$result)
    dimnames(object$prob) <- list(colnames(object$W), colnames(object$X))
    object$iterations <- as.numeric(resulting_values$total_iterations)
    object$logLik <- as.numeric(resulting_values$log_likelihood[length(resulting_values$log_likelihood)])
    object$time <- resulting_values$total_time
    object$message <- resulting_values$stopping_reason
    object$status <- as.integer(resulting_values$finish_id)
    object$cond_prob <- resulting_values$q
    object$cond_prob <- aperm(resulting_values$q, perm = c(3, 2, 1)) # Correct dimensions
    dimnames(object$cond_prob) <- list(
        colnames(object$X),
        colnames(object$W),
        rownames(object$X)
    )
    # Add function arguments
    object$maxiter <- maxiter
    object$maxiter <- maxtime
    object$param_threshold <- param_threshold
    object$ll_threshold <- ll_threshold
    object$initial_prob <- initial_prob

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
#' @param seed An optional integer indicating the random seed for the randomized algorithms. This argument is only applicable if `initial_prob = "random"` or `method` is either `"mcmc"` or `"mvn_cdf"`. Aditionally, it sets the random draws of the ballot boxes.
#'
#' @param ... Additional arguments passed to the [run_em] function that will execute the EM algorithm.
#'
#' @inherit run_em note
#'
#'
#' @seealso The [eim] object and [run_em] implementation.
#'
#' @return
#' Returns an `eim` object with the `sd` field containing the estimated standard deviations of the probabilities and the amount of iterations that were made. If an `eim` object is provided, its attributes (see [run_em]) are retained in the returned object.
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
#'     verbose = FALSE,
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
#'     param_threshold = 0.01,
#'     allow_mismatch = FALSE
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
#' @aliases bootstrap()
#' @export
bootstrap <- function(object = NULL,
                      X = NULL,
                      W = NULL,
                      json_path = NULL,
                      nboot = 50,
                      allow_mismatch = TRUE,
                      seed = NULL,
                      ...) {
    # Retrieve the default values from run_em() as a list
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint # It would validate nboot too.

    # Initialize eim object if needed
    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("Bootstrap: The object must be initialized with the `eim()` function.")
    }

    # I need to define the method before on this case
    method <- if (!is.null(all_params$method)) all_params$method else "mult"
    # Note: Mismatch restricted methods are checked inside .validate_compute
    if (!allow_mismatch) {
        mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))

        if (length(mismatch_rows) > 0) {
            stop(
                "bootstrap: Row-wise mismatch in vote totals detected.\n",
                "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
                "To allow mismatches, set `allow_mismatch = TRUE`."
            )
        }
    } else {
        if (method == "exact") stop("run_em: Exact method isn't supported with mismatch")
    }
    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)

    # Extract parameters with defaults if missing
    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.01
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE

    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    step_size <- 0L
    samples <- 0L
    mc_method <- ""
    mc_samples <- 0L
    mc_error <- 0.0

    if (method == "mcmc") {
        step_size <- if (!is.null(all_params$step_size)) all_params$step_size else 3000
        samples <- if (!is.null(all_params$samples)) all_params$samples else 1000
        mc_method <- ""
        mc_samples <- 0L
        mc_error <- 0.0
    } else if (method == "mvn_cdf") {
        mc_method <- if (!is.null(all_params$mc_method)) all_params$method else "genz2"
        mc_samples <- if (!is.null(all_params$mc_samples)) all_params$mc_samples else 5000
        mc_error <- if (!is.null(all_params$mc_error)) all_params$mc_error else 1e-6
        step_size <- 0L
        samples <- 0L
    }

    # Call C bootstrap function
    result <- bootstrapAlg(
        t(object$X),
        object$W,
        as.integer(nboot),
        as.character(method),
        as.character(initial_prob),
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(verbose),
        as.integer(step_size),
        as.integer(samples),
        as.character(mc_method),
        as.double(mc_error),
        as.integer(mc_samples)
    )

    object$sd <- result
    dimnames(object$sd) <- list(colnames(object$W), colnames(object$X))
    object$nboot <- nboot
    object$sd[object$sd == 9999] <- Inf

    class(object) <- "eim"
    return(object)
}

#' Runs the EM algorithm aggregating adjacent groups, maximizing the variability of macro-group allocation in ballot boxes.
#'
#' This function estimates the voting probabilities (computed using [run_em]) aggregating adjacent groups so that the estimated probabilities' standard deviation (computed using [bootstrap]) is below a given threshold. See **Details** for more information.
#'
#' Groups need to have an order relation so that adjacent groups can be merged. Groups of consecutive column indices in the matrix W are considered adjacent. For example, consider the following seven groups defined by voters' age ranges: 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, and 80+. A possible group aggregation can be a macro-group composed of the three following age ranges: 20-39, 40-59, and 60+. Since there are multiple group aggregations, even for a fixed number of macro-groups, a Dynamic Program (DP) mechanism is used to find the group aggregation that maximizes the sum of the standard deviation of the macro-groups proportions among ballot boxes for a specific number of macro-groups. If no group aggregation standard deviation statistic meets the threshold condition, `NULL` is returned.
#'
#' To find the best group aggregation, the function runs the DP iteratively, starting with all groups (this case is trivial since the group aggregation is such that all macro-groups match exactly the original groups). If the standard deviation statistic (`sd_statistic`) is below the threshold (`sd_threshold`), it stops. Otherwise, it runs the DP such that the number of macro-groups is one unit less than the original number of macro-groups. If the standard deviation statistic is below the threshold, it stops. This continues until either the algorithm stops, or until no group aggregation obtained by the DP satisfies the threshold condition. If the former holds, then the last group aggregation obtained (before stopping) is returned; while if the latter holds, then no output is returned unless the user sets the input parameter `feasible=FALSE`, in which case it returns the group aggregation that has the least standard deviation statistic, among the group-aggregations obtained from the DP.
#'
#' @param object An object of class `eim`, which can be created using the [eim] function. This parameter should not be used if either (i) `X` and `W` matrices or (ii) `json_path` is supplied. See **Note** in [run_em].
#'
#' @param sd_statistic String indicates the statistic for the standard deviation `(g x c)` matrix for the stopping condition, i.e., the algorithm stops when the statistic is below the threshold. It can take the value `maximum`, in which case computes the maximum over the standard deviation matrix, or `average`, in which case computes the average.
#'
#' @param sd_threshold Numeric with the value to use as a threshold for the statistic (`sc_statistic`) of the standard deviation of the estimated probabilities. Defaults to 0.05.
#'
#' @param method An optional string specifying the method used for estimating the E-step. Valid
#'   options are:
#' - `mult`: The default method, using a single sum of Multinomial distributions.
#' - `mvn_cdf`: Uses a Multivariate Normal CDF distribution to approximate the conditional probability.
#' - `mvn_pdf`: Uses a Multivariate Normal PDF distribution to approximate the conditional probability.
#' - `mcmc`: Uses MCMC to sample vote outcomes. This is used to estimate the conditional probability of the E-step.
#' - `exact`: Solves the E-step using the Total Probability Law.
#'
#' @param feasible Logical indicating whether the returned matrix must strictly satisfy the `sd_threshold`.
#' If `TRUE`, no output is returned if the method does not find a group aggregation whose standard deviation statistic is below the threshold. If `FALSE` and the latter holds, it returns the group aggregation obtained from the DP with the the lowest standard deviation statistic. See **Details** for more information. Default is `TRUE`.
#'
#' @inheritParams bootstrap
#'
#' @param ... Additional arguments passed to the [run_em] function that will execute the EM algorithm.
#'
#' @seealso The [eim] object and [run_em] implementation.
#'
#' @return
#' It returns an eim object with the same attributes as the output of [run_em], plus the attributes:
#'
#' - **sd**: A `(a x c)` matrix with the standard deviation of the estimated probabilities computed with bootstrapping. Note that `a` denotes the number of macro-groups of the resulting group aggregation, it should be between `1` and `g`.
#' - **nboot**: Number of samples used for the [bootstrap] method.
#' - **seed**: Random seed used (if specified).
#' - **sd_statistic**: The statistic used as input.
#' - **sd_threshold**: The threshold used as input.
#' - **is_feasible**:  Boolean indicating whether the statistic of the standard deviation matrix is below the threshold.
#' - **group_agg**: Vector with the resulting group aggregation. See **Examples** for more details.
#'
#' Additionally, it will create the `W_agg` attribute with the aggregated groups, along with the attributes corresponding to running [run_em] with the aggregated groups.
#'
#' @examples
#' # Example 1: Using a simulated instance
#' simulations <- simulate_election(
#'     num_ballots = 400,
#'     num_candidates = 3,
#'     num_groups = 6,
#'     group_proportions = c(0.4, 0.1, 0.1, 0.1, 0.2, 0.1),
#'     lambda = 0.7,
#'     seed = 42
#' )
#'
#' result <- get_agg_proxy(
#'     X = simulations$X,
#'     W = simulations$W,
#'     sd_threshold = 0.015,
#'     seed = 42
#' )
#'
#' result$group_agg # c(2 6)
#' # This means that the resulting group aggregation is conformed by
#' # two macro-groups: one that has the original groups 1 and 2; and
#' # a second that has the original groups 3, 4, 5, and 6:
#' # {[1, 2], [3, 6]}
#'
#' # Example 2: Using the chilean election results
#' data(chile_election_2021)
#'
#' niebla_df <- chile_election_2021[chile_election_2021$ELECTORAL.DISTRICT == "NIEBLA", ]
#'
#' # Create the X matrix with selected columns
#' X <- as.matrix(niebla_df[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7")])
#'
#' # Create the W matrix with selected columns
#' W <- as.matrix(niebla_df[, c(
#'     "X18.19", "X20.29",
#'     "X30.39", "X40.49",
#'     "X50.59", "X60.69",
#'     "X70.79", "X80."
#' )])
#'
#' solution <- get_agg_proxy(
#'     X = X, W = W,
#'     allow_mismatch = TRUE, sd_threshold = 0.03,
#'     sd_statistic = "average", nboot = 100, seed = 42
#' )
#'
#' solution$group_agg # c(3, 4, 5, 6, 8)
#' # This means that the resulting group aggregation consists of
#' # five macro-groups: one that includes the original groups 1, 2, and 3;
#' # three singleton groups (4, 5, and 6); and one macro-group that includes groups 7 and 8.
#' # {[1, 2, 3], [4], [5], [6], [7, 8]}
#'
#' @export
get_agg_proxy <- function(object = NULL,
                          X = NULL,
                          W = NULL,
                          json_path = NULL,
                          sd_statistic = "maximum",
                          sd_threshold = 0.05,
                          method = "mult",
                          feasible = TRUE,
                          nboot = 50,
                          allow_mismatch = TRUE,
                          seed = NULL, ...) {
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

    # I need to define the method before
    method <- if (!is.null(all_params$method)) all_params$method else "mult"
    # Note: Mismatch restricted methods are checked inside .validate_compute
    if (!allow_mismatch) {
        mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))

        if (length(mismatch_rows) > 0) {
            stop(
                "get_agg_proxy: Row-wise mismatch in vote totals detected.\n",
                "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
                "To allow mismatches, set `allow_mismatch = TRUE`."
            )
        }
    } else {
        if (method == "exact") stop("run_em: Exact method isn't supported with mismatch")
    }

    # Extract parameters with defaults if missing
    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.01
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE

    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    step_size <- 0L
    samples <- 0L
    mc_method <- ""
    mc_samples <- 0L
    mc_error <- 0.0

    if (method == "mcmc") {
        step_size <- if (!is.null(all_params$step_size)) all_params$step_size else 3000
        samples <- if (!is.null(all_params$samples)) all_params$samples else 1000
        mc_method <- ""
        mc_samples <- 0L
        mc_error <- 0.0
    } else if (method == "mvn_cdf") {
        mc_method <- if (!is.null(all_params$mc_method)) all_params$method else "genz2"
        mc_samples <- if (!is.null(all_params$mc_samples)) all_params$mc_samples else 5000
        mc_error <- if (!is.null(all_params$mc_error)) all_params$mc_error else 1e-6
        step_size <- 0L
        samples <- 0L
    }

    result <- groupAgg(
        as.character(sd_statistic),
        as.double(sd_threshold),
        as.logical(feasible),
        t(object$X),
        object$W,
        as.integer(nboot),
        as.character(method),
        as.character(initial_prob),
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(verbose),
        as.integer(step_size),
        as.integer(samples),
        as.character(mc_method),
        as.double(mc_error),
        as.integer(mc_samples)
    )

    # If the returned matrix isn't the best non-feasible result
    # best_result <- TRUE if it's unfeasible
    if (!result$best_result || !feasible) {
        # Convert the 'W' matrix by merging columns
        # We add '2' to indices since it's originally 0-based.
        col_groups <- split(seq_len(ncol(object$W)), findInterval(seq_len(ncol(object$W)), c(1, result$indices + 2)))
        # Lambda function to add the columns, if there wasn't a group aggregation, return object$W
        # object$W_agg <- as.matrix(sapply(col_groups, function(cols) rowSums(object$W[, cols, drop = FALSE])))
        object$W_agg <- as.matrix(do.call(cbind, lapply(col_groups, function(cols) rowSums(object$W[, cols, drop = FALSE]))))
        rownames(object$W_agg) <- rownames(object$W) # Ballot boxes
        # if (result$indices[1] != -1) {
        object$group_agg <- unique(result$indices + 1)
        # } # Use R's index system
        object$sd <- result$bootstrap_result
        dimnames(object$sd) <- list(colnames(object$W_agg), colnames(object$X))
        object$sd[object$sd == 9999] <- Inf
        object$sd_statistic <- sd_statistic
        object$sd_threshold <- sd_threshold
        object$is_feasible <- !result$best_result

        # Add EM arguments aswell
        # core_args <- all_params[!names(all_params) %in% c("object", "X", "W", "X2", "W2", "json_path", "verbose")]
        # final_args <- c(core_args, list(object = NULL, json_path = NULL, X = object$X, W = object$W_agg, verbose = FALSE))
        em_results <- run_em(X = object$X, W = object$W_agg, method = method, allow_mismatch = TRUE)
        object <- c(
            object,
            em_results[setdiff(names(em_results), names(object))]
        )
    }

    class(object) <- "eim"
    return(object)
}

#' Runs the EM algorithm **over all possible group aggregating**, returning the one with higher likelihood while constraining the standard deviation of the probabilities.
#'
#' This function estimates the voting probabilities (computed using [run_em]) by trying all group aggregations (of adjacent groups), choosing
#' the one that achieves the higher likelihood as long as the standard deviation (computed using [bootstrap]) of the estimated probabilities
#' is below a given threshold. See **Details** for more informacion on adjacent groups.
#'
#' Groups of consecutive column indices in the matrix `W` are considered adjacent. For example, consider the following seven groups defined by voters' age
#' ranges: 20-29, 30-39, 40-49, 50-59, 60-69, 70-79, and 80+. A possible group aggregation can be a macro-group composed of the three following age
#' ranges: 20-39, 40-59, and 60+. Since there are multiple group aggregations, the method evaluates all possible group aggregations (merging only adjacent groups).
#'
#' @inheritParams get_agg_proxy
#'
#' @param ... Additional arguments passed to the [run_em] function that will execute the EM algorithm.
#'
#' @return
#' It returns an eim object with the same attributes as the output of [run_em], plus the attributes:
#'
#' - **sd**: A `(a x c)` matrix with the standard deviation of the estimated probabilities computed with bootstrapping. Note that `a` denotes the number of macro-groups of the resulting group aggregation, it should be between `1` and `g`.
#' - **nboot**: Number of samples used for the [bootstrap] method.
#' - **seed**: Random seed used (if specified).
#' - **sd_statistic**: The statistic used as input.
#' - **sd_threshold**: The threshold used as input.
#' - **group_agg**: Vector with the resulting group aggregation. See **Examples** for more details.
#'
#' Additionally, it will create the `W_agg` attribute with the aggregated groups, along with the attributes corresponding to running [run_em] with the aggregated groups.
#'
#' @examples
#' # Example 1: Using a simulated instance
#' simulations <- simulate_election(
#'     num_ballots = 20,
#'     num_candidates = 3,
#'     num_groups = 8,
#'     seed = 42
#' )
#'
#' result <- get_agg_opt(
#'     X = simulations$X,
#'     W = simulations$W,
#'     sd_threshold = 0.05,
#'     seed = 42
#' )
#'
#' result$group_agg # c(3,8)
#' # This means that the resulting group aggregation consists of
#' # two macro-groups: one that includes the original groups 1, 2, and 3;
#' # the remaining one with groups 4, 5, 6, 7 and 8.
#' # {[1, 2, 3], [4, 5, 6, 7, 8]}
#' \dontrun{
#' # Example 2: Getting an unfeasible result
#' result2 <- get_agg_opt(
#'     X = simulations$X,
#'     W = simulations$W,
#'     sd_threshold = 0.001
#' )
#'
#' result2$group_agg # Error
#' result2$X # Input candidates' vote matrix
#' result2$W # Input group-level voter matrix
#' }
#' @export
get_agg_opt <- function(object = NULL,
                        X = NULL,
                        W = NULL,
                        json_path = NULL,
                        sd_statistic = "maximum",
                        sd_threshold = 0.05,
                        method = "mult",
                        nboot = 50,
                        allow_mismatch = TRUE,
                        seed = NULL,
                        ...) {
    # Retrieve the default values from run_em() as a list
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint # It would validate nboot too.

    # Retrieve default values from run_em() and update with user parameters
    run_em_defaults <- formals(run_em)
    run_em_args <- modifyList(as.list(run_em_defaults), all_params)
    run_em_args <- run_em_args[names(run_em_args) != "..."] # Remove ellipsis

    # Initialize eim object if needed
    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("get_agg_opt: The object must be initialized with the `eim()` function.")
    }

    # Note: Mismatch restricted methods are checked inside .validate_compute
    # Method needs to be defined before
    method <- if (!is.null(all_params$method)) all_params$method else "mult"
    if (!allow_mismatch) {
        mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))

        if (length(mismatch_rows) > 0) {
            stop(
                "get_agg_opt: Row-wise mismatch in vote totals detected.\n",
                "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
                "To allow mismatches, set `allow_mismatch = TRUE`."
            )
        }
    } else {
        if (method == "exact") stop("run_em: Exact method isn't supported with mismatch")
    }

    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.01
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE

    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    step_size <- 0L
    samples <- 0L
    mc_method <- ""
    mc_samples <- 0L
    mc_error <- 0.0

    if (method == "mcmc") {
        step_size <- if (!is.null(all_params$step_size)) all_params$step_size else 3000
        samples <- if (!is.null(all_params$samples)) all_params$samples else 1000
        mc_method <- ""
        mc_samples <- 0L
        mc_error <- 0.0
    } else if (method == "mvn_cdf") {
        mc_method <- if (!is.null(all_params$mc_method)) all_params$mc_method else "genz2"
        mc_samples <- if (!is.null(all_params$mc_samples)) all_params$mc_samples else 5000
        mc_error <- if (!is.null(all_params$mc_error)) all_params$mc_error else 1e-6
        step_size <- 0L
        samples <- 0L
    }

    result <- groupAggGreedy(
        as.character(sd_statistic),
        as.double(sd_threshold),
        t(object$X),
        object$W,
        as.integer(nboot),
        as.character(method),
        as.character(initial_prob),
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(verbose),
        as.integer(step_size),
        as.integer(samples),
        as.character(mc_method),
        as.double(mc_error),
        as.integer(mc_samples)
    )

    if (result$indices[[1]] == -1) {
        return(object)
    }

    # Convert the 'W' matrix by merging columns
    # We add '2' to indices since it's originally 0-based.
    col_groups <- split(seq_len(ncol(object$W)), findInterval(seq_len(ncol(object$W)), c(1, result$indices + 2)))
    # Lambda function to add the columns
    # object$W_agg <- as.matrix(sapply(col_groups, function(cols) rowSums(object$W[, cols, drop = FALSE])))
    object$W_agg <- do.call(cbind, lapply(col_groups, function(cols) rowSums(object$W[, cols, drop = FALSE])))
    rownames(object$W_agg) <- rownames(object$W)
    object$group_agg <- result$indices + 1 # Use R's index system
    object$prob <- as.matrix(result$probabilities)
    dimnames(object$prob) <- list(NULL, colnames(object$X))
    object$iterations <- as.numeric(result$total_iterations)
    object$logLik <- as.numeric(result$log_likelihood)
    object$time <- result$total_time
    object$message <- result$stopping_reason
    object$status <- as.integer(result$finish_id)
    object$cond_prob <- result$q
    object$sd <- result$bootstrap_sol
    dimnames(object$sd) <- dimnames(object$prob)
    object$sd[object$sd == 9999] <- Inf
    object$cond_prob <- aperm(result$q, perm = c(2, 3, 1)) # Correct dimensions
    dimnames(object$cond_prob) <- list(
        NULL,
        colnames(object$X),
        rownames(object$X)
    )
    object$method <- method
    object$ll_threshold <- ll_threshold
    object$param_threshold <- param_threshold
    object$maxtime <- maxtime
    object$maxiter <- maxiter
    object$initial_prob <- initial_prob

    if (method == "mvn_cdf") {
        object$mc_error <- mc_error
        object$mc_samples <- mc_samples
        object$mc_method <- mc_method
    } else if (method == "mcmc") {
        object$step_size <- step_size
        object$samples <- samples
    }

    class(object) <- "eim"
    return(object)
}

#' Performs a matrix-wise Welch's t-test for two eim objects
#'
#' This function compares two `eim` objects (or sets of matrices that can be converted to such objects) by computing a Welch's t-test on each component
#' of their estimated probability matrices (`p`). The Welch test is applied using bootstrap-derived standard deviations, and the result is a matrix
#' of p-values corresponding to each group-candidate combination.
#'
#' It uses Welch's t-test to analyze if there is a significant difference between the estimated probabilities between a treatment and a control set. The test is performed independently for each component of the probability matrix.
#'
#' @inheritParams bootstrap
#'
#' @param object An `eim` object, as returned by [eim].
#' @param object2 A second `eim` object to compare with `object`.
#' @param X A `(b x c)` matrix representing candidate votes per ballot box.
#' @param W A `(b x g)` matrix representing group votes per ballot box.
#' @param X2 A second `(b x c)` matrix to compare with `X`.
#' @param W2 A second `(b x g)` matrix to compare with `W`.
#' @param nboot Integer specifying how many times to run the EM algorithm per object.
#' @param alternative Character string specifying the type of alternative hypothesis to test. Must be one of `"two.sided"` (default), `"greater"`, or `"less"`. If `"two.sided"`, the test checks for any difference in estimated probabilities. If `"greater"`, it tests whether the first object has systematically higher probabilities than the second. If `"less"`, it tests whether the first has systematically lower probabilities.
#' @param ... Additional arguments passed to [bootstrap] and [run_em].
#'
#' @return A list with components:
#'   - `pvals`: a numeric matrix of p-values with the same dimensions as the estimated probability matrices (`pvals`) from the input objects.
#'   - `statistic`: a numeric matrix of t-statistics with the same dimensions as the estimated probability matrices (`pvals`).
#'   - `eim1` and `eim2`: the original `eim` objects used for comparison.
#'
#' Each entry in the pvals matrix is the p-value from Welch's t-test between the corresponding
#' entries of the two estimated probability matrices.
#'
#' @details
#' The user must provide either of the following (but not both):
#' - Two `eim` objects via `object` and `object2`, or
#' - Four matrices: `X`, `W`, `X2`, and `W2`, which will be converted into `eim` objects internally.
#'
#' The Welch test is computed using the formula:
#'
#' \deqn{
#' t_{ij} = \frac{p_{1,ij} - p_{2,ij}}{\sqrt{(s_{1,ij}^2 + s_{2,ij}^2) / n}},
#' }
#' In this expression, \eqn{s_{1,ij}^2} and \eqn{s_{2,ij}^2} represent the bootstrap sample variances for the treatment and control sets, respectively, while \eqn{p_{1,ij}} and \eqn{p_{2,ij}} are the corresponding estimated probability matrices obtained via the EM algorithm. The number of bootstrap samples is denoted by \eqn{n}, and the degrees of freedom for each component are calculated using the Welch-Satterthwaite equation
#'
#' @examples
#' sim <- simulate_election(num_ballots = 100, num_candidates = 3, num_groups = 5, seed = 123)
#' sim2 <- simulate_election(num_ballots = 100, num_candidates = 3, num_groups = 5, seed = 124)
#' eim1 <- eim(sim$X, sim$W)
#' eim2 <- eim(sim2$X, sim2$W)
#'
#' result <- welchtest(object = eim1, object2 = eim2, nboot = 100)
#'
#' # Check which entries are significantly different
#' which(result$pvals < 0.05, arr.ind = TRUE)
#'
#' @export
welchtest <- function(object = NULL,
                      object2 = NULL,
                      X = NULL,
                      W = NULL,
                      X2 = NULL,
                      W2 = NULL,
                      nboot = 50,
                      seed = NULL,
                      alternative = "two.sided",
                      ...) {
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint # It would validate nboot too.

    # Retrieve default values from bootstrap() and update with user parameters
    bootstrap_defaults <- formals(bootstrap)
    bootstrap_args <- modifyList(as.list(bootstrap_defaults), all_params)
    bootstrap_args <- bootstrap_args[names(bootstrap_args) != "..."] # Remove ellipsis

    # Create boolean flags to check which set of inputs is being used
    using_objects <- !is.null(object) && !is.null(object2)
    using_matrices <- !is.null(X) && !is.null(W) && !is.null(X2) && !is.null(W2)

    # Count how many valid input modes are being used
    input_modes <- sum(using_objects, using_matrices)

    if (input_modes == 0) {
        stop("You must provide either (1) two objects, or (2) four matrices (X, X2, W, W2)")
    } else if (input_modes > 1) {
        stop("Please provide only one input mode: either objects or matrices - not multiple.")
    }

    if (using_matrices) {
        object <- eim(X, W)
        object2 <- eim(X2, W2)
    }

    if (ncol(object$X) != ncol(object2$X) || ncol(object$W) != ncol(object2$W)) {
        stop("Column dimensions must be the same for both 'eim' objects")
    }

    # First object
    if (!is.null(all_params$verbose) && all_params$verbose) {
        message("Obtaining the values of the first object.\n")
    }
    boot1 <- do.call(bootstrap, c(
        list(object = object),
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X", "X2", "W", "W2", "json_path")]
    ))
    em1 <- do.call(run_em, c(
        list(object = object),
        bootstrap_args[!names(bootstrap_args) %in% c("verbose", "object", "nboot", "object2", "X", "X2", "W", "W2", "json_path")],
        list(verbose = FALSE)
    ))

    if (!is.null(all_params$verbose) && all_params$verbose) {
        message("Obtaining the values of the second object.\n")
    }
    # Second object
    boot2 <- do.call(bootstrap, c(
        list(object = object2),
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X", "X2", "W", "W2", "json_path")]
    ))
    em2 <- do.call(run_em, c(
        list(object = object2),
        bootstrap_args[!names(bootstrap_args) %in% c("verbose", "nboot", "object", "object2", "X", "X2", "W", "W2", "json_path")],
        list(verbose = FALSE)
    ))

    # Matrix-wise p-values
    var1 <- boot1$sd^2
    var2 <- boot2$sd^2

    se_diff <- sqrt((var1 + var2) / nboot)
    t_stat <- (em1$prob - em2$prob) / se_diff

    df_num <- ((var1 + var2) / nboot)^2
    df_den <- ((var1 / nboot)^2) / (nboot - 1) + ((var2 / nboot)^2) / (nboot - 1)
    df <- df_num / df_den


    pvals <- switch(alternative,
        "two.sided" = 2 * pt(-abs(t_stat), df),
        "greater"   = pt(-t_stat, df),
        "less"      = pt(t_stat, df)
    )
    em1$sd <- boot1$sd
    em2$sd <- boot2$sd

    result <- list()
    result$pvals <- pvals
    result$statistic <- t_stat
    class(em1) <- "eim"
    class(em2) <- "eim"
    result$eim1 <- em1
    result$eim2 <- em2

    return(result)
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
print.eim <- function(x, ...) {
    object <- x
    cat("eim ecological inference model\n")
    # Determine if truncation is needed
    truncated_X <- (nrow(object$X) > 5)
    truncated_W <- (nrow(object$W) > 5)
    is_aggregated <- !is.null(object$W_agg)
    is_method <- !is.null(object$method)
    dim <- if (is_aggregated) "a" else "g"

    cat("Candidates' vote matrix (X) [b x c]:\n")
    print(object$X[1:min(5, nrow(object$X)), ], drop = FALSE) # nolint
    if (truncated_X) cat(".\n.\n.\n") else cat("\n")

    cat("Group-level voter matrix (W) [b x g]:\n")
    print(object$W[1:min(5, nrow(object$W)), , drop = FALSE]) # nolint
    if (truncated_W) cat(".\n.\n.\n") else cat("\n")

    if (is_aggregated) {
        cat("Macro group-level voter matrix (W_agg) [b x a]:\n")
        print(object$W_agg[1:min(5, nrow(object$W_agg)), , drop = FALSE]) # nolint
        truncated_W_agg <- (nrow(object$W_agg) > 5)
        if (truncated_W_agg) cat(".\n.\n.\n") else cat("\n")
    }

    if (is_method) {
        cat(sprintf("Estimated probability [%s x c]:\n", dim))
        truncated_P <- (nrow(object$prob) > 5)
        print(round(object$prob[1:min(5, nrow(object$prob)), ], 3)) # nolint
        if (truncated_P) cat(".\n.\n.\n") else cat("\n")
    }
    # Consider showing matrices first
    if (!is.null(object$sd)) {
        cat(sprintf("Standard deviation of the estimated probability matrix [%s x c]:\n", dim))
        truncated_boot <- (nrow(object$sd) > 5)
        print(round(object$sd[1:min(5, nrow(object$sd)), ], 3)) # nolint
        if (truncated_boot) cat(".\n.\n.\n") else cat("\n")
    }
    if (is_method) {
        cat("Method:\t", object$method, "\n")
        if (object$method == "mcmc") {
            cat("Step size (M):", object$step_size, "\n")
            cat("Samples (S):", object$samples, "\n")
        } else if (object$method == "mvn_cdf") {
            cat("Montecarlo method:", object$mc_method, "\n")
            cat("Montecarlo iterations:", object$mc_samples, "\n")
            cat("Montecarlo error:", object$mc_error, "\n")
        }
        cat("Total Iterations:", object$iterations, "\n")
        cat("Total Time (s):", object$time, "\n")
        cat("Log-likelihood:", tail(object$logLik, 1), "\n")
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
        ballots = nrow(object$X)
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
as.matrix.eim <- function(x, ...) {
    object <- x
    if (is.null(object$prob)) {
        stop(paste0(
            "Probability matrix not available. Run 'run_em()'."
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
#' \dontrun{
#' model <- eim(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#'
#' model <- run_em(model)
#'
#' # Save as RDS
#' save_eim(model, "model_results.rds")
#'
#' # Save as JSON
#' save_eim(model, "model_results.json")
#'
#' # Save as CSV
#' save_eim(model, "model_results.csv")
#' }
#'
#' @name save_eim
#' @aliases save_eim()
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

        jsonlite::write_json(json_data, filename, pretty = TRUE, auto_unbox = TRUE, digits = 10)
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
# write.csv.eim <- function(object, filename, ...) {
#    if (!inherits(object, "eim")) {
#        stop("The object must be initialized with the `eim()` function.")
#    }
#    if (!is.character(filename) || length(filename) != 1) {
#        stop("Invalid filename. Please provide a valid file path as a character string.")
#    }
#
#    # Get file extension
#    file_ext <- tools::file_ext(filename)
#
#    if (file_ext != "csv") {
#        stop("The filepath provided must end with '.csv'")
#    }
#
#    if (!is.null(object$prob)) {
#        write.csv(as.matrix(object$prob), filename, row.names = TRUE)
#        message("Probability matrix saved as CSV: ", filename)
#    } else {
#        stop("The `run_em()` method must be called for saving a '.csv' file.")
#    }
# }

#' @noRd
#' @export
# dput.eim <- function(object, filename, ...) {
#    if (!inherits(object, "eim")) {
#        stop("The object must be initialized with the `eim()` function.")
#    }
#    if (!is.character(filename) || length(filename) != 1) {
#        stop("Invalid filename. Please provide a valid file path as a character string.")
#    }
#
#    # Get file extension
#    file_ext <- tools::file_ext(filename)
#
#    if (file_ext != "rda") {
#        stop("The filepath provided must end with '.rda'")
#    }
#    saveRDS(object, file = filename)
#    message("Results saved as RDS: ", filename)
# }


#' @title Extract log-likelihood
#' @description
#'   Return the log-likelihood of the last EM iteration
#'
#' @param object An `eim` object
#' @param ... Additional parameters that will be ignored
#'
#' @return A numeric value with the log-likelihood from the last iteration.
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
