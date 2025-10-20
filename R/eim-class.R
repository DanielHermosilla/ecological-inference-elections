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
            "(1)\tA json path\n",
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
            "mcmc_samples",
            "mcmc_stepsize",
            "mvncdf_method",
            "mvncdf_samples",
            "mvncdf_error",
            "miniter",
            "adjust_prob_cond_method",
            "adjust_prob_cond_every"
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
#' @param allow_mismatch Boolean, if `TRUE`, allows a mismatch between the voters and votes for each ballot-box. If `FALSE`, throws an error if there is a mismatch. By default it is `TRUE`. See **Notes** for more details.
#'
#' @param maxiter An optional integer indicating the maximum number of EM iterations.
#'   The default value is `1000`.
#'
#' @param miniter An optional integer indicating the minimum number of EM iterations. The default value is `0`.
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
#' @param compute_ll An optional boolean indicating whether to compute the log-likelihood at each iteration. The default value is `TRUE`.
#'
#' @param adjust_prob_cond_method An optional string indicating the method to adjust the conditional probability so that for each candidate, the sum product of voters and conditional probabilities across groups equals the votes obtained by the candidate. It can take values: `""` if no adjusting is made, `"lp"` if the adjustment is based on a linear programming that penalizes with L1-norm, `"project_lp"` if the adjustment is performed using projection and linear programming (this is the default)
#'
#' @param adjust_prob_cond_every An optional boolean indicating whether to adjust the conditional probability on every iteration (if `TRUE`), or only at the conditional probabilities obtained at the end of the EM algorithm (if `FALSE`, this is the default). This parameter applies only if `adjust_prob_conditional_method` is `lp` or `project_lp`.
#'
#' @param verbose An optional boolean indicating whether to print informational messages during the EM
#'   iterations. The default value is `FALSE`.
#'
#' @param seed An optional integer indicating the random seed for the randomized algorithms. This argument is only applicable if `initial_prob = "random"` or `method` is either `"mcmc"` or `"mvn_cdf"`.
#'
#' @param group_agg An optional vector of increasing integers from 1 to the number of columns in `W`, specifying how to aggregate groups in `W` before running the EM algorithm. Each value represents the highest column index included in each aggregated group. For example, if `W` has four columns, `group_agg = c(2, 4)` indicates that columns 1 and 2 should be combined into one group, and columns 3 and 4 into another. Defaults to `NULL`, in which case no group aggregation is performed.
#'
#' @param mcmc_stepsize An optional integer specifying the step size for the `mcmc`
#'   algorithm. This parameter is only applicable when `method = "mcmc"` and will
#'   be ignored otherwise. The default value is `3000`.
#'
#' @param mcmc_samples An optional integer indicating the number of samples to generate for the
#'   **MCMC** method. This parameter is only relevant when `method = "mcmc"`.
#'   The default value is `1000`.
#'
#' @param mvncdf_method An optional string specifying the method used to estimate the `mvn_cdf` method
#'   via a Monte Carlo simulation. Accepted values are `genz` and `genz2`, with `genz`
#'   set as the default. This parameter is only applicable when `method = "mvn_cdf"`. See **References** for more details.
#'
#' @param mvncdf_error An optional numeric value defining the error threshold for the Monte Carlo
#'   simulation when estimating the `mvn_cdf` method. The default value is `1e-6`. This parameter is only relevant
#' when `method = "mvn_cdf"`.
#'
#' @param mvncdf_samples An optional integer specifying the number of Monte Carlo
#'   samples for the `mvn_cdf` method. The default value is `5000`. This argument is only applicable when `method = "mvn_cdf"`.
#'
#' @param scale_factor An optional numeric value used to scale down the `X` and `W` matrices before executing the EM algorithm. This scaling can help improve performance when dealing with large vote counts. For example if `scale_factor = 2` all elements of `X` and `W` are divided by two and rounded. The default value is `1`, which means no scaling is applied. In case the scaling results in mismatch between `W` and `X`, ensure that `allow_mismatch = TRUE`.
#'
#' @param ... Added for compability
#'
#' @references
#' [Thraves, C., Ubilla, P. and Hermosilla, D.: *"Fast Ecological Inference Algorithm for the RxC Case"*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834). Aditionally, the MVN CDF is computed by the methods introduced in [Genz, A. (2000). Numerical computation of multivariate normal probabilities. *Journal of Computational and Graphical Statistics*](https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities)
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
#' If there are ballot-boxes with mismatch between `W` and `X`, and `allow_mismatch = TRUE`, then: if `method = "exact"`, at each ballot-box with mismatch D'Hont is applied to add or remove the necessary voters from (`W`) so that its total match the total number of votes (`X`); if method is "`mvn_pdf`", "`mvn_cdf`" or "`mcmc`", the number of voters (`W`) of the ballot-box with mismatch is scaled to match its total number of votes (`X`).
#'
#' @seealso The [eim] object implementation.
#'
#' @return
#' The function returns an `eim` object with the function arguments and the following attributes:
#' \describe{
#'   \item{prob}{The estimated probability matrix `(g x c)`.}
#' 	 \item{cond_prob}{A `(g x c x b)` 3d-array with the probability that a at each ballot-box a voter of each group voted for each candidate, given the observed outcome at the particular ballot-box.}
#' 	 \item{expected_outcome}{A `(g x c x b)` 3d-array with the expected votes cast for each ballot box.}
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
#' Aditionally, it will create `mcmc_samples` and `mcmc_stepsize` parameters if the specified `method = "mcmc"`, or `mvncdf_method`, `mvncdf_error` and `mvncdf_samples` if `method = "mvn_cdf"`.
#'
#' Also, if the eim object supplied is created with the function [simulate_election], it also returns the real probability and unobserved votes with the name `real_prob` and `outcome` respectively. See [simulate_election].
#'
#' If `group_agg` is different than `NULL`, two values are returned: `W_agg` a `(b x a)` matrix with the number of voters of each aggregated group o each ballot-box, and `group_agg` the same input vector.
#'
#' @examples
#' \donttest{
#' # Example 1: Compute the Expected-Maximization with default settings
#' simulations <- simulate_election(
#'     num_ballots = 300,
#'     num_candidates = 5,
#'     num_groups = 3,
#' )
#' model <- eim(simulations$X, simulations$W)
#' model <- run_em(model) # Returns the object with updated attributes
#'
#' # Example 2: Compute the Expected-Maximization using the mvn_pdf method
#' model <- run_em(
#'     object = model,
#'     method = "mvn_pdf",
#' )
#'
#' # Example 3: Run the mvn_cdf method with default settings
#' model <- run_em(object = model, method = "mvn_cdf")
#' }
#' \dontrun{
#' # Example 4: Perform an Exact estimation using user-defined parameters
#'
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
                   miniter = 0,
                   maxtime = 3600,
                   param_threshold = 0.001,
                   ll_threshold = as.double(-Inf),
                   compute_ll = TRUE,
                   seed = NULL,
                   verbose = FALSE,
                   group_agg = NULL,
                   mcmc_samples = 1000,
                   mcmc_stepsize = 3000,
                   mvncdf_method = "genz",
                   mvncdf_error = 1e-3,
                   mvncdf_samples = 5000,
                   adjust_prob_cond_method = "project_lp",
                   adjust_prob_cond_every = FALSE,
                   scale_factor = 1,
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


    # Applies a scaling
    if (scale_factor != 1) {
        object$X <- round(object$X / all_params$scale_factor)
        object$W <- round(object$W / all_params$scale_factor)
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
        W <- .dhondt_correction(object$W, object$X)
        message("Applying a D'Hondt correction for correcting mismatches in W")
        # stop("run_em: Exact method isn't supported with mismatch")
    }

    # Handle the group aggregation, if provided
    if (!is.null(group_agg)) {
        sizes <- diff(c(0, group_agg))
        rep_labels <- rep(seq_along(sizes), sizes)
        groups <- split(seq_len(ncol(object$W)), rep_labels)
        Wagg <- do.call(
            cbind,
            lapply(groups, function(cols) rowSums(object$W[, cols, drop = FALSE]))
        )
        rownames(Wagg) <- rownames(object$W)
        object$W_agg <- Wagg
        object$group_agg <- group_agg
    }

    object$method <- method

    # Default values
    if (method == "mcmc") {
        # Step size
        object$mcmc_stepsize <- as.integer(if ("mcmc_stepsize" %in% names(all_params)) all_params$mcmc_stepsize else 3000)
        # Samples
        object$mcmc_samples <- as.integer(if ("mcmc_samples" %in% names(all_params)) all_params$mcmc_samples else 1000)
        # Burn in
        object$burn_in <- as.integer(if ("burn_in" %in% names(all_params)) all_params$burn_in else 10000)
    } else if (method == "mvn_cdf") {
        # Montecarlo method
        object$mvncdf_method <- if ("mvncdf_method" %in% names(all_params)) all_params$mvncdf_method else "genz"
        # Montecarlo samples
        object$mvncdf_samples <- if ("mvncdf_samples" %in% names(all_params)) all_params$mvncdf_samples else 5000
        # Montecarlo error
        object$mvncdf_error <- if ("mvncdf_error" %in% names(all_params)) all_params$mvncdf_error else 1e-3
    }

    W <- if (is.null(object$W_agg)) object$W else object$W_agg
    # RsetParameters(t(object$X), W)

    resulting_values <- EMAlgorithmFull(
        t(object$X),
        W,
        method,
        if (is.character(initial_prob)) initial_prob else "custom",
        maxiter,
        maxtime,
        param_threshold,
        ll_threshold,
        compute_ll,
        verbose,
        as.integer(if (!is.null(object$mcmc_stepsize)) object$mcmc_stepsize else 3000),
        as.integer(if (!is.null(object$mcmc_samples)) object$mcmc_samples else 1000),
        if (!is.null(object$mvncdf_method)) object$mvncdf_method else "genz",
        as.numeric(if (!is.null(object$mvncdf_error)) object$mvncdf_error else 1e-3),
        as.numeric(if (!is.null(object$mvncdf_samples)) object$mvncdf_samples else 5000),
        miniter,
        adjust_prob_cond_method,
        adjust_prob_cond_every,
        if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1)
    )
    # ---------- ... ---------- #

    object$cond_prob <- resulting_values$q
    # object$cond_prob <- aperm(resulting_values$q, perm = c(2, 3, 1)) # Correct dimensions
    dimnames(object$cond_prob) <- list(
        colnames(W),
        colnames(object$X),
        rownames(object$X)
    )
    object$expected_outcome <- resulting_values$expected_outcome
    dimnames(object$expected_outcome) <- list(
        colnames(W),
        colnames(object$X),
        rownames(object$X)
    )
    object$prob <- as.matrix(resulting_values$result)
    dimnames(object$prob) <- list(colnames(W), colnames(object$X))
    object$iterations <- as.numeric(resulting_values$total_iterations)
    if (compute_ll) {
        object$logLik <- as.numeric(resulting_values$log_likelihood[length(resulting_values$log_likelihood)])
    }
    object$time <- resulting_values$total_time
    object$message <- resulting_values$stopping_reason
    object$status <- as.integer(resulting_values$finish_id)
    # Add function arguments
    object$miniter <- miniter
    object$maxiter <- maxiter
    object$maxtime <- maxtime
    object$param_threshold <- param_threshold
    object$ll_threshold <- ll_threshold
    object$initial_prob <- initial_prob
    object$adjust_prob_cond_method <- adjust_prob_cond_method
    object$adjust_prob_cond_every <- adjust_prob_cond_every

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
#' \donttest{
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
#'     method = "mult",
#'     maxiter = 500,
#'     verbose = FALSE,
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
#'     nboot = 15,
#'     method = "mvn_pdf",
#'     maxiter = 100,
#'     maxtime = 5,
#'     param_threshold = 0.01,
#'     allow_mismatch = FALSE
#' )
#'
#' print(model$sd)
#' }
#'
#' # Example 3: Using a JSON file as input
#'
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
                      nboot = 100,
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

    # Handle the group aggregation, if provided
    if (!is.null(all_params$group_agg)) {
        sizes <- diff(c(0, all_params$group_agg))
        rep_labels <- rep(seq_along(sizes), sizes)
        groups <- split(seq_len(ncol(object$W)), rep_labels)
        Wagg <- do.call(
            cbind,
            lapply(groups, function(cols) rowSums(object$W[, cols, drop = FALSE]))
        )
        rownames(Wagg) <- rownames(object$W)
        object$W_agg <- Wagg
        object$group_agg <- all_params$group_agg
    }

    # I need to define the method before on this case
    method <- if (!is.null(all_params$method)) all_params$method else "mult"

    # Applies a scaling
    if (!is.null(all_params$scale_factor) && all_params$scale_factor != 1) {
        object$X <- round(object$X / all_params$scale_factor)
        object$W <- round(object$W / all_params$scale_factor)
    }

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
        if (method == "exact") {
            W <- .dhondt_correction(object$W, object$X)
        }
        message("Applying a D'Hondt correction for correcting mismatches in W")

        # stop("run_em: Exact method isn't supported with mismatch")
    }
    # Set seed for reproducibility
    if (!is.null(seed)) set.seed(seed)

    # Extract parameters with defaults if missing
    W <- if (is.null(object$W_agg)) object$W else object$W_agg
    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    miniter <- if (!is.null(all_params$miniter)) all_params$miniter else 0
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.001
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE
    compute_ll <- if (!is.null(all_params$compute_ll)) all_params$compute_ll else TRUE
    adjust_prob_cond_method <- if (!is.null(all_params$adjust_prob_cond_method)) all_params$adjust_prob_cond_method else ""
    adjust_prob_cond_every <- if (!is.null(all_params$adjust_prob_cond_every)) all_params$adjust_prob_cond_every else FALSE

    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    mcmc_stepsize <- 0L
    mcmc_samples <- 0L
    mvncdf_method <- ""
    mvncdf_samples <- 0L
    mvncdf_error <- 0.0

    if (method == "mcmc") {
        mcmc_stepsize <- if (!is.null(all_params$mcmc_stepsize)) all_params$mcmc_stepsize else 3000
        mcmc_samples <- if (!is.null(all_params$mcmc_samples)) all_params$mcmc_samples else 1000
        burn_in <- if (!is.null(all_params$burn_in)) all_params$burn_in else 10000
    } else if (method == "mvn_cdf") {
        mvncdf_method <- if (!is.null(all_params$mvncdf_method)) all_params$method else "genz"
        mvncdf_samples <- if (!is.null(all_params$mvncdf_samples)) all_params$mvncdf_samples else 5000
        mvncdf_error <- if (!is.null(all_params$mvncdf_error)) all_params$mvncdf_error else 1e-6
    } # Call C bootstrap function

    result <- bootstrapAlg(
        t(object$X),
        W,
        as.integer(nboot),
        as.character(method),
        if (is.character(initial_prob)) as.character(initial_prob) else "custom",
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(compute_ll),
        as.logical(verbose),
        as.integer(mcmc_stepsize),
        as.integer(mcmc_samples),
        as.character(mvncdf_method),
        as.double(mvncdf_error),
        as.integer(mvncdf_samples),
        as.integer(miniter),
        as.character(adjust_prob_cond_method),
        as.logical(adjust_prob_cond_every),
        if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1)
    )

    object$sd <- result
    dimnames(object$sd) <- list(colnames(W), colnames(object$X))
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
#' @param allow_mismatch Boolean, if `TRUE`, allows a mismatch between the voters and votes for each ballot-box. If `FALSE`, throws an error if there is a mismatch. By default it is `TRUE`. See **Notes** in [run_em] for more details.
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
#' # a second that has the original groups 3, 4, 5, and 6.
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
#'     sd_statistic = "average", seed = 42
#' )
#'
#' solution$group_agg # c(3, 4, 5, 6, 8)
#' # This means that the resulting group aggregation consists of
#' # five macro-groups: one that includes the original groups 1, 2, and 3;
#' # three singleton groups (4, 5, and 6); and one macro-group that includes groups 7 and 8.
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
                          nboot = 100,
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
    if (!is.null(all_params$scale_factor) && all_params$scale_factor != 1) {
        object$X <- round(object$X / all_params$scale_factor)
        object$W <- round(object$W / all_params$scale_factor)
    }

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
        if (method == "exact") {
            W <- .dhondt_correction(object$W, object$X)
        }
        message("Applying a D'Hondt correction for correcting mismatches in W")

        # stop("run_em: Exact method isn't supported with mismatch")
    }

    # Extract parameters with defaults if missing
    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    miniter <- if (!is.null(all_params$miniter)) all_params$miniter else 0
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.01
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE
    compute_ll <- if (!is.null(all_params$compute_ll)) all_params$compute_ll else TRUE
    adjust_prob_cond_method <- if (!is.null(all_params$adjust_prob_cond_method)) all_params$adjust_prob_cond_method else ""
    adjust_prob_cond_every <- if (!is.null(all_params$adjust_prob_cond_every)) all_params$adjust_prob_cond_every else FALSE


    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    mcmc_stepsize <- 0L
    mcmc_samples <- 0L
    mvncdf_method <- ""
    mvncdf_samples <- 0L
    mvncdf_error <- 0.0

    if (method == "mcmc") {
        mcmc_stepsize <- if (!is.null(all_params$mcmc_stepsize)) all_params$mcmc_stepsize else 3000
        mcmc_samples <- if (!is.null(all_params$mcmc_samples)) all_params$mcmc_samples else 1000
        burn_in <- if (!is.null(all_params$burn_in)) all_params$burn_in else 10000
    } else if (method == "mvn_cdf") {
        mvncdf_method <- if (!is.null(all_params$mvncdf_method)) all_params$method else "genz"
        mvncdf_samples <- if (!is.null(all_params$mvncdf_samples)) all_params$mvncdf_samples else 5000
        mvncdf_error <- if (!is.null(all_params$mvncdf_error)) all_params$mvncdf_error else 1e-6
    }

    result <- groupAgg(
        as.character(sd_statistic),
        as.double(sd_threshold),
        as.logical(feasible),
        t(object$X),
        object$W,
        as.integer(nboot),
        as.character(method),
        if (is.character(initial_prob)) as.character(initial_prob) else "custom",
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(compute_ll),
        as.logical(verbose),
        as.integer(mcmc_stepsize),
        as.integer(mcmc_samples),
        as.character(mvncdf_method),
        as.double(mvncdf_error),
        as.integer(mvncdf_samples),
        as.integer(miniter),
        as.character(adjust_prob_cond_method),
        as.logical(adjust_prob_cond_every),
        if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1)
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
#' @param allow_mismatch Boolean, if `TRUE`, allows a mismatch between the voters and votes for each ballot-box. If `FALSE`, throws an error if there is a mismatch. By default it is `TRUE`. See **Notes** in [run_em] for more details.
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
#'
#' \donttest{
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
                        nboot = 100,
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

    if (!is.null(seed)) set.seed(seed) # Set seed for reproducibility
    # Initialize eim object if needed
    if (is.null(object)) {
        object <- eim(X, W, json_path)
    } else if (!inherits(object, "eim")) {
        stop("get_agg_opt: The object must be initialized with the `eim()` function.")
    }

    # Note: Mismatch restricted methods are checked inside .validate_compute
    # Method needs to be defined before

    if (!is.null(all_params$scale_factor) && all_params$scale_factor != 1) {
        object$X <- round(object$X / all_params$scale_factor)
        object$W <- round(object$W / all_params$scale_factor)
    }

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
        if (method == "exact") {
            W <- .dhondt_correction(object$W, object$X)
        }
        message("Applying a D'Hondt correction for correcting mismatches in W")

        # stop("run_em: Exact method isn't supported with mismatch")
    }

    initial_prob <- if (!is.null(all_params$initial_prob)) all_params$initial_prob else "group_proportional"
    maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
    miniter <- if (!is.null(all_params$miniter)) all_params$miniter else 0
    maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
    param_threshold <- if (!is.null(all_params$param_threshold)) all_params$param_threshold else 0.01
    verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE
    compute_ll <- if (!is.null(all_params$compute_ll)) all_params$compute_ll else TRUE
    adjust_prob_cond_method <- if (!is.null(all_params$adjust_prob_cond_method)) all_params$adjust_prob_cond_method else ""
    adjust_prob_cond_every <- if (!is.null(all_params$adjust_prob_cond_every)) all_params$adjust_prob_cond_every else FALSE

    # R does a subtle type conversion when handing -Inf. Hence, we'll use a direct assignment
    if ("ll_threshold" %in% names(all_params)) {
        ll_threshold <- all_params$ll_threshold
    } else {
        ll_threshold <- as.double(-Inf)
    }

    # Handle method-specific defaults
    mcmc_stepsize <- 0L
    mcmc_samples <- 0L
    mvncdf_method <- ""
    mvncdf_samples <- 0L
    mvncdf_error <- 0.0

    if (method == "mcmc") {
        mcmc_stepsize <- if (!is.null(all_params$mcmc_stepsize)) all_params$mcmc_stepsize else 3000
        mcmc_samples <- if (!is.null(all_params$mcmc_samples)) all_params$mcmc_samples else 1000
        burn_in <- if (!is.null(all_params$burn_in)) all_params$burn_in else 10000
    } else if (method == "mvn_cdf") {
        mvncdf_method <- if (!is.null(all_params$mvncdf_method)) all_params$mvncdf_method else "genz"
        mvncdf_samples <- if (!is.null(all_params$mvncdf_samples)) all_params$mvncdf_samples else 5000
        mvncdf_error <- if (!is.null(all_params$mvncdf_error)) all_params$mvncdf_error else 1e-6
    }

    result <- groupAggGreedy(
        as.character(sd_statistic),
        as.double(sd_threshold),
        t(object$X),
        object$W,
        as.integer(nboot),
        as.character(method),
        if (is.character(initial_prob)) as.character(initial_prob) else "group_proportional",
        as.integer(maxiter),
        as.double(maxtime),
        as.double(param_threshold),
        as.double(ll_threshold),
        as.logical(compute_ll),
        as.logical(verbose),
        as.integer(mcmc_stepsize),
        as.integer(mcmc_samples),
        as.character(mvncdf_method),
        as.double(mvncdf_error),
        as.integer(mvncdf_samples),
        as.integer(miniter),
        as.character(adjust_prob_cond_method),
        as.logical(adjust_prob_cond_every)
        # if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1)
    )

    if (result$indices[[1]] == -1) {
        return(object)
    }

    # Convert the 'W' matrix by merging columns
    # We add '2' to indices since it's originally 0-based.
    col_groups <- split(seq_len(ncol(object$W)), findInterval(seq_len(ncol(object$W)), c(1, result$indices + 2)))
    # Lambda function to add the columns
    # object$W_agg <- as.matrix(sapply(col_groups, function(cols) rowSums(object$W[, cols, drop = FALSE])))
    object$cond_prob <- result$q
    # object$cond_prob <- aperm(result$q, perm = c(2, 3, 1)) # Correct dimensions
    dimnames(object$cond_prob) <- list(
        NULL,
        colnames(object$X),
        rownames(object$X)
    )
    object$expected_outcome <- result$expected_outcome
    dimnames(object$expected_outcome) <- list(
        NULL,
        colnames(object$X),
        rownames(object$X)
    )
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
    object$sd <- result$bootstrap_sol
    dimnames(object$sd) <- dimnames(object$prob)
    object$sd[object$sd == 9999] <- Inf
    object$method <- method
    object$ll_threshold <- ll_threshold
    object$param_threshold <- param_threshold
    object$maxtime <- maxtime
    object$maxiter <- maxiter
    object$miniter <- miniter
    object$initial_prob <- initial_prob
    object$adjust_prob_cond_method <- adjust_prob_cond_method
    object$adjust_prob_cond_every <- adjust_prob_cond_every

    if (method == "mvn_cdf") {
        object$mvncdf_error <- mvncdf_error
        object$mvncdf_samples <- mvncdf_samples
        object$mvncdf_method <- mvncdf_method
    } else if (method == "mcmc") {
        object$mcmc_stepsize <- mcmc_stepsize
        object$mcmc_samples <- mcmc_samples
    }

    class(object) <- "eim"
    return(object)
}

#' Performs a matrix-wise Wald test for two eim objects
#'
#' This function compares two `eim` objects (or sets of matrices that can be converted to such objects) by computing a Wald test on each component
#' of their estimated probability matrices. The Wald test is applied using bootstrap-derived standard deviations, and the result is a matrix
#' of p-values corresponding to each group-candidate combination.
#'
#' It uses Wald test to analyze if there is a significant difference between the estimated probabilities between a treatment and a control set. The test is performed independently for each component of the probability matrix.
#'
#' @inheritParams bootstrap
#'
#' @param object1 An `eim` object, as returned by [eim].
#' @param object2 A second `eim` object to compare with `object`.
#' @param X1 A `(b x c)` matrix representing candidate votes per ballot box.
#' @param W1 A `(b x g)` matrix representing group votes per ballot box.
#' @param X2 A second `(b x c)` matrix to compare with `X`.
#' @param W2 A second `(b x g)` matrix to compare with `W`.
#' @param nboot Integer specifying how many times to run the EM algorithm per object.
#' @param alternative Character string specifying the type of alternative hypothesis to test. Must be one of `"two.sided"` (default), `"greater"`, or `"less"`. If `"two.sided"`, the test checks for any difference in estimated probabilities. If `"greater"`, it tests whether the first object has systematically higher probabilities than the second. If `"less"`, it tests whether the first has systematically lower probabilities.
#' @param ... Additional arguments passed to [bootstrap] and [run_em].
#'
#' @return A list with components:
#'   - `pvals`: a numeric matrix of p-values with the same dimensions as the estimated probability matrices (`pvals`) from the input objects.
#'   - `statistic`: a numeric matrix of z-statistics with the same dimensions as the estimated probability matrices (`pvals`).
#'   - `eim1` and `eim2`: the original `eim` objects used for comparison.
#'
#' Each entry in the pvals matrix is the p-value from Wald test between the corresponding
#' entries of the two estimated probability matrices.
#'
#' @details
#' The user must provide either of the following (but not both):
#' - Two `eim` objects via `object1` and `object2`, or
#' - Four matrices: `X1`, `W1`, `X2`, and `W2`, which will be converted into `eim` objects internally.
#'
#' The Wald test is computed using the formula:
#'
#' \deqn{
#' z_{ij} = \frac{p_{1,ij}-p_{2,ij}}{\sqrt{s_{1,ij}^2+s_{2,ij}^2}}
#' }
#' In this expression, \eqn{s_{1,ij}^2} and \eqn{s_{2,ij}^2} represent the bootstrap sample variances for the treatment and control sets, respectively, while \eqn{p_{1,ij}} and \eqn{p_{2,ij}} are the corresponding estimated probability matrices obtained via the EM algorithm.
#'
#' @examples
#' sim1 <- simulate_election(num_ballots = 100, num_candidates = 3, num_groups = 5, seed = 123)
#' sim2 <- simulate_election(num_ballots = 100, num_candidates = 3, num_groups = 5, seed = 124)
#'
#' result <- waldtest(sim1, sim2, nboot = 100)
#'
#' # Check which entries are significantly different
#' which(result$pvals < 0.05, arr.ind = TRUE)
#'
#' @export
waldtest <- function(object1 = NULL,
                     object2 = NULL,
                     X1 = NULL,
                     W1 = NULL,
                     X2 = NULL,
                     W2 = NULL,
                     nboot = 100,
                     seed = NULL,
                     alternative = "two.sided",
                     ...) {
    object <- object1
    X <- X1
    W <- W1
    provided <- c(!is.null(X1), !is.null(W1), !is.null(X2), !is.null(W2))
    invalidMat <- any(provided) && !all(provided)
    using_objects <- !is.null(object1) && !is.null(object2)
    using_matrices <- all(provided)
    input_modes <- sum(using_objects, using_matrices)

    if (input_modes == 0) {
        stop("Invalid input: must supply two objects or four matrices.")
    } else if (input_modes > 1 || invalidMat) {
        stop("Invalid input: you must provide either: two eim objects (object1, object2), or four matrices (X1, X2, W1, W2), but not both.")
    }

    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    # .validate_compute(all_params) # nolint # It would validate nboot too.

    # Retrieve default values from bootstrap() and update with user parameters
    bootstrap_defaults <- formals(bootstrap)
    bootstrap_args <- modifyList(as.list(bootstrap_defaults), all_params)
    bootstrap_args <- bootstrap_args[names(bootstrap_args) != "..."] # Remove ellipsis

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
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X", "X1", "X2", "W", "W1", "W2", "json_path")],
        list(verbose = FALSE)
    ))
    em1 <- do.call(run_em, c(
        list(object = object),
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X1", "X", "X2", "W1", "W", "W2", "json_path")],
        list(verbose = FALSE)
    ))

    if (!is.null(all_params$verbose) && all_params$verbose) {
        message("Obtaining the values of the second object.\n")
    }
    # Second object
    boot2 <- do.call(bootstrap, c(
        list(object = object2),
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X", "X1", "X2", "W", "W1", "W2", "json_path")],
        list(verbose = FALSE)
    ))
    em2 <- do.call(run_em, c(
        list(object = object2),
        bootstrap_args[!names(bootstrap_args) %in% c("object", "object2", "X", "X1", "X2", "W", "W1", "W2", "json_path")],
        list(verbose = FALSE)
    ))

    # Matrix-wise p-values
    var1 <- boot1$sd^2
    var2 <- boot2$sd^2

    delta <- em1$prob - em2$prob
    se_delta <- sqrt(var1 + var2)

    z <- delta / se_delta

    pvals <- switch(alternative,
        "two.sided" = 2 * pnorm(-abs(z)),
        "greater"   = pnorm(-z),
        "less"      = pnorm(z)
    )

    em1$sd <- boot1$sd
    em2$sd <- boot2$sd

    result <- list()
    result$pvals <- pvals
    result$statistic <- z
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
            cat("Step size (M):", object$mcmc_stepsize, "\n")
            cat("Samples (S):", object$mcmc_samples, "\n")
        } else if (object$method == "mvn_cdf") {
            cat("Montecarlo method:", object$mvncdf_method, "\n")
            cat("Montecarlo iterations:", object$mvncdf_samples, "\n")
            cat("Montecarlo error:", object$mvncdf_error, "\n")
        }
        cat("Total Iterations:", object$iterations, "\n")
        cat("Total Time (s):", object$time, "\n")
        if (!is.null(object$logLik)) {
            cat("Log-likelihood:", tail(object$logLik, 1), "\n")
        }
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
#' \donttest{
#' model <- eim(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))
#'
#' model <- run_em(model)
#'
#' td <- tempdir()
#' out_rds <- file.path(td, "model_results.rds")
#' out_json <- file.path(td, "model_results.json")
#' out_csv <- file.path(td, "model_results.csv")
#'
#' # Save as RDS
#' save_eim(model, filename = out_rds)
#'
#' # Save as JSON
#' save_eim(model, filename = out_json)
#'
#' # Save as CSV
#' save_eim(model, filename = out_csv)
#'
#' # Remove the files
#' files <- c(out_rds, out_json, out_csv)
#' file.remove(files)
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
            val <- object[[name]]

            # if it's our 3-D array cond_prob, swap dim 1 ↔ dim 3
            if (identical(name, "cond_prob") &&
                is.array(val) &&
                length(dim(val)) == 3) {
                val <- aperm(val, perm = c(3, 1, 2))
            }

            json_data[[name]] <- val
        }
        # Add the names of ballot boxes
        if (!is.null(object$X) && !is.null(rownames(object$X))) {
            json_data$ballotbox_id <- rownames(object$X)
        }
        # Add the names of candidates
        if (!is.null(object$X) && !is.null(colnames(object$X))) {
            json_data$candidates_id <- colnames(object$X)
        }
        # Add the names of ballot boxes
        if (!is.null(object$W) && !is.null(colnames(object$W))) {
            json_data$group_id <- colnames(object$W)
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

#' @export
exactLL <- function(object, scale_factor = 1) {
    if (is.null(object$X) || is.null(object$W) || is.null(object$prob)) {
        stop("The object must contain X, W and prob matrices.")
    }

    # Applies a scaling
    if (scale_factor != 1) {
        object$X <- round(object$X / scale_factor)
        object$W <- round(object$W / scale_factor)
        object$W <- .dhondt_correction(object$W, object$X)
    }

    ll <- computeExactLL(t(object$X), object$W, object$prob)
    ll
}
