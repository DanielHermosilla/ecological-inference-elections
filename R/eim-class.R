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
#' @param V Optional `(b x a)` matrix with the attributes for each ballot box. This is only used for parametric models.
#'
#' @param json_path A path to a JSON file containing `X`, `W` (and optionally `V`) fields, stored as nested arrays. It may contain additional fields with other attributes, which will be added to the returned object.
#'
#' @details
#' If `X` and `W` are directly supplied, they must match the
#' dimensions of ballot boxes `(b)`. Alternatively, if `json_path` is provided, the function expects
#' the JSON file to contain elements named `"X"` and `"W"` (and optionally `"V"`) under the
#' top-level object. This two approaches are **mutually exclusable**, yielding an error otherwise.
#'
#' When `V` is supplied, the object is treated as parametric and includes the `V` attribute.
#'
#' Internally, this function also initializes the corresponding instance within
#' the low-level (C-based) API, ensuring the data is correctly registered for
#' further processing by the EM algorithm.
#'
#' @return A list of class `eim` containing:
#' \describe{
#'   \item{\code{X}}{The candidate votes matrix \code{(b x c)}.}
#'   \item{\code{W}}{The group votes matrix \code{(b x g)}.}
#'   \item{\code{V}}{The parametric covariates matrix \code{(b x a)}, when provided.}
#' }
#'
#' @note
#' A way to generate synthetic data for `X` and `W` is by using the [simulate_election] function. See Example 2 below.
#' This constructor can be used for both non-parametric and parametric models (by providing `V`).
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
#' 	 \item \code{\link{plot.eim}} - Plots the probability matrix.
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
eim <- function(X = NULL, W = NULL, V = NULL, json_path = NULL) {
    x_provided <- !is.null(X)
    w_provided <- !is.null(W)
    v_provided <- !is.null(V)
    xw_provided <- x_provided || w_provided
    xwv_provided <- xw_provided || v_provided
    json_provided <- !is.null(json_path)

    if (sum(x_provided, w_provided) == 1) {
        stop("eim: If providing a matrix, 'X' and 'W' must be provided.")
    }

    if (v_provided && !xw_provided) {
        stop("eim: If providing 'V', you must also provide 'X' and 'W'.")
    }

    if (sum(xwv_provided, json_provided) != 1) {
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
        if (!is.null(matrices$V)) {
            V <- as.matrix(matrices$V)
        }
        allowed_params <- c(
            "prob",
            "avg_prob",
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
            "HET",
            "AE",
            "W_agg",
            "beta",
            "alpha",
            "sd_beta",
            "sd_alpha",
            "mcmc_samples",
            "mcmc_stepsize",
            "mvncdf_method",
            "mvncdf_samples",
            "mvncdf_error",
            "miniter",
            "adjust_prob_cond_method",
            "adjust_prob_cond_every",
            "prob_inv",
            "cond_prob_inv",
            "expected_outcome_inv",
            "LL_ind",
            "LL_rev_ind",
            "dLL",
            "dLL_rev",
            "symmetric_weights",
            "symmetric_weight_method"
        )
        extra_params <- matrices[names(matrices) %in% allowed_params] # TODO: Validate them
    }

    # Perform matricial validation
    .validate_eim(X, W, V) # nolint

    # Create the S3 object
    obj <- list(
        X = X,
        W = W
    )
    if (!is.null(V)) {
        obj$V <- V
    }

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
#' It supports both non-covariate and covariate models; the covariate mode is enabled by providing `V`.
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
#' - `saddlepoint`: Uses a saddlepoint approximation of the PMF to estimate the conditional probability.
#' - `mcmc`: Uses MCMC to sample vote outcomes. This is used to estimate the conditional probability of the E-step.
#' - `exact`: Solves the E-step using the Total Probability Law.
#'
#' When `V` is supplied (covariate mode), only `mult` is supported.
#'
#' For a detailed description of each method, see [fastei-package] and **References**.
#'
#' @param initial_prob An optional string specifying the method used to obtain the initial
#'   probability. Accepted values are:
#' - `uniform`: Assigns equal probability to every candidate within each group.
#' - `proportional`: Assigns probabilities to each group based on the proportion of candidates votes.
#' - `group_proportional`: Computes the probability matrix by taking into account both group and candidate proportions. This is the default method.
#' - `random`: Use randomized values to fill the probability matrix.
#' This argument is ignored if `V` is supplied (covariate mode), as the initial probabilities are computed with `alpha_init` and `beta_init`.
#'
#' @param mixture Positive integer indicating the number of latent voting profiles.
#'   If `mixture = 1` (default), the standard EM is used. If `mixture > 1`,
#'   a finite-mixture EM extension is applied with that number of components.
#'
#' @param H Positive integer indicating the number of latent profiles at the group
#'   level (row-level mixture). If `H = 1` (default), row-level mixture is disabled.
#'   If `H > 1`, the non-parametric EM uses a row-level finite-mixture backend in C.
#'   If `H > 1`, then `K` is implied by `H` and `mixture` is forced to `H`.
#'
#' @param HET Optional non-negative numeric threshold (in percentage points) used to
#'   adaptively increase the number of mixture components. If supplied, the model computes
#'   `HET = 50 * sum(|z_hat - z_tilde|) / sum_b I_b`,
#'   where `z_hat` is `expected_outcome` and
#'   `z_tilde[g,c,b] = W[b,g] * sum_k p[g,c,k] * r[b,k]`,
#'   with `r[b,k]` the posterior responsibility of ballot-box `b` for component `k`.
#'   If `H > 1`, the adaptive search is run over `H = 1, ..., 7` and the
#'   returned `K` tracks the selected row-level profile count (`K = H` in that path).
#'   Otherwise, the adaptive search is run over `K = 1, ..., 7` (with `H = 1`).
#'   While `HET >= alpha` (with `alpha = HET`), the model is re-estimated with the next
#'   value in the selected sequence. If no value up to 7 satisfies the threshold, the fit
#'   with the smallest observed `HET` is returned. If `HET = 0`, all values from 1 to 7 are
#'   always evaluated and the fit with minimum observed `HET` is returned.
#'   The returned object always includes the final `K`, `H`, and achieved `HET`.
#'
#' @param AE Optional non-negative numeric selector to activate adaptive search
#'   using an absolute-error criterion inspired by `mae_inverse`. For each fit, it computes
#'   `AE = tau = sum(|X - W %*% p|)`. If `H > 1`, the adaptive search is run over
#'   `H = 1, ..., 7` (with `K = H`); otherwise it is run over `K = 1, ..., 7` (with `H = 1`).
#'   When supplied, all candidate values in the selected sequence are evaluated and
#'   the fit with minimum observed `AE` is returned. This argument is mutually exclusive
#'   with `HET`.
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
#' @param beta_init Optional `g x (c-1)` matrix of initial group coefficients. Ignored if no covariates are provided (i.e., `V = NULL`).
#'
#' @param alpha_init Optional `(c-1) x a` matrix of initial attribute coefficients used for initialization. Ignored if no covariates are provided (`V = NULL`).
#'
#' @param maxnewton Maximum number of Newton iterations used in the parametric M-step. Default is 1. Ignored if no covariates are provided (i.e., `V = NULL`).
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
#' @param symmetric A boolean indicating whether to perform a symmetric estimation. If `TRUE`, the algorithm runs twice: first estimating the probabilities of candidates given groups, and then estimating the probabilities of groups given candidates. The final probabilities are obtained by combining the expected outcomes from both runs (equal average by default; see `symmetric_weight_method`). This approach can provide a more balanced estimation in certain scenarios. The default value is `FALSE`.
#'
#' @param symmetric_weight_method Character string indicating how to combine both directions when `symmetric = TRUE`. Valid values are `"average"` (default, equal weights `0.5/0.5`) and `"delta_ll"` (weights proportional to `dLL = LL - LL_ind` and `dLL_rev = LL_rev - LL_rev_ind`).
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
#'   \item{prob}{
#'     If `V` is `NULL` (non-covariate), the estimated global probability matrix `(g x c)`.
#'     If `V` is supplied (covariates), a `(g x c x b)` 3d-array of probabilities for each ballot-box.
#'   }
#' 	 \item{cond_prob}{
#'     A `(g x c x b)` 3d-array with the probability that at each ballot-box a voter of each group voted for each candidate,
#'     given the observed outcome at the particular ballot-box.
#'   }
#' 	 \item{expected_outcome}{
#'     A `(g x c x b)` 3d-array with the expected votes cast for each ballot-box.
#'   }
#'   \item{logLik}{The log-likelihood value from the last iteration.}
#'   \item{iterations}{The total number of iterations performed by the EM algorithm.}
#'   \item{time}{The total execution time of the algorithm in seconds.}
#'   \item{K}{The final number of latent profiles used by the fitted model.}
#'   \item{H}{The final number of row-level latent profiles used by the fitted model.}
#'   \item{HET}{The final heterogeneity metric achieved by the fitted model.}
#'   \item{AE}{The final absolute-error metric achieved by the fitted model.}
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
#' Furthermore, if `symmetric = TRUE`, the following additional attributes are included:
#' \describe{
#' 		\item{prob_inv}{The estimated probability matrix `(c x g)`, obtained by swapping `X` and `W`.}
#' 		\item{cond_prob_inv}{A `(c x g x b)` 3d-array with the probability that at each ballot-box a voter of each candidate voted for each group, given the observed outcome at the particular ballot-box.}
#' }
#' If `symmetric_weight_method = "delta_ll"` and the model is non-parametric with a single profile (`K = H = 1`), the object also includes `LL_ind`, `LL_rev_ind`, `dLL`, `dLL_rev`, and `symmetric_weights`.
#' Under this argument, the conditional probabilities will be obtained by dividing new expected outcomes by `W`. The probabilities will be calculated by performing an M-step.
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
                   V = NULL,
                   json_path = NULL,
                   method = "mult",
                   initial_prob = "group_proportional",
                   mixture = 1,
                   H = 1,
                   HET = NULL,
                   AE = NULL,
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
                   maxnewton = 1,
                   beta_init = NULL,
                   alpha_init = NULL,
                   scale_factor = 1,
                   symmetric = FALSE,
                   symmetric_weight_method = "average",
                   ...) {
    base_call <- match.call()
    caller_env <- parent.frame()
    # Debug toggle: if `mixture` is explicitly passed (even as 1), route through mixture backend.
    mixture_explicit <- "mixture" %in% names(as.list(base_call))
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint

    if (!is.null(seed)) {
        set.seed(seed)
    }

    if (is.null(object)) {
        object <- eim(X = X, W = W, V = V, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("run_em: The object must be initialized with the eim() function.")
    }
    mixture <- as.integer(mixture)
    H <- as.integer(H)
    if (!is.null(HET)) {
        HET <- as.numeric(HET)
    }
    if (!is.null(AE)) {
        AE <- as.numeric(AE)
    }
    if (!is.null(V)) {
        object$V <- V
    }


    # Applies a scaling
    if (scale_factor != 1) {
        object$X <- round(object$X / all_params$scale_factor)
        object$W <- round(object$W / all_params$scale_factor)
    }

    compute_HET_metric <- function(fit_object) {
        if (is.null(fit_object$expected_outcome)) {
            return(NA_real_)
        }

        W_fit <- if (is.null(fit_object$W_agg)) fit_object$W else fit_object$W_agg
        if (is.null(W_fit)) {
            return(NA_real_)
        }

        denom <- sum(W_fit)
        if (!is.finite(denom) || denom <= 0) {
            return(NA_real_)
        }

        z_hat <- fit_object$expected_outcome
        if (!is.array(z_hat) || length(dim(z_hat)) != 3) {
            return(NA_real_)
        }

        z_tilde <- NULL
        if (!is.null(fit_object$component_prob) &&
            !is.null(fit_object$responsibilities)) {
            component_prob <- fit_object$component_prob
            resp <- as.matrix(fit_object$responsibilities)
            if (length(dim(component_prob)) == 3 &&
                ncol(resp) == dim(component_prob)[3]) {
                G <- dim(component_prob)[1]
                C <- dim(component_prob)[2]
                B <- nrow(W_fit)
                z_tilde <- array(0, dim = c(G, C, B))
                for (b in seq_len(B)) {
                    for (g in seq_len(G)) {
                        for (c in seq_len(C)) {
                            z_tilde[g, c, b] <- W_fit[b, g] * sum(component_prob[g, c, ] * resp[b, ])
                        }
                    }
                }
            }
        }

        if (is.null(z_tilde) && is.array(fit_object$prob) && length(dim(fit_object$prob)) == 3) {
            prob_arr <- fit_object$prob
            G <- dim(prob_arr)[1]
            C <- dim(prob_arr)[2]
            B <- dim(prob_arr)[3]
            z_tilde <- array(0, dim = c(G, C, B))
            for (b in seq_len(B)) {
                for (g in seq_len(G)) {
                    z_tilde[g, , b] <- W_fit[b, g] * prob_arr[g, , b]
                }
            }
        }

        if (is.null(z_tilde) && is.matrix(fit_object$prob)) {
            p_tilde <- fit_object$prob
            G <- nrow(p_tilde)
            C <- ncol(p_tilde)
            B <- nrow(W_fit)
            z_tilde <- array(0, dim = c(G, C, B))
            for (b in seq_len(B)) {
                for (g in seq_len(G)) {
                    z_tilde[g, , b] <- W_fit[b, g] * p_tilde[g, ]
                }
            }
        }

        if (is.null(z_tilde)) {
            return(NA_real_)
        }

        if (!identical(dim(z_hat), dim(z_tilde))) {
            return(NA_real_)
        }

        het <- 50 * sum(abs(z_hat - z_tilde)) / denom
        if (!is.finite(het)) {
            return(NA_real_)
        }
        het
    }

    compute_AE_metric <- function(fit_object) {
        W_fit <- if (is.null(fit_object$W_agg)) fit_object$W else fit_object$W_agg
        X_fit <- fit_object$X
        if (is.null(W_fit) || is.null(X_fit)) {
            return(NA_real_)
        }

        prob_forward <- fit_object$prob
        if (is.array(prob_forward) && length(dim(prob_forward)) == 3) {
            G <- dim(prob_forward)[1]
            C <- dim(prob_forward)[2]
            prob_mat <- matrix(0, nrow = G, ncol = C)
            for (g in seq_len(G)) {
                den_g <- sum(W_fit[, g])
                if (is.finite(den_g) && den_g > 0) {
                    prob_mat[g, ] <- colSums(t(prob_forward[g, , ]) * W_fit[, g]) / den_g
                } else {
                    prob_mat[g, ] <- rep(1 / C, C)
                }
            }
            prob_forward <- prob_mat
        }
        if (!is.matrix(prob_forward)) {
            return(NA_real_)
        }

        x_hat <- as.matrix(W_fit) %*% prob_forward
        tau <- sum(abs(as.matrix(X_fit) - x_hat))
        if (!is.finite(tau)) {
            return(NA_real_)
        }
        tau
    }

    normalize_prob_rows_run_em <- function(prob_matrix) {
        prob_matrix[!is.finite(prob_matrix) | prob_matrix < 0] <- 0
        row_sums <- rowSums(prob_matrix)
        valid_rows <- is.finite(row_sums) & row_sums > 0
        if (any(valid_rows)) {
            prob_matrix[valid_rows, ] <- sweep(prob_matrix[valid_rows, , drop = FALSE], 1, row_sums[valid_rows], "/")
        }
        if (any(!valid_rows)) {
            prob_matrix[!valid_rows, ] <- rep(1 / ncol(prob_matrix), ncol(prob_matrix))
        }
        prob_matrix
    }

    zero_vote_state <- NULL
    restore_zero_vote_columns <- function(fit_object) {
        if (is.null(zero_vote_state)) {
            return(fit_object)
        }

        original_X <- zero_vote_state$original_X
        keep_cols <- zero_vote_state$keep_cols
        total_candidates <- ncol(original_X)
        candidate_names <- colnames(original_X)

        restore_candidate_second_matrix <- function(mat) {
            if (is.null(mat) || !is.matrix(mat)) {
                return(mat)
            }
            restored <- matrix(0, nrow = nrow(mat), ncol = total_candidates)
            restored[, keep_cols] <- mat
            rownames(restored) <- rownames(mat)
            colnames(restored) <- candidate_names
            restored
        }

        restore_candidate_first_matrix <- function(mat) {
            if (is.null(mat) || !is.matrix(mat)) {
                return(mat)
            }
            restored <- matrix(0, nrow = total_candidates, ncol = ncol(mat))
            restored[keep_cols, ] <- mat
            rownames(restored) <- candidate_names
            colnames(restored) <- colnames(mat)
            restored
        }

        restore_candidate_second_array <- function(arr) {
            if (is.null(arr) || !is.array(arr) || length(dim(arr)) != 3) {
                return(arr)
            }
            arr_dimnames <- dimnames(arr)
            restored <- array(0, dim = c(dim(arr)[1], total_candidates, dim(arr)[3]))
            restored[, keep_cols, ] <- arr
            dimnames(restored) <- list(arr_dimnames[[1]], candidate_names, arr_dimnames[[3]])
            restored
        }

        restore_candidate_first_array <- function(arr) {
            if (is.null(arr) || !is.array(arr) || length(dim(arr)) != 3) {
                return(arr)
            }
            arr_dimnames <- dimnames(arr)
            restored <- array(0, dim = c(total_candidates, dim(arr)[2], dim(arr)[3]))
            restored[keep_cols, , ] <- arr
            dimnames(restored) <- list(candidate_names, arr_dimnames[[2]], arr_dimnames[[3]])
            restored
        }

        fit_object$X <- original_X
        fit_object$prob <- if (is.array(fit_object$prob) && length(dim(fit_object$prob)) == 3) {
            restore_candidate_second_array(fit_object$prob)
        } else {
            restore_candidate_second_matrix(fit_object$prob)
        }
        fit_object$cond_prob <- restore_candidate_second_array(fit_object$cond_prob)
        fit_object$expected_outcome <- restore_candidate_second_array(fit_object$expected_outcome)
        fit_object$component_prob <- restore_candidate_second_array(fit_object$component_prob)
        fit_object$prob_inv <- if (is.array(fit_object$prob_inv) && length(dim(fit_object$prob_inv)) == 3) {
            restore_candidate_first_array(fit_object$prob_inv)
        } else {
            restore_candidate_first_matrix(fit_object$prob_inv)
        }
        fit_object$cond_prob_inv <- restore_candidate_first_array(fit_object$cond_prob_inv)
        fit_object$expected_outcome_inv <- restore_candidate_first_array(fit_object$expected_outcome_inv)

        if (is.matrix(fit_object$initial_prob)) {
            fit_object$initial_prob <- restore_candidate_second_matrix(fit_object$initial_prob)
        }

        if (isTRUE(zero_vote_state$can_restore_coefficients) && is.matrix(fit_object$beta)) {
            coef_keep_cols <- keep_cols[-length(keep_cols)]
            restored_beta <- matrix(0, nrow = nrow(fit_object$beta), ncol = total_candidates - 1)
            if (length(coef_keep_cols) > 0 && ncol(fit_object$beta) > 0) {
                restored_beta[, coef_keep_cols] <- fit_object$beta
            }
            rownames(restored_beta) <- rownames(fit_object$beta)
            colnames(restored_beta) <- candidate_names[-total_candidates]
            fit_object$beta <- restored_beta
        }

        if (isTRUE(zero_vote_state$can_restore_coefficients) && is.matrix(fit_object$alpha)) {
            coef_keep_cols <- keep_cols[-length(keep_cols)]
            restored_alpha <- matrix(0, nrow = total_candidates - 1, ncol = ncol(fit_object$alpha))
            if (length(coef_keep_cols) > 0 && nrow(fit_object$alpha) > 0) {
                restored_alpha[coef_keep_cols, ] <- fit_object$alpha
            }
            rownames(restored_alpha) <- candidate_names[-total_candidates]
            colnames(restored_alpha) <- colnames(fit_object$alpha)
            fit_object$alpha <- restored_alpha
        }

        fit_object
    }

    finalize_fit_object <- function(fit_object) {
        fit_object <- restore_zero_vote_columns(fit_object)
        fit_object$HET <- compute_HET_metric(fit_object)
        fit_object$AE <- compute_AE_metric(fit_object)
        invisible(fit_object)
    }

    is_parametric <- !is.null(object$V)

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

    full_candidate_object <- object

    zero_vote_cols <- which(colSums(object$X) == 0)
    if (length(zero_vote_cols) > 0) {
        keep_cols <- which(colSums(object$X) > 0)
        if (length(keep_cols) == 0) {
            stop("run_em: 'X' must contain at least one candidate with a positive vote total.")
        }

        zero_vote_state <- list(
            original_X = object$X,
            keep_cols = keep_cols,
            can_restore_coefficients = is_parametric && tail(keep_cols, 1) == ncol(object$X)
        )
        object$X <- object$X[, keep_cols, drop = FALSE]

        if (is.matrix(initial_prob)) {
            initial_prob <- normalize_prob_rows_run_em(initial_prob[, keep_cols, drop = FALSE])
            all_params$initial_prob <- initial_prob
        }

        if (is_parametric && zero_vote_state$can_restore_coefficients) {
            coef_keep_cols <- keep_cols[-length(keep_cols)]
            if (!is.null(beta_init)) {
                beta_init <- beta_init[, coef_keep_cols, drop = FALSE]
            }
            if (!is.null(alpha_init)) {
                alpha_init <- alpha_init[coef_keep_cols, , drop = FALSE]
            }
        } else if (is_parametric && (!is.null(beta_init) || !is.null(alpha_init))) {
            beta_init <- NULL
            alpha_init <- NULL
        }
    }

    W_effective <- if (is.null(object$W_agg)) object$W else object$W_agg
    G_effective <- ncol(W_effective)

    # Row-level mixture implies K = H.
    if (H > 1L) {
        if (mixture != H) {
            warning(sprintf(
                "run_em: `H > 1` implies `K = H`. Setting `mixture` from %d to %d.",
                mixture,
                H
            ))
        }
        mixture <- as.integer(H)
    }

    object$method <- method

    if (is_parametric) {
        if (method != "mult") {
            stop("run_em: Parametric mode only supports method = \"mult\".")
        }
        if (mixture > 1) {
            stop("run_em: `mixture > 1` is currently supported only for non-parametric models (`V = NULL`).")
        }
        if (!is.null(HET)) {
            stop("run_em: `HET` is currently supported only for non-parametric models (`V = NULL`).")
        }
        if (!is.null(AE)) {
            stop("run_em: `AE` is currently supported only for non-parametric models (`V = NULL`).")
        }
        if (H > 1) {
            stop("run_em: `H > 1` is currently supported only for non-parametric models (`V = NULL`).")
        }

        W <- if (is.null(object$W_agg)) object$W else object$W_agg
        V <- object$V
        num_candidates <- ncol(object$X)
        num_groups <- ncol(W)
        num_attributes <- ncol(V)

        if (is.null(beta_init)) {
            beta <- matrix(0, nrow = num_groups, ncol = num_candidates - 1)
        } else {
            beta <- beta_init
        }
        if (is.null(alpha_init)) {
            alpha <- matrix(0, nrow = num_candidates - 1, ncol = num_attributes)
        } else {
            alpha <- alpha_init
        }

        if (!is.matrix(beta) || nrow(beta) != num_groups || ncol(beta) != num_candidates - 1) {
            stop("run_em: 'beta' must be a matrix with dimensions (g x (c-1)).")
        }
        if (!is.matrix(alpha) || nrow(alpha) != num_candidates - 1 || ncol(alpha) != num_attributes) {
            stop("run_em: 'alpha' must be a matrix with dimensions ((c-1) x a).")
        }

        resulting_values <- EMAlgorithmParametric(
            as.matrix(object$X),
            as.matrix(W),
            as.matrix(V),
            as.matrix(beta),
            as.matrix(alpha),
            maxiter,
            maxtime,
            ll_threshold,
            maxnewton,
            verbose,
            adjust_prob_cond_method,
            adjust_prob_cond_every
        )

        object$prob <- resulting_values$prob
        dimnames(object$prob) <- list(
            colnames(W),
            colnames(object$X),
            rownames(object$X)
        )
        object$cond_prob <- resulting_values$cond_prob
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
        object$beta <- resulting_values$beta
        if (!is.null(colnames(W))) {
            rownames(object$beta) <- colnames(W)
        }
        if (!is.null(colnames(object$X))) {
            colnames(object$beta) <- colnames(object$X)[-ncol(object$X)]
        }
        object$alpha <- resulting_values$alpha
        if (!is.null(colnames(object$X))) {
            rownames(object$alpha) <- colnames(object$X)[-ncol(object$X)]
        }
        if (!is.null(colnames(V))) {
            colnames(object$alpha) <- colnames(V)
        }
        object$iterations <- as.numeric(resulting_values$iter)
        object$logLik <- as.numeric(resulting_values$logLik)
        object$time <- resulting_values$time

        object$maxiter <- maxiter
        object$maxtime <- maxtime
        object$ll_threshold <- ll_threshold
        object$maxnewton <- maxnewton
        object$adjust_prob_cond_method <- adjust_prob_cond_method
        object$adjust_prob_cond_every <- adjust_prob_cond_every
        object$mixture <- as.integer(mixture)
        object$K <- as.integer(mixture)
        object$H <- as.integer(H)

        if (symmetric) {
            W_sym <- if (is.null(object$W_agg)) object$W else object$W_agg
            object <- .run_em_symmetric_helper(
                object = object,
                base_call = base_call,
                parent_env = parent.frame(),
                W_sym = W_sym,
                all_params = all_params,
                include_V = TRUE,
                reset_parametric_init = TRUE,
                transform_initial_prob = FALSE,
                prob_from_cond_prob = TRUE
            )
        }

        return(finalize_fit_object(object))
    }

    adaptive_search_metric <- function(metric_name, threshold) {
        metric_fun <- switch(metric_name,
            HET = compute_HET_metric,
            AE = compute_AE_metric,
            stop("run_em: Unsupported adaptive metric.")
        )

        K_max <- 7L
        best_fit <- NULL
        best_metric <- Inf
        last_fit <- NULL
        search_best_only <- isTRUE(all.equal(threshold, 0)) || identical(metric_name, "AE")
        search_over_H <- H > 1L

        if (search_over_H) {
            H_limit <- as.integer(K_max)
            H_current <- 1L
            repeat {
                base_call_metric <- base_call
                base_call_metric$object <- full_candidate_object
                base_call_metric$X <- NULL
                base_call_metric$W <- NULL
                base_call_metric$V <- NULL
                base_call_metric$json_path <- NULL
                base_call_metric$mixture <- H_current
                base_call_metric$H <- H_current
                base_call_metric$HET <- NULL
                base_call_metric$AE <- NULL

                fit_metric <- eval(base_call_metric, caller_env)
                last_fit <- fit_metric
                metric_value <- metric_fun(fit_metric)
                fit_metric[[metric_name]] <- metric_value
                fit_metric$K <- as.integer(H_current)
                fit_metric$mixture <- as.integer(H_current)
                fit_metric$H <- as.integer(H_current)

                if (is.finite(metric_value) && metric_value < best_metric) {
                    best_metric <- metric_value
                    best_fit <- fit_metric
                }

                if (!search_best_only && is.finite(metric_value) && metric_value < threshold) {
                    return(fit_metric)
                }

                if (H_current >= H_limit) {
                    if (!is.null(best_fit)) {
                        if (!search_best_only) {
                            warning(sprintf(
                                "run_em: No H <= %d satisfied the %s threshold. Returning the fit with minimum %s.",
                                H_limit,
                                metric_name,
                                metric_name
                            ))
                        }
                        return(best_fit)
                    }
                    warning(sprintf(
                        "run_em: No finite %s found for H <= %d. Returning the last fitted model.",
                        metric_name,
                        H_limit
                    ))
                    return(last_fit)
                }

                H_current <- H_current + 1L
            }
        } else {
            K_current <- 1L
            repeat {
                base_call_metric <- base_call
                base_call_metric$object <- full_candidate_object
                base_call_metric$X <- NULL
                base_call_metric$W <- NULL
                base_call_metric$V <- NULL
                base_call_metric$json_path <- NULL
                base_call_metric$mixture <- K_current
                base_call_metric$H <- 1L
                base_call_metric$HET <- NULL
                base_call_metric$AE <- NULL

                fit_metric <- eval(base_call_metric, caller_env)
                last_fit <- fit_metric
                metric_value <- metric_fun(fit_metric)
                fit_metric[[metric_name]] <- metric_value
                fit_metric$K <- as.integer(K_current)
                fit_metric$mixture <- as.integer(K_current)
                fit_metric$H <- 1L

                if (is.finite(metric_value) && metric_value < best_metric) {
                    best_metric <- metric_value
                    best_fit <- fit_metric
                }

                if (!search_best_only && is.finite(metric_value) && metric_value < threshold) {
                    return(fit_metric)
                }

                if (K_current >= K_max) {
                    if (!is.null(best_fit)) {
                        if (!search_best_only) {
                            warning(sprintf(
                                "run_em: No K <= %d satisfied the %s threshold. Returning the fit with minimum %s.",
                                K_max,
                                metric_name,
                                metric_name
                            ))
                        }
                        return(best_fit)
                    }
                    warning(sprintf(
                        "run_em: No finite %s found for K <= %d. Returning the last fitted model.",
                        metric_name,
                        K_max
                    ))
                    return(last_fit)
                }

                K_current <- K_current + 1L
            }
        }
    }

    if (!is.null(HET)) {
        return(invisible(adaptive_search_metric("HET", HET)))
    }
    if (!is.null(AE)) {
        return(invisible(adaptive_search_metric("AE", AE)))
    }

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

    K <- as.integer(mixture)
    H_group <- as.integer(H)
    use_mixture <- (K > 1L) || (K == 1L && mixture_explicit)
    use_row_mixture <- H_group > 1

    if (use_row_mixture) {
        resulting_values <- EMAlgorithmRowMixture(
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
            if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1),
            H_group
        )

        object$cond_prob <- resulting_values$q
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
            object$logLik <- as.numeric(resulting_values$log_likelihood)
        }
        object$time <- resulting_values$total_time
        object$message <- resulting_values$stopping_reason
        object$status <- as.integer(resulting_values$finish_id)
        object$mixture <- as.integer(H_group)
        object$K <- as.integer(H_group)
        object$H <- as.integer(H_group)

        object$phi <- as.matrix(resulting_values$phi)
        dimnames(object$phi) <- list(colnames(W), paste0("H", seq_len(H_group)))

        object$responsibilities <- as.matrix(resulting_values$responsibilities)
        assignment_count <- ncol(object$responsibilities)
        dimnames(object$responsibilities) <- list(rownames(object$X), paste0("A", seq_len(assignment_count)))

        object$component_prob <- resulting_values$component_prob
        dimnames(object$component_prob) <- list(colnames(W), colnames(object$X), paste0("H", seq_len(H_group)))

        object$miniter <- miniter
        object$maxiter <- maxiter
        object$maxtime <- maxtime
        object$param_threshold <- param_threshold
        object$ll_threshold <- ll_threshold
        object$initial_prob <- initial_prob
        object$adjust_prob_cond_method <- adjust_prob_cond_method
        object$adjust_prob_cond_every <- adjust_prob_cond_every

        if (symmetric) {
            W_sym <- if (is.null(object$W_agg)) object$W else object$W_agg
            object <- .run_em_symmetric_helper(
                object = object,
                base_call = base_call,
                parent_env = parent.frame(),
                W_sym = W_sym,
                all_params = all_params,
                include_V = FALSE,
                reset_parametric_init = FALSE,
                transform_initial_prob = TRUE,
                prob_from_cond_prob = FALSE
            )
        }

        return(finalize_fit_object(object))
    }

    if (use_mixture) {
        resulting_values <- EMAlgorithmMixture(
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
            if (is.matrix(initial_prob)) initial_prob else matrix(-1, nrow = 1, ncol = 1),
            K
        )

        object$cond_prob <- resulting_values$q
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
            object$logLik <- as.numeric(resulting_values$log_likelihood)
        }
        object$time <- resulting_values$total_time
        object$message <- resulting_values$stopping_reason
        object$status <- as.integer(resulting_values$finish_id)
        object$mixture <- as.integer(mixture)
        object$K <- as.integer(K)
        object$H <- as.integer(H_group)
        object$phi <- as.numeric(resulting_values$phi)
        names(object$phi) <- paste0("H", seq_len(K))
        object$responsibilities <- as.matrix(resulting_values$responsibilities)
        dimnames(object$responsibilities) <- list(rownames(object$X), paste0("H", seq_len(K)))
        object$component_prob <- resulting_values$component_prob
        dimnames(object$component_prob) <- list(colnames(W), colnames(object$X), paste0("H", seq_len(K)))

        object$miniter <- miniter
        object$maxiter <- maxiter
        object$maxtime <- maxtime
        object$param_threshold <- param_threshold
        object$ll_threshold <- ll_threshold
        object$initial_prob <- initial_prob
        object$adjust_prob_cond_method <- adjust_prob_cond_method
        object$adjust_prob_cond_every <- adjust_prob_cond_every

        if (symmetric) {
            W_sym <- if (is.null(object$W_agg)) object$W else object$W_agg
            object <- .run_em_symmetric_helper(
                object = object,
                base_call = base_call,
                parent_env = parent.frame(),
                W_sym = W_sym,
                all_params = all_params,
                include_V = FALSE,
                reset_parametric_init = FALSE,
                transform_initial_prob = TRUE,
                prob_from_cond_prob = FALSE
            )
        }

        return(finalize_fit_object(object))
    }

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
    object$mixture <- as.integer(mixture)
    object$K <- as.integer(mixture)
    object$H <- as.integer(H_group)
    object$adjust_prob_cond_method <- adjust_prob_cond_method
    object$adjust_prob_cond_every <- adjust_prob_cond_every

    if (symmetric) {
        object <- .run_em_symmetric_helper(
            object = object,
            base_call = base_call,
            parent_env = parent.frame(),
            W_sym = object$W,
            all_params = all_params,
            include_V = FALSE,
            reset_parametric_init = FALSE,
            transform_initial_prob = TRUE,
            prob_from_cond_prob = FALSE
        )
    }

    finalize_fit_object(object)
}

.run_em_symmetric_helper <- function(object,
                                     base_call,
                                     parent_env,
                                     W_sym,
                                     all_params = NULL,
                                     include_V = FALSE,
                                     reset_parametric_init = FALSE,
                                     transform_initial_prob = FALSE,
                                     prob_from_cond_prob = FALSE) {
    base_call_sym <- base_call
    base_call_sym$symmetric <- FALSE
    base_call_sym$X <- object$W
    base_call_sym$W <- object$X
    base_call_sym$json_path <- NULL
    base_call_sym$object <- NULL
    base_call_sym$scale_factor <- 1

    if (include_V) {
        base_call_sym$V <- object$V
    }

    if (reset_parametric_init) {
        base_call_sym$beta_init <- NULL
        base_call_sym$alpha_init <- NULL
    }

    if (transform_initial_prob && !is.null(all_params)) {
        ip <- all_params$initial_prob
        if (is.matrix(ip)) {
            col_tot_X <- colSums(object$X)
            num <- sweep(ip, 2, col_tot_X, "*")
            denom <- rowSums(num)
            base_call_sym$initial_prob <- sweep(num, 1, denom, "/")
            base_call_sym$initial_prob <- t(base_call_sym$initial_prob)
        }
    }

    inverse <- eval(base_call_sym, parent_env)
    object$cond_prob_inv <- inverse$cond_prob
    object$prob_inv <- inverse$prob
    object$expected_outcome_inv <- inverse$expected_outcome
    object$time <- object$time + inverse$time
    object$iterations <- object$iterations + inverse$iterations

    normalize_prob_rows <- function(prob_matrix) {
        prob_matrix[!is.finite(prob_matrix) | prob_matrix < 0] <- 0
        rs <- rowSums(prob_matrix)
        valid <- is.finite(rs) & rs > 0
        if (any(valid)) {
            prob_matrix[valid, ] <- sweep(prob_matrix[valid, , drop = FALSE], 1, rs[valid], "/")
        }
        if (any(!valid)) {
            prob_matrix[!valid, ] <- rep(1 / ncol(prob_matrix), ncol(prob_matrix))
        }
        prob_matrix
    }

    normalize_cube_last_dim <- function(arr3) {
        B <- dim(arr3)[1]
        R <- dim(arr3)[2]
        K <- dim(arr3)[3]
        for (b in seq_len(B)) {
            for (r in seq_len(R)) {
                row_vals <- arr3[b, r, ]
                row_vals[!is.finite(row_vals) | row_vals < 0] <- 0
                row_sum <- sum(row_vals)
                if (is.finite(row_sum) && row_sum > 0) {
                    arr3[b, r, ] <- row_vals / row_sum
                } else {
                    arr3[b, r, ] <- rep(1 / K, K)
                }
            }
        }
        arr3
    }

    mstep_from_q <- function(q_array, W_matrix) {
        q_bgc <- aperm(q_array, c(3, 1, 2))
        weighted <- sweep(q_bgc, c(1, 2), W_matrix, "*")
        num <- apply(weighted, c(2, 3), sum)
        den <- colSums(W_matrix)
        prob <- sweep(num, 1, den, "/")
        normalize_prob_rows(prob)
    }

    em_loglik_from_prob <- function(X_matrix, W_matrix, prob_matrix, method_name) {
        if (!is.matrix(prob_matrix)) {
            return(NA_real_)
        }
        if (ncol(X_matrix) != ncol(prob_matrix) || ncol(W_matrix) != nrow(prob_matrix)) {
            return(NA_real_)
        }
        ll_result <- EMAlgorithmFull(
            t(X_matrix),
            W_matrix,
            as.character(method_name),
            "custom",
            as.integer(0),
            as.double(0),
            as.double(0),
            as.double(-Inf),
            as.logical(TRUE),
            as.logical(FALSE),
            as.integer(3000),
            as.integer(1000),
            "genz",
            as.double(1e-3),
            as.integer(5000),
            as.integer(0),
            "",
            as.logical(FALSE),
            prob_matrix
        )
        as.numeric(ll_result$log_likelihood)
    }

    weight_method <- "average"
    if (!is.null(all_params) && !is.null(all_params$symmetric_weight_method)) {
        weight_method <- as.character(all_params$symmetric_weight_method)
    }
    em_method <- "mult"
    if (!is.null(all_params) && !is.null(all_params$method)) {
        em_method <- as.character(all_params$method)
    }

    weight_original <- 0.5
    weight_reverse <- 0.5
    expected_orig <- object$expected_outcome
    A2 <- aperm(inverse$expected_outcome, c(2, 1, 3))

    single_profile <- (is.null(object$K) || as.integer(object$K) == 1L) &&
        (is.null(object$H) || as.integer(object$H) == 1L)
    can_use_delta <- identical(weight_method, "delta_ll") &&
        !include_V &&
        single_profile &&
        is.matrix(object$X) &&
        is.matrix(W_sym) &&
        is.matrix(object$W) &&
        is.matrix(object$prob) &&
        is.matrix(inverse$prob) &&
        nrow(object$X) == nrow(W_sym) &&
        nrow(object$W) == nrow(object$X) &&
        ncol(W_sym) == ncol(object$W) &&
        identical(dim(expected_orig), dim(A2))

    if (can_use_delta) {
        q_orig_bgc <- aperm(object$cond_prob, c(3, 1, 2))
        z_from_orig_bgc <- sweep(q_orig_bgc, c(1, 2), W_sym, "*")
        q_rev_ind_bcg <- sweep(aperm(z_from_orig_bgc, c(1, 3, 2)), c(1, 2), object$X, "/")
        q_rev_ind_bcg <- normalize_cube_last_dim(q_rev_ind_bcg)
        q_rev_ind <- aperm(q_rev_ind_bcg, c(2, 3, 1))
        p_rev_ind <- mstep_from_q(q_rev_ind, object$X)

        q_rev_bcg <- aperm(inverse$cond_prob, c(3, 1, 2))
        z_from_rev_bgc <- aperm(sweep(q_rev_bcg, c(1, 2), object$X, "*"), c(1, 3, 2))
        q_ind_bgc <- sweep(z_from_rev_bgc, c(1, 2), W_sym, "/")
        q_ind_bgc <- normalize_cube_last_dim(q_ind_bgc)
        q_ind <- aperm(q_ind_bgc, c(2, 3, 1))
        p_ind <- mstep_from_q(q_ind, W_sym)

        LL <- suppressWarnings(as.numeric(object$logLik))
        if (!is.finite(LL)) {
            LL <- em_loglik_from_prob(object$X, W_sym, object$prob, em_method)
        }
        LL_rev <- suppressWarnings(as.numeric(inverse$logLik))
        if (!is.finite(LL_rev)) {
            LL_rev <- em_loglik_from_prob(object$W, object$X, inverse$prob, em_method)
        }

        LL_ind <- em_loglik_from_prob(object$X, W_sym, p_ind, em_method)
        LL_rev_ind <- em_loglik_from_prob(object$W, object$X, p_rev_ind, em_method)

        dLL <- LL - LL_ind
        dLL_rev <- LL_rev - LL_rev_ind

        dLL_pos <- max(0, dLL)
        dLL_rev_pos <- max(0, dLL_rev)
        den_dLL <- dLL_pos + dLL_rev_pos
        if (is.finite(den_dLL) && den_dLL > 0) {
            weight_original <- dLL_pos / den_dLL
            weight_reverse <- dLL_rev_pos / den_dLL
        }

        object$LL_ind <- LL_ind
        object$LL_rev_ind <- LL_rev_ind
        object$dLL <- dLL
        object$dLL_rev <- dLL_rev
        object$symmetric_weight_method <- "delta_ll"
    } else {
        object$symmetric_weight_method <- "average"
    }
    object$symmetric_weights <- c(original = weight_original, reverse = weight_reverse)

    object$expected_outcome <- weight_original * expected_orig + weight_reverse * A2
    dimnames(object$expected_outcome) <- list(
        colnames(W_sym),
        colnames(object$X),
        rownames(object$X)
    )

    E_bgc <- aperm(object$expected_outcome, c(3, 1, 2))
    Q_bgc <- sweep(E_bgc, c(1, 2), W_sym, "/")
    B <- dim(Q_bgc)[1]
    G <- dim(Q_bgc)[2]
    C <- dim(Q_bgc)[3]
    for (b in seq_len(B)) {
        for (g in seq_len(G)) {
            row_vals <- Q_bgc[b, g, ]
            row_vals[!is.finite(row_vals) | row_vals < 0] <- 0
            row_sum <- sum(row_vals)
            if (is.finite(row_sum) && row_sum > 0) {
                Q_bgc[b, g, ] <- row_vals / row_sum
            } else {
                Q_bgc[b, g, ] <- rep(1 / C, C)
            }
        }
    }
    object$cond_prob <- aperm(Q_bgc, c(2, 3, 1))
    dimnames(object$cond_prob) <- list(
        colnames(W_sym),
        colnames(object$X),
        rownames(object$X)
    )

    if (prob_from_cond_prob) {
        object$prob <- object$cond_prob
        dimnames(object$prob) <- list(
            colnames(W_sym),
            colnames(object$X),
            rownames(object$X)
        )
    } else {
        num <- apply(object$expected_outcome, c(1, 2), sum)
        den <- colSums(W_sym)
        object$prob <- sweep(num, 1, den, "/")
        object$prob[!is.finite(object$prob) | object$prob < 0] <- 0
        rs <- rowSums(object$prob)
        valid <- is.finite(rs) & rs > 0
        if (any(valid)) {
            object$prob[valid, ] <- sweep(object$prob[valid, , drop = FALSE], 1, rs[valid], "/")
        }
        if (any(!valid)) {
            object$prob[!valid, ] <- rep(1 / ncol(object$prob), ncol(object$prob))
        }
        dimnames(object$prob) <- list(colnames(W_sym), colnames(object$X))
    }

    object
}

#' Runs a Bootstrap to Estimate the **Standard Deviation** of Predicted Probabilities
#'
#' @description
#' This function computes the Expected-Maximization (EM) algorithm "`nboot`" times. It then computes the standard deviation from the `nboot` estimated probability matrices on each component.
#' It supports both non-parametric and parametric models; the parametric mode is enabled by providing `V` and only supports `method = "mult"`.
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
#' @seealso The [eim] object and [run_em] implementation.
#'
#' @return
#' Returns an `eim` object with the `sd` field containing the estimated standard deviations of the probabilities and the `avg_prob` field with the average bootstrapped probability matrix. If an `eim` object is provided, its attributes (see [run_em]) are retained in the returned object.
#'
#' For parametric models, it returns `sd_beta` and `sd_alpha` instead of `sd` and `avg_prob`.
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
                      V = NULL,
                      json_path = NULL,
                      nboot = 100,
                      allow_mismatch = TRUE,
                      seed = NULL,
                      maxnewton = 1,
                      ...) {
    # Retrieve the default values from run_em() as a list
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint # It would validate nboot too.

    # Initialize eim object if needed
    if (is.null(object)) {
        object <- eim(X = X, W = W, V = V, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("Bootstrap: The object must be initialized with the `eim()` function.")
    }
    if (!is.null(V)) {
        object$V <- V
    }

    is_parametric <- !is.null(object$V)

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

    if (is_parametric) {
        if (method != "mult") {
            stop("bootstrap: Parametric mode only supports method = \"mult\".")
        }

        # Applies a scaling
        if (!is.null(all_params$scale_factor) && all_params$scale_factor != 1) {
            object$X <- round(object$X / all_params$scale_factor)
            object$W <- round(object$W / all_params$scale_factor)
        }

        if (!allow_mismatch) {
            mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))
            if (length(mismatch_rows) > 0) {
                stop(
                    "bootstrap: Row-wise mismatch in vote totals detected.\n",
                    "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
                    "To allow mismatches, set `allow_mismatch = TRUE`."
                )
            }
        }

        # Set seed for reproducibility
        if (!is.null(seed)) set.seed(seed)

        W <- if (is.null(object$W_agg)) object$W else object$W_agg
        V <- object$V
        num_candidates <- ncol(object$X)
        num_groups <- ncol(W)
        num_attributes <- ncol(V)

        if (is.null(all_params$init_beta)) {
            beta <- matrix(0, nrow = num_groups, ncol = num_candidates - 1)
        }
        if (is.null(all_params$init_alpha)) {
            alpha <- matrix(0, nrow = num_candidates - 1, ncol = num_attributes)
        }

        if (!is.matrix(beta) || nrow(beta) != num_groups || ncol(beta) != num_candidates - 1) {
            stop("bootstrap: 'init_beta' must be a matrix with dimensions (g x (c-1)).")
        }
        if (!is.matrix(alpha) || nrow(alpha) != num_candidates - 1 || ncol(alpha) != num_attributes) {
            stop("bootstrap: 'init_alpha' must be a matrix with dimensions ((c-1) x a).")
        }

        maxiter <- if (!is.null(all_params$maxiter)) all_params$maxiter else 1000
        maxtime <- if (!is.null(all_params$maxtime)) all_params$maxtime else 3600
        verbose <- if (!is.null(all_params$verbose)) all_params$verbose else FALSE
        adjust_prob_cond_method <- if (!is.null(all_params$adjust_prob_cond_method)) all_params$adjust_prob_cond_method else "project_lp"
        adjust_prob_cond_every <- if (!is.null(all_params$adjust_prob_cond_every)) all_params$adjust_prob_cond_every else FALSE

        if ("ll_threshold" %in% names(all_params)) {
            ll_threshold <- all_params$ll_threshold
        } else {
            ll_threshold <- as.double(-Inf)
        }

        result <- bootstrapParametricAlg(
            as.matrix(object$X),
            as.matrix(W),
            as.matrix(V),
            as.matrix(beta),
            as.matrix(alpha),
            as.integer(maxiter),
            as.integer(nboot),
            as.double(maxtime),
            as.double(ll_threshold),
            as.integer(maxnewton),
            as.logical(verbose),
            as.character(adjust_prob_cond_method),
            as.logical(adjust_prob_cond_every)
        )

        object$sd_beta <- result$sd_beta
        if (!is.null(colnames(W))) {
            rownames(object$sd_beta) <- colnames(W)
        }
        if (!is.null(colnames(object$X))) {
            colnames(object$sd_beta) <- colnames(object$X)[-ncol(object$X)]
        }

        object$sd_alpha <- result$sd_alpha
        if (!is.null(colnames(object$X))) {
            rownames(object$sd_alpha) <- colnames(object$X)[-ncol(object$X)]
        }
        if (!is.null(colnames(V))) {
            colnames(object$sd_alpha) <- colnames(V)
        }

        object$nboot <- nboot
        class(object) <- "eim"
        return(object)
    }

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

    object$sd <- result$sd
    dimnames(object$sd) <- list(colnames(W), colnames(object$X))
    object$avg_prob <- result$avg_prob
    dimnames(object$avg_prob) <- list(colnames(W), colnames(object$X))
    object$nboot <- nboot
    object$sd[object$sd == 9999] <- Inf
    object$avg_prob[object$avg_prob == 9999] <- Inf

    class(object) <- "eim"
    return(object)
}

#' Reduce Parametric Covariates with PCA
#'
#' @description
#' Applies a Principal Component Analysis (PCA) to the covariates matrix `V` and
#' replaces it with a lower dimensional representation. This function is intended
#' for parametric workflows and requires a valid `V` matrix.
#'
#' @param object An object of class `eim`, which can be created using the [eim] function.
#'
#' @param X A `(b x c)` matrix representing candidate votes per ballot box.
#'
#' @param W A `(b x g)` matrix representing group votes per ballot box.
#'
#' @param V A `(b x a)` matrix with parametric covariates.
#'
#' @param json_path A path to a JSON file containing `X`, `W`, and `V` fields.
#'
#' @param components Integer specifying the number of principal components to keep.
#'
#' @param sd_threshold Numeric in `(0, 1]` indicating the minimum cumulative proportion
#'   of variance explained by the retained components.
#'
#' @param center Logical indicating whether to center the columns of `V` before PCA.
#'
#' @param scale Logical indicating whether to scale the columns of `V` before PCA.
#'
#' @return
#' Returns an `eim` object with the `V` matrix replaced by its PCA scores. The
#' columns of `V` are renamed as `PCA 1`, `PCA 2`, ..., up to the chosen number
#' of components.
#'
#' @examples
#' sim <- simulate_election(
#'     num_ballots = 50,
#'     num_candidates = 3,
#'     num_groups = 2,
#'     ballot_voters = 40,
#'     num_covariates = 10,
#'     num_districts = 2,
#'     seed = 1
#' )
#'
#' sim_pca <- PCA(sim, components = 2)
#' sim_pca$V
#'
#' @name PCA
#' @aliases PCA()
#' @export
PCA <- function(object = NULL,
                X = NULL,
                W = NULL,
                V = NULL,
                json_path = NULL,
                components = NULL,
                sd_threshold = NULL,
                center = TRUE,
                scale = TRUE) {
    all_params <- lapply(as.list(match.call(expand.dots = TRUE)), eval, parent.frame())
    .validate_compute(all_params) # nolint

    if (is.null(object)) {
        object <- eim(X = X, W = W, V = V, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("PCA: The object must be initialized with the `eim()` function.")
    }

    if (!is.null(V)) {
        object$V <- V
    }

    if (is.null(object$V)) {
        stop("PCA: This function requires a parametric object with a V matrix.")
    }

    if (!is.null(components) && !is.null(sd_threshold)) {
        stop("PCA: Provide either 'components' or 'sd_threshold', not both.")
    }

    if (is.null(components) && is.null(sd_threshold)) {
        stop("PCA: Provide either 'components' or 'sd_threshold'.")
    }

    Vmat <- as.matrix(object$V)
    if (ncol(Vmat) < 1) {
        stop("PCA: 'V' must have at least 1 column.")
    }
    if (any(!is.finite(Vmat))) {
        stop("PCA: 'V' cannot contain missing values or infinite values.")
    }

    pca <- stats::prcomp(Vmat, center = center, scale. = scale)

    if (!is.null(sd_threshold)) {
        if (!is.numeric(sd_threshold) || length(sd_threshold) != 1 || is.na(sd_threshold) ||
            sd_threshold <= 0 || sd_threshold > 1) {
            stop("PCA: 'sd_threshold' must be a numeric value in (0, 1].")
        }
        var_ratio <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
        components <- which(var_ratio >= sd_threshold)[1]
    }

    if (!is.numeric(components) || length(components) != 1 || is.na(components) ||
        components < 1 || as.integer(components) != components) {
        stop("PCA: 'components' must be a positive integer.")
    }

    components <- as.integer(components)
    if (components > ncol(Vmat)) {
        stop("PCA: 'components' cannot exceed the number of columns in V.")
    }

    scores <- pca$x[, seq_len(components), drop = FALSE]
    colnames(scores) <- paste("PCA", seq_len(components))
    rownames(scores) <- rownames(Vmat)

    object$V <- scores
    class(object) <- "eim"
    return(object)
}

#' Runs the EM algorithm aggregating adjacent groups, maximizing the variability of macro-group allocation in ballot boxes.
#'
#' This function estimates the voting probabilities (computed using [run_em]) aggregating adjacent groups so that the estimated probabilities' standard deviation (computed using [bootstrap]) is below a given threshold. See **Details** for more information.
#'
#' @note
#' This function only supports non-parametric models. Parametric objects (with `V`) are not supported.
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
        object <- eim(X = X, W = W, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("Bootstrap: The object must be initialized with the `eim()` function.")
    }
    if (!is.null(object$V)) {
        stop("get_agg_proxy: Parametric models are not supported.")
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
#' @note
#' This function only supports non-parametric models. Parametric objects (with `V`) are not supported.
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
        object <- eim(X = X, W = W, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("get_agg_opt: The object must be initialized with the `eim()` function.")
    }
    if (!is.null(object$V)) {
        stop("get_agg_opt: Parametric models are not supported.")
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
#' @note
#' This function does not support covariate models (i.e., `eim` objects with non-NULL `V` attributes).
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

    if (!is.null(object$V) || !is.null(object2$V)) {
        stop("waldtest: Parametric models are not supported.")
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


#' Save an `eim` object to a file
#'
#' This function saves an `eim` object to a specified file format. Supported formats are
#' **RDS**, **JSON**, and **CSV**. The function dynamically extracts and saves all available
#' attributes when exporting to JSON. If the `prob` field exists, it is saved when using CSV;
#' otherwise, it yields an error.
#'
#' @note
#' This function supports both non-parametric and parametric models. For parametric probabilities, the CSV output is a flattened matrix where rows correspond to ballot-box and group pairs.
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
#'   - For parametric probabilities, the 3D array is flattened into a 2D matrix with rows for each ballot-box/group pair.
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

            # if it's our 3-D array cond_prob, swap dim 1 <-> dim 3
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
            if (is.array(object$prob) && length(dim(object$prob)) == 3) {
                B <- dim(object$prob)[3]
                G <- dim(object$prob)[1]
                C <- dim(object$prob)[2]
                p_bgc <- aperm(object$prob, c(3, 1, 2)) # B x G x C
                flat <- matrix(0, nrow = B * G, ncol = C)
                row_labels <- character(B * G)
                idx <- 1
                ballot_names <- if (!is.null(rownames(object$X))) rownames(object$X) else paste0("B", seq_len(B))
                group_names <- if (!is.null(colnames(object$W))) colnames(object$W) else paste0("G", seq_len(G))
                for (b in seq_len(B)) {
                    for (g in seq_len(G)) {
                        flat[idx, ] <- p_bgc[b, g, ]
                        row_labels[idx] <- paste0(ballot_names[b], "::", group_names[g])
                        idx <- idx + 1
                    }
                }
                colnames(flat) <- colnames(object$X)
                rownames(flat) <- row_labels
                write.csv(flat, filename, row.names = TRUE)
                message("Probability matrix saved as CSV: ", filename)
                return(invisible(NULL))
            }
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


# Maybe for a future patch, if it's needed, add the option to get the eaxct log-likelihood
# exactLL <- function(object, scale_factor = 1) {
#     if (is.null(object$X) || is.null(object$W) || is.null(object$prob)) {
#         stop("The object must contain X, W and prob matrices.")
#     }
#
#     # Applies a scaling
#     if (scale_factor != 1) {
#         object$X <- round(object$X / scale_factor)
#         object$W <- round(object$W / scale_factor)
#         object$W <- .dhondt_correction(object$W, object$X)
#     }
#
#     ll <- computeExactLL(t(object$X), object$W, object$prob)
#     ll
# }
