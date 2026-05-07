#' Internal function!
#' Validates all of the 'run_em' arguments
#' @noRd
.validate_compute <- function(args) {
    # General checks: Vectors aren't accepted
    # if (any(sapply(args, function(x) length(x) > 1))) {
    #    stop("run_em:\tInvalid input: no vector inputs allowed")
    # }

    object_provided <- "object" %in% names(args) || "object1" %in% names(args)
    x_provided <- "X" %in% names(args) || "X1" %in% names(args)
    w_provided <- "W" %in% names(args) || "W1" %in% names(args)
    xw_provided <- x_provided || w_provided
    json_provided <- "json_path" %in% names(args)

    if (x_provided + w_provided == 1) {
        stop("If providing a matrix, 'X' and 'W' must be provided.")
    }

    if (sum(object_provided, xw_provided, json_provided) != 1) {
        stop(
            "You must provide exactly one of the following:\n",
            "(1)\tan `eim` object (initialized with `eim`)\n",
            "(2)\t`X` and `W`\n",
            "(3)\ta `json_path`"
        )
    }

    # Mismatch argument
    if ("allow_mismatch" %in% names(args)) {
        if (!is.logical(args$allow_mismatch)) {
            stop("run_em: Invalid 'allow_mismatch'. It has to be a boolean")
        }
    }

    # Method argument
    valid_methods <- c("mcmc", "exact", "mvn_cdf", "mvn_pdf", "mult", "metropolis")
    if ("method" %in% names(args) &&
        (!is.character(args$method) || length(args$method) != 1 || !(args$method %in% valid_methods))) {
        stop("Invalid 'method'. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    valid_symmetric_weight_methods <- c("average", "delta_ll", "mae_inverse", "joint")
    if ("symmetric_weight_method" %in% names(args) &&
        (!is.character(args$symmetric_weight_method) ||
            length(args$symmetric_weight_method) != 1 ||
            !(args$symmetric_weight_method %in% valid_symmetric_weight_methods))) {
        stop(
            "Invalid 'symmetric_weight_method'. Must be one of: ",
            paste(valid_symmetric_weight_methods, collapse = ", ")
        )
    }

    # Initial prob argument
    # valid_p_methods <- c("group_proportional", "proportional", "uniform", "random", "mult", "mcmc", "mvn_cdf", "mvn_pdf", "exact")
    # if ("initial_prob" %in% names(args) && (!is.matrix(args$initial_prob) ||
    #     (!is.character(args$initial_prob) || length(args$initial_prob) != 1 || !(args$initial_prob %in% valid_p_methods)))) {
    #     stop("Invalid 'initial_prob'. Must be one of: ", paste(valid_p_methods, collapse = ", "))
    # }

    if ("maxiter" %in% names(args)) {
        if (!is.numeric(args$maxiter) || as.integer(args$maxiter) != args$maxiter || args$maxiter < 1) { # Infinite are valid, skip this case
            stop("Invalid 'maxiter'. Must be a positive integer.")
        }
    }

    if ("ll_threshold" %in% names(args)) {
        if (!is.infinite(args$ll_threshold) && (!is.numeric(args$ll_threshold) || args$ll_threshold < 0)) { # Infinite are valid, skip this case
            stop("Invalid 'll_threshold'. Must be a positive numeric or infinite value.")
        }
    }

    # Maxtime argument
    if ("maxtime" %in% names(args) &&
        (!is.numeric(args$maxtime) || args$maxtime < 0)) {
        stop("Invalid 'maxtime'. Must be positive.")
    }

    # Stop threshold argument
    if ("param_threshold" %in% names(args)) {
        if (!is.infinite(args$param_threshold) && (!is.numeric(args$param_threshold) || args$param_threshold < 0)) {
            stop("run_em: Invalid 'param_threshold'. Must be a positive numeric or an infinite value.")
        }
        if (args$param_threshold >= 1) {
            warning("Warning: A 'param_threshold' greater or equal than one will always be true after the first iteration.")
        }
        if ("compute_ll" %in% names(args) && !args$compute_ll && is.infinite(args$param_threshold)) {
            stop("You must provide a parameter threshold if 'compute_ll' is FALSE.")
        }
    }

    # Verbose argument
    if ("verbose" %in% names(args) && !is.logical(args$verbose)) {
        stop("run_em: Invalid 'verbose'. It has to be a boolean.")
    }

    # mcmc: mcmc_stepsize argument
    if ("mcmc_stepsize" %in% names(args)) {
        if (!is.numeric(args$mcmc_stepsize) || as.integer(args$mcmc_stepsize) != args$mcmc_stepsize || args$mcmc_stepsize < 0) {
            stop("run_em: Invalid 'mcmc_stepsize'. Must be a positive integer.")
        }
        if (args$mcmc_stepsize < 15) {
            warning("Warning: A small 'mcmc_stepsize' could lead to highly correlated samples.")
        }
    }

    # mcmc: Samples argument
    if ("samples" %in% names(args) &&
        (!is.numeric(args$mcmc_samples) || as.integer(args$mcmc_samples) != args$mcmc_samples || args$mcmc_samples < 0)) {
        stop("run_em: Invalid 'mcmc_samples'. Must be a positive integer.")
    }

    # metropolis: Samples argument
    if ("metropolis_iter" %in% names(args) &&
        (!is.numeric(args$metropolis_iter) || as.integer(args$metropolis_iter) != args$metropolis_iter || args$metropolis_iter < 0)) {
        stop("run_em: Invalid 'metropolis_iter'. Must be a positive integer.")
    }

    # CDF: Mc_method argument
    valid_cdf_methods <- c("genz", "genz2")
    if ("mvncdf_method
" %in% names(args) &&
        (!is.character(args$mvncdf_method) || !args$mvncdf_method
            %in% valid_cdf_methods)) {
        stop("run_em: Invalid 'mvncdf_method
'. Must be one of: ", paste(valid_cdf_methods, collapse = ", "))
    }

    # CDF: Mc_error argument
    if ("mvncdf_error" %in% names(args) &&
        (!is.numeric(args$mvncdf_error) || args$mvncdf_error <= 0)) {
        stop("run_em: Invalid 'mvncdf_error'. Must be a positive number.")
    }

    # CDF: Mc_error argument
    if ("mvncdf_samples" %in% names(args) &&
        (!is.numeric(args$mvncdf_samples) || as.integer(args$mvncdf_samples) != args$mvncdf_samples || args$mvncdf_samples < 0)) {
        stop("run_em: Invalid 'mvncdf_samples'. Must be a positive integer.")
    }

    # Check mismatch
    if ("mismatch" %in% names(args)) {
        if (!is.logical(args$mismatch)) {
            stop("run_em: Invalid 'mismatch'. Must be a boolean value.")
        }
        # if ("method" %in% names(args) && "method" %in% c("exact")) {
        #    stop("run_em: Mismatched results are not supported when using 'exact'.")
        # }
    }

    # Include nboot aswell if bootstrapping is provided
    if ("nboot" %in% names(args) &&
        (!is.numeric(args$nboot) || as.integer(args$nboot) != args$nboot || args$nboot < 0)) {
        stop("Bootstrap: Invalid 'nboot'. Must be a positive integer.")
    }

    valid_sd_methods <- c("maximum", "average")
    if ("sd_statistic" %in% names(args) &&
        (!is.character(args$sd_statistic) || length(args$sd_statistic) != 1 || !(args$sd_statistic %in% valid_sd_methods))) {
        stop("Invalid 'sd_statistic'. Must be one of: ", paste(valid_sd_methods, collapse = ", "))
    }

    if ("sd_threshold" %in% names(args) &&
        (!is.numeric(args$sd_threshold) || args$sd_threshold <= 0)) {
        stop("Invalid 'sd_threshold'. Must be a positive number.")
    }

    if ("alternative" %in% names(args) && !args$alternative %in% c("two.sided", "greater", "less")) {
        stop("Invalid 'alternative'. Must be one of: two.sided, greater, less")
    }

    valid_lp_methods <- c("", "lp", "project_lp")
    if ("adjust_prob_cond_method" %in% names(args) &&
        (!is.character(args$adjust_prob_cond_method) || !(args$adjust_prob_cond_method %in% valid_lp_methods))) {
        stop("Invalid 'adjust_prob_cond_method'. Must be one of: ", paste(valid_lp_methods, collapse = ", "))
    }

    if ("adjust_prob_cond_every" %in% names(args)) {
        if (!is.logical(args$adjust_prob_cond_every)) {
            stop("Invalid 'adjust_prob_cond_every'. Must be a boolean value.")
        }
        if ("adjust_prob_cond_method" %in% names(args) && args$adjust_prob_cond_method == "") {
            warning("You provided 'adjust_prob_cond_every' but not 'adjust_prob_cond_method'. The former will be ignored.")
        }
    }
}


#' Internal function!
#'
#' Validate the 'eim' object inputs
#'
#' @param X A matrix representing candidate votes per ballot box.
#' @param W A matrix representing group votes per ballot box.
#' @return Stops execution if validation fails.
#' @noRd
.validate_eim <- function(X, W) {
    # Ensure X and W are provided
    if (is.null(X) || is.null(W)) {
        stop("Either provide X and W matrices, or a valid JSON path containing them.")
    }

    if (!is.matrix(X) || !is.matrix(W)) {
        stop("'X' and 'W' must be matrices.")
    }

    # Ensure they are matrices
    X <- as.matrix(X)
    W <- as.matrix(W)

    # Check matching dimensions
    if (nrow(X) != nrow(W)) {
        stop(
            "Mismatch in the number of ballot boxes: 'X' has ", nrow(X),
            " rows, but 'W' has ", nrow(W), " rows."
        )
    }

    # Check minimum column constraints
    if (ncol(X) < 2) {
        stop("Candidate matrix 'X' must have at least 2 columns.")
    }
    if (ncol(W) < 2) {
        # stop("Group matrix 'W' must have at least 2 columns.")
    }

    # Check for missing values
    if (any(is.na(X)) || any(is.na(W))) {
        stop("Matrices 'X' and 'W' cannot contain missing values (NA).")
    }

    TRUE
}

#' Internal function!
#'
#' Validate the 'eim' object JSON path
#'
#' @param json_path A path to a JSON file containing `"X"` and `"W"`.
#' @return A list with the `"X"` and `"W"` matrices. Stops execution if validation fails.
#' @noRd
.validate_json_eim <- function(json_path) {
    if (!file.exists(json_path)) {
        stop("The specified JSON file does not exist: ", json_path)
    }

    data <- tryCatch(
        jsonlite::fromJSON(json_path),
        error = function(e) stop("Failed to read JSON file: ", e$message)
    )

    # Validate JSON contents
    if (!all(c("X", "W") %in% names(data))) {
        stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")
    }

    if (is.null(data$X) || is.null(data$W)) {
        stop("'X' and 'W' cannot be NULL in the JSON file")
    }

    result <- list(
        X = as.matrix(data$X),
        W = as.matrix(data$W)
    )

    result
}

#' Internal function!
#'
#' Normalize each row of a probability matrix so every row sums to one.
#'
#' @param prob_matrix Numeric matrix with one probability vector per row.
#' @return A matrix with finite, non-negative rows normalized to one.
#' @noRd
.normalize_prob_rows <- function(prob_matrix) {
    prob_matrix[!is.finite(prob_matrix) | prob_matrix < 0] <- 0
    row_sums <- rowSums(prob_matrix)
    valid_rows <- is.finite(row_sums) & row_sums > 0

    if (any(valid_rows)) {
        prob_matrix[valid_rows, ] <- sweep(
            prob_matrix[valid_rows, , drop = FALSE],
            1,
            row_sums[valid_rows],
            "/"
        )
    }

    if (any(!valid_rows)) {
        prob_matrix[!valid_rows, ] <- rep(1 / ncol(prob_matrix), ncol(prob_matrix))
    }

    prob_matrix
}

#' Internal function!
#'
#' Normalize the last dimension of a 3d-array into valid probability vectors.
#'
#' @param arr3 Numeric 3d-array where the last dimension is normalized.
#' @return A 3d-array with finite, non-negative slices summing to one.
#' @noRd
.normalize_cube_last_dim <- function(arr3) {
    num_ballots <- dim(arr3)[1]
    num_rows <- dim(arr3)[2]
    num_cols <- dim(arr3)[3]

    for (ballot in seq_len(num_ballots)) {
        for (row in seq_len(num_rows)) {
            values <- arr3[ballot, row, ]
            values[!is.finite(values) | values < 0] <- 0
            values_sum <- sum(values)

            if (is.finite(values_sum) && values_sum > 0) {
                arr3[ballot, row, ] <- values / values_sum
            } else {
                arr3[ballot, row, ] <- rep(1 / num_cols, num_cols)
            }
        }
    }

    arr3
}

#' Internal function!
#'
#' Perform the M-step from ballot-box conditional probabilities.
#'
#' @param q_array A `(g x c x b)` array with conditional probabilities.
#' @param W_matrix A `(b x g)` matrix with group totals.
#' @return A `(g x c)` matrix with updated probabilities.
#' @noRd
.mstep_from_q <- function(q_array, W_matrix) {
    q_bgc <- aperm(q_array, c(3, 1, 2))
    weighted_q <- sweep(q_bgc, c(1, 2), W_matrix, "*")
    numerator <- apply(weighted_q, c(2, 3), sum)
    denominator <- colSums(W_matrix)
    probabilities <- sweep(numerator, 1, denominator, "/")

    .normalize_prob_rows(probabilities)
}

#' Internal function!
#'
#' Evaluate the log-likelihood of a fixed probability matrix under an EM method.
#'
#' @param object An `eim` object providing method-specific controls.
#' @param X_matrix A candidate-vote matrix.
#' @param W_matrix A group-vote matrix.
#' @param prob_matrix A probability matrix to evaluate.
#' @param method Character string with the EM method name.
#' @param miniter Minimum number of iterations stored in the temporary object.
#' @param adjust_prob_cond_method Character string controlling probability adjustment.
#' @param adjust_prob_cond_every Boolean indicating whether to project at every iteration.
#' @return A numeric scalar with the evaluated log-likelihood.
#' @noRd
.run_em_loglik_from_prob <- function(object,
                                     X_matrix,
                                     W_matrix,
                                     prob_matrix,
                                     method,
                                     miniter,
                                     adjust_prob_cond_method,
                                     adjust_prob_cond_every) {
    ll_object <- list(
        X = as.matrix(X_matrix),
        W = as.matrix(W_matrix),
        method = method,
        mcmc_stepsize = if (!is.null(object$mcmc_stepsize)) object$mcmc_stepsize else 3000,
        mcmc_samples = if (!is.null(object$mcmc_samples)) object$mcmc_samples else 1000,
        mvncdf_method = if (!is.null(object$mvncdf_method)) object$mvncdf_method else "genz",
        mvncdf_error = if (!is.null(object$mvncdf_error)) object$mvncdf_error else 1e-3,
        mvncdf_samples = if (!is.null(object$mvncdf_samples)) object$mvncdf_samples else 5000,
        miniter = miniter,
        adjust_prob_cond_method = adjust_prob_cond_method,
        adjust_prob_cond_every = adjust_prob_cond_every
    )
    class(ll_object) <- "eim"

    as.numeric(logLik(ll_object, prob = prob_matrix, method = method))
}

#' Internal function!
#'
#' Build the standard dimnames used by `run_em()` 3d outputs.
#'
#' @param object An `eim` object with candidate and ballot-box names.
#' @param W_matrix A group matrix providing group names.
#' @return A list of dimnames ordered as groups, candidates, ballot-boxes.
#' @noRd
.run_em_dimnames <- function(object, W_matrix) {
    list(
        colnames(W_matrix),
        colnames(object$X),
        rownames(object$X)
    )
}

#' Internal function!
#'
#' Retrieve the group matrix used internally by `run_em()`.
#'
#' @param object An `eim` object.
#' @return `object$W_agg` when available, otherwise `object$W`.
#' @noRd
.run_em_working_group_matrix <- function(object) {
    if (is.null(object$W_agg)) {
        object$W
    } else {
        object$W_agg
    }
}

#' Internal function!
#'
#' Apply a group aggregation specification to the group matrix of an `eim` object.
#'
#' @param object An `eim` object.
#' @param group_agg A vector describing how to aggregate `W` columns.
#' @return The updated `eim` object, including `W_agg` and `group_agg` when requested.
#' @noRd
.run_em_apply_group_agg <- function(object, group_agg) {
    if (is.null(group_agg)) {
        return(object)
    }

    sizes <- diff(c(0, group_agg))
    rep_labels <- rep(seq_along(sizes), sizes)
    groups <- split(seq_len(ncol(object$W)), rep_labels)
    W_agg <- do.call(
        cbind,
        lapply(groups, function(cols) rowSums(object$W[, cols, drop = FALSE]))
    )
    rownames(W_agg) <- rownames(object$W)

    object$W_agg <- W_agg
    object$group_agg <- group_agg
    object
}

#' Internal function!
#'
#' Prepare the input object used by `run_em()` before dispatching the EM routine.
#'
#' @param object Optional `eim` object provided by the user.
#' @param X Candidate-vote matrix.
#' @param W Group-vote matrix.
#' @param json_path Optional JSON path used to build the object.
#' @param scale_factor Numeric scaling factor applied to `X` and `W`.
#' @param method Character string with the selected EM method.
#' @param allow_mismatch Boolean indicating whether row mismatches are allowed.
#' @param group_agg Optional aggregation specification for `W`.
#' @return A prepared `eim` object ready to run EM.
#' @noRd
.run_em_prepare_object <- function(object,
                                   X,
                                   W,
                                   json_path,
                                   scale_factor,
                                   method,
                                   allow_mismatch,
                                   group_agg) {
    if (is.null(object)) {
        object <- eim(X = X, W = W, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("run_em: The object must be initialized with the eim() function.")
    }

    if (scale_factor != 1) {
        object$X <- round(object$X / scale_factor)
        object$W <- round(object$W / scale_factor)
    }

    mismatch_rows <- which(rowSums(object$X) != rowSums(object$W))
    if (!allow_mismatch && length(mismatch_rows) > 0) {
        stop(
            "run_em: Row-wise mismatch in vote totals detected.\n",
            "Rows with mismatches: ", paste(mismatch_rows, collapse = ", "), "\n",
            "To allow mismatches, set `allow_mismatch = TRUE`."
        )
    }

    if (method == "exact" && length(mismatch_rows) > 0) {
        .dhondt_correction(object$W, object$X)
        message("Applying a D'Hondt correction for correcting mismatches in W")
    }

    object <- .run_em_apply_group_agg(object, group_agg)
    object$method <- method
    object
}

#' Internal function!
#'
#' Populate method-specific defaults for EM runs.
#'
#' @param object An `eim` object.
#' @param method Character string with the selected EM method.
#' @param all_params List of evaluated `run_em()` arguments.
#' @return The updated `eim` object with method-specific defaults.
#' @noRd
.run_em_apply_method_defaults <- function(object, method, all_params) {
    if (method == "mcmc") {
        object$mcmc_stepsize <- as.integer(
            if ("mcmc_stepsize" %in% names(all_params)) all_params$mcmc_stepsize else 3000
        )
        object$mcmc_samples <- as.integer(
            if ("mcmc_samples" %in% names(all_params)) all_params$mcmc_samples else 1000
        )
        object$burn_in <- as.integer(
            if ("burn_in" %in% names(all_params)) all_params$burn_in else 10000
        )
    } else if (method == "mvn_cdf") {
        object$mvncdf_method <- if ("mvncdf_method" %in% names(all_params)) all_params$mvncdf_method else "genz"
        object$mvncdf_samples <- if ("mvncdf_samples" %in% names(all_params)) all_params$mvncdf_samples else 5000
        object$mvncdf_error <- if ("mvncdf_error" %in% names(all_params)) all_params$mvncdf_error else 1e-3
    }

    object
}

#' Internal function!
#'
#' Build the recursive call used for the reverse symmetric run.
#'
#' @param base_call Original `run_em()` call.
#' @param object The forward-run `eim` object.
#' @param all_params List of evaluated `run_em()` arguments.
#' @return A modified call object for the reverse run.
#' @noRd
.run_em_inverse_call <- function(base_call, object, all_params) {
    base_call_sym <- base_call
    base_call_sym$symmetric <- FALSE
    base_call_sym$X <- object$W
    base_call_sym$W <- object$X
    base_call_sym$json_path <- NULL
    base_call_sym$object <- NULL
    base_call_sym$scale_factor <- 1

    initial_prob <- all_params$initial_prob
    if (is.matrix(initial_prob)) {
        col_totals_x <- colSums(object$X)
        numerator <- sweep(initial_prob, 2, col_totals_x, "*")
        denominator <- rowSums(numerator)
        base_call_sym$initial_prob <- sweep(numerator, 1, denominator, "/")
        base_call_sym$initial_prob <- t(base_call_sym$initial_prob)
    }

    base_call_sym
}

#' Internal function!
#'
#' Copy EM results back into an `eim` object.
#'
#' @param object An `eim` object.
#' @param resulting_values Output list returned by `EMAlgorithmFull`.
#' @param W_matrix Group matrix used in the fit.
#' @param control List with the active `run_em()` controls.
#' @return The updated `eim` object.
#' @noRd
.run_em_assign_results <- function(object, resulting_values, W_matrix, control) {
    object$cond_prob <- resulting_values$q
    dimnames(object$cond_prob) <- .run_em_dimnames(object, W_matrix)
    object$expected_outcome <- resulting_values$expected_outcome
    dimnames(object$expected_outcome) <- .run_em_dimnames(object, W_matrix)
    object$prob <- as.matrix(resulting_values$result)
    dimnames(object$prob) <- list(colnames(W_matrix), colnames(object$X))
    object$iterations <- as.numeric(resulting_values$total_iterations)

    if (control$compute_ll) {
        object$logLik <- as.numeric(resulting_values$log_likelihood[length(resulting_values$log_likelihood)])
    }

    object$time <- resulting_values$total_time
    object$message <- resulting_values$stopping_reason
    object$status <- as.integer(resulting_values$finish_id)
    object$miniter <- control$miniter
    object$maxiter <- control$maxiter
    object$maxtime <- control$maxtime
    object$param_threshold <- control$param_threshold
    object$ll_threshold <- control$ll_threshold
    object$initial_prob <- control$initial_prob
    object$adjust_prob_cond_method <- control$adjust_prob_cond_method
    object$adjust_prob_cond_every <- control$adjust_prob_cond_every

    object
}

#' Internal function!
#'
#' Finalize the direct joint symmetric path without storing inverse outputs.
#'
#' @param object An `eim` object returned by the joint symmetric run.
#' @return The updated `eim` object with joint symmetric metadata.
#' @noRd
.run_em_finalize_joint <- function(object) {
    object$cond_prob_inv <- NULL
    object$prob_inv <- NULL
    object$expected_outcome_inv <- NULL
    object$symmetric_weight_method <- "joint"
    object$symmetric_weights <- c(original = 0.5, reverse = 0.5)

    object
}

#' Internal function!
#'
#' Compute the weights used to combine forward and reverse symmetric runs.
#'
#' @param object The forward-run `eim` object.
#' @param inverse The reverse-run `eim` object.
#' @param W_sym Group matrix used to rebuild symmetric probabilities.
#' @param control List with the active `run_em()` controls.
#' @return A list with the updated object and the two symmetric weights.
#' @noRd
.run_em_symmetric_weights <- function(object, inverse, W_sym, control) {
    weight_original <- 0.5
    weight_reverse <- 0.5

    if (identical(control$symmetric_weight_method, "delta_ll")) {
        forward_ll <- suppressWarnings(as.numeric(object$logLik))
        reverse_ll <- suppressWarnings(as.numeric(inverse$logLik))

        if (is.finite(forward_ll) && is.finite(reverse_ll)) {
            q_orig_bgc <- aperm(object$cond_prob, c(3, 1, 2))
            z_from_orig_bgc <- sweep(q_orig_bgc, c(1, 2), W_sym, "*")
            q_rev_ind_bcg <- sweep(aperm(z_from_orig_bgc, c(1, 3, 2)), c(1, 2), object$X, "/")
            q_rev_ind <- aperm(.normalize_cube_last_dim(q_rev_ind_bcg), c(2, 3, 1))
            p_rev_ind <- .mstep_from_q(q_rev_ind, object$X)

            q_rev_bcg <- aperm(inverse$cond_prob, c(3, 1, 2))
            z_from_rev_bgc <- aperm(sweep(q_rev_bcg, c(1, 2), object$X, "*"), c(1, 3, 2))
            q_ind_bgc <- sweep(z_from_rev_bgc, c(1, 2), W_sym, "/")
            q_ind <- aperm(.normalize_cube_last_dim(q_ind_bgc), c(2, 3, 1))
            p_ind <- .mstep_from_q(q_ind, W_sym)

            LL_ind <- .run_em_loglik_from_prob(
                object,
                object$X,
                W_sym,
                p_ind,
                control$method,
                control$miniter,
                control$adjust_prob_cond_method,
                control$adjust_prob_cond_every
            )
            LL_rev_ind <- .run_em_loglik_from_prob(
                object,
                object$W,
                object$X,
                p_rev_ind,
                control$method,
                control$miniter,
                control$adjust_prob_cond_method,
                control$adjust_prob_cond_every
            )
            dLL <- forward_ll - LL_ind
            dLL_rev <- reverse_ll - LL_rev_ind

            tau <- if (abs(forward_ll) > .Machine$double.eps) max(0, dLL / abs(forward_ll)) else 0
            tau_rev <- if (abs(reverse_ll) > .Machine$double.eps) max(0, dLL_rev / abs(reverse_ll)) else 0
            denominator <- tau + tau_rev

            if (is.finite(denominator) && denominator > 0) {
                weight_original <- tau_rev / denominator
                weight_reverse <- tau / denominator
            }

            object$LL_ind <- LL_ind
            object$LL_rev_ind <- LL_rev_ind
            object$dLL <- dLL
            object$dLL_rev <- dLL_rev
            object$nu <- tau
            object$nu_rev <- tau_rev
            object$symmetric_weight_method <- "delta_ll"
        } else {
            object$symmetric_weight_method <- "average"
        }
    } else if (identical(control$symmetric_weight_method, "mae_inverse")) {
        prob_forward <- .mstep_from_q(object$cond_prob, W_sym)
        prob_reverse <- as.matrix(inverse$prob)

        if (!all(dim(prob_reverse) == c(ncol(object$X), ncol(W_sym)))) {
            prob_reverse <- .mstep_from_q(inverse$cond_prob, object$X)
        }

        if (all(dim(prob_reverse) == c(ncol(object$X), ncol(W_sym)))) {
            x_hat_forward <- W_sym %*% prob_forward
            w_hat_reverse <- object$X %*% prob_reverse
            tau <- sum(abs(object$X - x_hat_forward))
            tau_rev <- sum(abs(W_sym - w_hat_reverse))
            denominator <- tau + tau_rev

            if (is.finite(denominator) && denominator > 0) {
                weight_original <- tau_rev / denominator
                weight_reverse <- 1 - weight_original
            }

            object$err_forward <- tau
            object$err_inverse <- tau_rev
            object$symmetric_weight_method <- "mae_inverse"
        } else {
            object$symmetric_weight_method <- "average"
        }
    } else {
        object$symmetric_weight_method <- "average"
    }

    list(
        object = object,
        weight_original = weight_original,
        weight_reverse = weight_reverse
    )
}

#' Internal function!
#'
#' Combine the forward and reverse runs into a symmetric result.
#'
#' @param object The forward-run `eim` object.
#' @param inverse The reverse-run `eim` object.
#' @param control List with the active `run_em()` controls.
#' @return The symmetric `eim` object.
#' @noRd
.run_em_apply_symmetry <- function(object, inverse, control) {
    object$cond_prob_inv <- inverse$cond_prob
    object$prob_inv <- inverse$prob
    object$expected_outcome_inv <- inverse$expected_outcome
    object$time <- object$time + inverse$time
    object$iterations <- object$iterations + inverse$iterations

    W_sym <- .run_em_working_group_matrix(object)
    reverse_expected_outcome <- aperm(inverse$expected_outcome, c(2, 1, 3))
    weights <- .run_em_symmetric_weights(object, inverse, W_sym, control)
    object <- weights$object
    object$symmetric_weights <- c(
        original = weights$weight_original,
        reverse = weights$weight_reverse
    )
    object$expected_outcome <- weights$weight_original * object$expected_outcome +
        weights$weight_reverse * reverse_expected_outcome

    expected_bgc <- aperm(object$expected_outcome, c(3, 1, 2))
    cond_prob_bgc <- sweep(expected_bgc, c(1, 2), W_sym, "/")
    object$cond_prob <- aperm(.normalize_cube_last_dim(cond_prob_bgc), c(2, 3, 1))

    numerator <- apply(object$expected_outcome, c(1, 2), sum)
    denominator <- colSums(W_sym)
    object$prob <- sweep(numerator, 1, denominator, "/")
    object$prob <- .normalize_prob_rows(object$prob)

    object
}

#' Internal function!
#'
#' Execute `run_em()`.
#'
#' @param object A prepared `eim` object.
#' @param control List with the active `run_em()` controls.
#' @return An updated `eim` object with EM results.
#' @noRd
.run_em_core <- function(object, control) {
    object <- .run_em_apply_method_defaults(object, control$method, control$all_params)
    W_matrix <- .run_em_working_group_matrix(object)

    resulting_values <- EMAlgorithmFull(
        t(object$X),
        W_matrix,
        control$method,
        if (is.character(control$initial_prob)) control$initial_prob else "custom",
        control$maxiter,
        control$maxtime,
        control$param_threshold,
        control$ll_threshold,
        control$compute_ll,
        control$verbose,
        as.integer(if (!is.null(object$mcmc_stepsize)) object$mcmc_stepsize else 3000),
        as.integer(if (!is.null(object$mcmc_samples)) object$mcmc_samples else 1000),
        if (!is.null(object$mvncdf_method)) object$mvncdf_method else "genz",
        as.numeric(if (!is.null(object$mvncdf_error)) object$mvncdf_error else 1e-3),
        as.numeric(if (!is.null(object$mvncdf_samples)) object$mvncdf_samples else 5000),
        control$miniter,
        control$adjust_prob_cond_method,
        control$adjust_prob_cond_every,
        if (is.matrix(control$initial_prob)) control$initial_prob else matrix(-1, nrow = 1, ncol = 1),
        control$symmetric,
        control$symmetric_weight_method
    )

    object <- .run_em_assign_results(object, resulting_values, W_matrix, control)
    if (control$symmetric && identical(control$symmetric_weight_method, "joint")) {
        return(.run_em_finalize_joint(object))
    }
    if (!control$symmetric) {
        return(object)
    }

    inverse_call <- .run_em_inverse_call(control$base_call, object, control$all_params)
    inverse <- eval(inverse_call, control$caller_env)

    .run_em_apply_symmetry(object, inverse, control)
}

#' Internal function!
#'
#' Randomly create a voting instance by defining an interval
#'
#' @description
#' Given a range of possible \strong{observed} outcomes (such as ballot boxes, number of candidates, etc.),
#' it creates a completely random voting instance, simulating the unobserved results as well.
#'
#' @param ballots_range (integer) A vector of size 2 with the lower and upper bound of ballot boxes.
#'
#' @param candidates_range (integer) A vector of size 2 with the lower and upper bound of candidates to draw.
#'
#' @param demographic_range (integer) A vector of size 2 with the lower and upper bound of demographic groups
#' to draw.
#'
#' @param voting_range (integer) A vector of size 2 with the lower and upper bound of votes per ballot box.
#'
#' @param seed \emph{(numeric(1)} Optional. If provided, it overrides the current global seed. (default: \code{NULL})
#'
#' @return A list with components:
#' \item{X}{A matrix (b x c) with candidate votes per ballot box.}
#' \item{W}{A matrix (b x g) with demographic votes per ballot box.}
#' \item{real_p}{A matrix (g x c) with the estimated \strong{(unobserved)} probabilities that a demographic group votes for a given candidate.}
#' \item{ballots}{The number of ballot boxes that were drawn.}
#' \item{candidates}{The number of candidates that were drawn.}
#' \item{groups}{The number of demographic groups that were drawn.}
#' \item{total_votes}{A vector with the number of total votes per ballot box.}
#'
#' @seealso [simulate_election()]
#' @examples
#'
#' bal_range <- c(30, 50)
#' can_range <- c(2, 4)
#' group_range <- c(2, 6)
#' voting_range <- c(50, 100)
#' results <- random_samples(bal_range, can_range, group_range, voting_range)
#'
#' # X matrix
#' results$X # A randomly generated matrix of dimension (b x c)
#' ncol(results$X <= can_range[2]) # Always TRUE
#' ncol(results$X >= can_range[1]) # Always TRUE
#' nrow(results$X <= bal_range[2]) # Always TRUE
#' nrow(results$X >= bal_range[1]) # Always TRUE
#'
#' # W matrix
#' results$W # A randomly generated matrix of dimension (b x g)
#' ncol(results$W <= group_range[2]) # Always TRUE
#' ncol(results$W >= group_range[1]) # Always TRUE
#' nrow(results$W <= bal_range[2]) # Always TRUE
#' nrow(results$W >= bal_range[1]) # Always TRUE
#'
#' # Probability matrix
#' results$real_p # A matrix (g x c) that summarizes the unobserved outcomes
#' ncol(results$real_p) == ncol(results$X) # Always TRUE
#' nrow(results$real_p) == ncol(results$W) # Always TRUE
#'
#' @noRd
.random_samples <- function(ballots_range, # Arguments must be vectors of size 2
                            candidates_range,
                            demographic_range,
                            voting_range,
                            seed = NULL) {
    param_list <- list(ballots_range, candidates_range, demographic_range, voting_range)
    if (!(all(sapply(param_list, length) == 2))) {
        stop("The vectors must be of size 2.")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # Randomly choose a ballot box
    num_ballots <- sample(ballots_range[1]:ballots_range[2], 1)
    # Randomly choose demographic groups
    num_groups <- sample(demographic_range[1]:demographic_range[2], 1)
    # Randomly choose candidates
    num_candidates <- sample(candidates_range[1]:candidates_range[2], 1)
    # Randomly choose the total amount of votes per ballot box
    total_votes <- sample(
        seq.int(voting_range[1], voting_range[2]),
        size = num_ballots,
        replace = TRUE
    )
    # Randomly choose the group proportions
    group_prop <- rgamma(num_groups, shape = 1, rate = 1)
    group_prop <- group_prop / sum(group_prop)

    choosen_values <- list(
        ballots = num_ballots,
        candidates = num_candidates,
        groups = num_groups,
        total_votes = total_votes
    )

    result <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = total_votes,
        seed = seed,
        group_proportions = group_prop
    )

    appended_list <- c(result, choosen_values)
    appended_list
}

#' Internal function!
#' Applies dhond't correction to W and X matrices
#' @noRd
.dhondt_correction <- function(W, X) {
    if (any(W < 0) || any(X < 0)) stop("W and X must be non-negative.")
    adjust_row <- function(w, target_sum) {
        w <- as.numeric(w)
        cur <- sum(w)

        # trivial cases
        if (target_sum == cur) {
            return(w)
        }

        # if row is all zeros and target > 0, start from ones via uniform weights
        if (cur == 0 && target_sum > 0) {
            seats <- integer(length(w))
            for (k in seq_len(target_sum)) {
                # uniform weights -> all quotients equal; break ties by first
                j <- which.max(rep(1, length(w)) / (seats + 1))
                seats[j] <- seats[j] + 1L
            }
            return(seats)
        }

        if (target_sum > cur) {
            # add (target_sum - cur) units by D’Hondt on base weights w
            add <- target_sum - cur
            seats <- integer(length(w))
            for (k in seq_len(add)) {
                quot <- w / (seats + 1)
                j <- which.max(quot)
                seats[j] <- seats[j] + 1L
            }
            return(w + seats)
        } else {
            # remove (cur - target_sum) units greedily from the largest current entries
            rem <- cur - target_sum
            out <- w
            for (k in seq_len(rem)) {
                j <- which.max(out)
                if (out[j] > 0) out[j] <- out[j] - 1L else break
            }
            return(out)
        }
    }

    targets <- rowSums(X)
    W_adj <- W
    for (i in seq_len(nrow(W))) {
        W_adj[i, ] <- adjust_row(W[i, ], targets[i])
    }

    W_adj
}
