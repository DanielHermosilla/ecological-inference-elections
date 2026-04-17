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
    v_provided <- "V" %in% names(args) || "V1" %in% names(args)
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

    if (v_provided && !is.null(args$V) && !is.matrix(args$V)) {
        stop("Invalid 'V'. It has to be a matrix.")
    }

    if ("alpha" %in% names(args) && !is.null(args$alpha) && !is.matrix(args$alpha)) {
        stop("Invalid 'alpha'. It has to be a matrix.")
    }

    if ("beta" %in% names(args) && !is.null(args$beta) && !is.matrix(args$beta)) {
        stop("Invalid 'beta'. It has to be a matrix.")
    }

    # Mismatch argument
    if ("allow_mismatch" %in% names(args)) {
        if (!is.logical(args$allow_mismatch)) {
            stop("run_em: Invalid 'allow_mismatch'. It has to be a boolean")
        }
    }

    # Method argument
    valid_methods <- c("mcmc", "exact", "mvn_cdf", "mvn_pdf", "saddlepoint", "mult", "metropolis")
    if ("method" %in% names(args) &&
        (!is.character(args$method) || length(args$method) != 1 || !(args$method %in% valid_methods))) {
        stop("Invalid 'method'. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    valid_symmetric_weight_methods <- c("average", "delta_ll", "EM_weight")
    if ("symmetric_weight_method" %in% names(args) &&
        (!is.character(args$symmetric_weight_method) ||
            length(args$symmetric_weight_method) != 1 ||
            !(args$symmetric_weight_method %in% valid_symmetric_weight_methods))) {
        stop(
            "Invalid 'symmetric_weight_method'. Must be one of: ",
            paste(valid_symmetric_weight_methods, collapse = ", ")
        )
    }

    if ("mixture" %in% names(args) && !is.null(args$mixture)) {
        if (!is.numeric(args$mixture) || length(args$mixture) != 1 ||
            !is.finite(args$mixture) || as.integer(args$mixture) != args$mixture ||
            args$mixture < 1) {
            stop("run_em: Invalid 'mixture'. Must be a positive integer.")
        }
    }

    if ("row_mixture" %in% names(args) && !is.null(args$row_mixture)) {
        if (!is.numeric(args$row_mixture) || length(args$row_mixture) != 1 ||
            !is.finite(args$row_mixture) || as.integer(args$row_mixture) != args$row_mixture ||
            args$row_mixture < 1) {
            stop("run_em: Invalid 'row_mixture'. Must be a positive integer.")
        }
    }

    if ("H" %in% names(args) && !is.null(args$H)) {
        if (!is.numeric(args$H) || length(args$H) != 1 ||
            !is.finite(args$H) || as.integer(args$H) != args$H ||
            args$H < 1) {
            stop("run_em: Invalid 'H'. Must be a positive integer.")
        }
    }

    if ("HET" %in% names(args) && !is.null(args$HET)) {
        if (!is.numeric(args$HET) || length(args$HET) != 1 ||
            !is.finite(args$HET) || args$HET < 0) {
            stop("run_em: Invalid 'HET'. Must be a non-negative numeric value.")
        }
    }

    if ("AE" %in% names(args) && !is.null(args$AE)) {
        if (!is.numeric(args$AE) || length(args$AE) != 1 ||
            !is.finite(args$AE) || args$AE < 0) {
            stop("run_em: Invalid 'AE'. Must be a non-negative numeric value.")
        }
    }

    if ("HET" %in% names(args) && !is.null(args$HET) &&
        "AE" %in% names(args) && !is.null(args$AE)) {
        stop("run_em: 'HET' and 'AE' cannot be provided simultaneously.")
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

    if ("maxnewton" %in% names(args)) {
        if (!is.numeric(args$maxnewton) || as.integer(args$maxnewton) != args$maxnewton || args$maxnewton < 1) {
            stop("Invalid 'maxnewton'. Must be a positive integer.")
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
.validate_eim <- function(X, W, V = NULL) {
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

    if (!is.null(V)) {
        if (!is.matrix(V)) {
            stop("'V' must be a matrix.")
        }
        if (nrow(V) != nrow(X)) {
            stop(
                "Mismatch in the number of ballot boxes: 'V' has ", nrow(V),
                " rows, but 'X' has ", nrow(X), " rows."
            )
        }
        if (ncol(V) < 1) {
            stop("Attribute matrix 'V' must have at least 1 column.")
        }
        if (any(is.na(V))) {
            stop("Matrix 'V' cannot contain missing values (NA).")
        }
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
#' @return A list with the `"X"` and `"W"` matrix (and `"V"` if provided). Stops execution if validation fails.
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

    if ("V" %in% names(data) && !is.null(data$V)) {
        result$V <- as.matrix(data$V)
    }

    result
}

.normalize_prob_rows_zero_columns <- function(prob_matrix) {
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

.prepare_zero_vote_columns <- function(object,
                                       initial_prob,
                                       beta_init,
                                       alpha_init,
                                       all_params,
                                       is_parametric) {
    zero_vote_cols <- which(colSums(object$X) == 0)
    if (length(zero_vote_cols) == 0) {
        return(list(
            object = object,
            initial_prob = initial_prob,
            beta_init = beta_init,
            alpha_init = alpha_init,
            all_params = all_params,
            state = NULL
        ))
    }

    keep_cols <- which(colSums(object$X) > 0)
    if (length(keep_cols) == 0) {
        stop("run_em: 'X' must contain at least one candidate with a positive vote total.")
    }

    state <- list(
        original_X = object$X,
        keep_cols = keep_cols,
        can_restore_coefficients = is_parametric && tail(keep_cols, 1) == ncol(object$X)
    )
    object$X <- object$X[, keep_cols, drop = FALSE]

    if (is.matrix(initial_prob)) {
        initial_prob <- .normalize_prob_rows_zero_columns(initial_prob[, keep_cols, drop = FALSE])
        all_params$initial_prob <- initial_prob
    }

    if (is_parametric && state$can_restore_coefficients) {
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

    list(
        object = object,
        initial_prob = initial_prob,
        beta_init = beta_init,
        alpha_init = alpha_init,
        all_params = all_params,
        state = state
    )
}

.restore_zero_vote_columns <- function(fit_object, state) {
    if (is.null(state)) {
        return(fit_object)
    }

    original_X <- state$original_X
    keep_cols <- state$keep_cols
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

    if (isTRUE(state$can_restore_coefficients) && is.matrix(fit_object$beta)) {
        coef_keep_cols <- keep_cols[-length(keep_cols)]
        restored_beta <- matrix(0, nrow = nrow(fit_object$beta), ncol = total_candidates - 1)
        if (length(coef_keep_cols) > 0 && ncol(fit_object$beta) > 0) {
            restored_beta[, coef_keep_cols] <- fit_object$beta
        }
        rownames(restored_beta) <- rownames(fit_object$beta)
        colnames(restored_beta) <- candidate_names[-total_candidates]
        fit_object$beta <- restored_beta
    }

    if (isTRUE(state$can_restore_coefficients) && is.matrix(fit_object$alpha)) {
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
#' Build the standard dimnames used by `run_em()` 3d outputs.
#'
#' @param object An `eim` object with candidate and ballot-box names.
#' @param W_matrix A group matrix providing group names.
#' @return A dimnames list for `(g x c x b)` arrays.
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
#' @return Either `object$W_agg` or `object$W`.
#' @noRd
.run_em_working_group_matrix <- function(object) {
    if (is.null(object$W_agg)) object$W else object$W_agg
}

#' Internal function!
#'
#' Apply adjacent group aggregation to an `eim` object.
#'
#' @param object An `eim` object.
#' @param group_agg Optional vector of cumulative group endpoints.
#' @return The updated object with `W_agg` and `group_agg` when requested.
#' @noRd
.run_em_apply_group_agg <- function(object, group_agg) {
    if (is.null(group_agg)) {
        return(object)
    }

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

    object
}

#' Internal function!
#'
#' Prepare the input object used by `run_em()` before dispatching the EM routine.
#'
#' @param object Optional `eim` object provided by the user.
#' @param X Candidate-vote matrix.
#' @param W Group-vote matrix.
#' @param V Optional attribute matrix for the parametric model.
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
                                   V,
                                   json_path,
                                   scale_factor,
                                   method,
                                   allow_mismatch,
                                   group_agg) {
    if (is.null(object)) {
        object <- eim(X = X, W = W, V = V, json_path = json_path)
    } else if (!inherits(object, "eim")) {
        stop("run_em: The object must be initialized with the eim() function.")
    }

    if (!is.null(V)) {
        object$V <- V
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
    } else if (method == "exact" && length(mismatch_rows) > 0) {
        .dhondt_correction(object$W, object$X)
        message("Applying a D'Hondt correction for correcting mismatches in W")
    }

    object <- .run_em_apply_group_agg(object, group_agg)
    object$method <- method
    object
}

#' Internal function!
#'
#' Populate method-specific defaults for non-parametric EM runs.
#'
#' @param object An `eim` object.
#' @param method Character string with the selected EM method.
#' @param all_params List of evaluated `run_em()` arguments.
#' @return The updated `eim` object with method-specific defaults.
#' @noRd
.run_em_apply_nonparametric_defaults <- function(object, method, all_params) {
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
#' Compute the mixture HET metric for a fitted object.
#'
#' @param fit_object An `eim` object returned by `run_em()`.
#' @return A numeric scalar, or `NA_real_` when not computable.
#' @noRd
.run_em_compute_HET_metric <- function(fit_object) {
    if (is.null(fit_object$expected_outcome)) {
        return(NA_real_)
    }

    W_fit <- .run_em_working_group_matrix(fit_object)
    if (is.null(W_fit)) {
        return(NA_real_)
    }

    denominator <- sum(W_fit)
    if (!is.finite(denominator) || denominator <= 0) {
        return(NA_real_)
    }

    z_hat <- fit_object$expected_outcome
    if (!is.array(z_hat) || length(dim(z_hat)) != 3) {
        return(NA_real_)
    }

    z_tilde <- NULL
    if (!is.null(fit_object$component_prob) && !is.null(fit_object$responsibilities)) {
        component_prob <- fit_object$component_prob
        responsibilities <- as.matrix(fit_object$responsibilities)
        if (length(dim(component_prob)) == 3 &&
            ncol(responsibilities) == dim(component_prob)[3]) {
            num_groups <- dim(component_prob)[1]
            num_candidates <- dim(component_prob)[2]
            num_ballots <- nrow(W_fit)
            z_tilde <- array(0, dim = c(num_groups, num_candidates, num_ballots))
            for (ballot in seq_len(num_ballots)) {
                for (group in seq_len(num_groups)) {
                    for (candidate in seq_len(num_candidates)) {
                        z_tilde[group, candidate, ballot] <-
                            W_fit[ballot, group] *
                                sum(component_prob[group, candidate, ] * responsibilities[ballot, ])
                    }
                }
            }
        }
    }

    if (is.null(z_tilde) && is.array(fit_object$prob) && length(dim(fit_object$prob)) == 3) {
        prob_arr <- fit_object$prob
        num_groups <- dim(prob_arr)[1]
        num_ballots <- dim(prob_arr)[3]
        z_tilde <- array(0, dim = dim(prob_arr))
        for (ballot in seq_len(num_ballots)) {
            for (group in seq_len(num_groups)) {
                z_tilde[group, , ballot] <- W_fit[ballot, group] * prob_arr[group, , ballot]
            }
        }
    }

    if (is.null(z_tilde) && is.matrix(fit_object$prob)) {
        prob_matrix <- fit_object$prob
        num_groups <- nrow(prob_matrix)
        num_candidates <- ncol(prob_matrix)
        num_ballots <- nrow(W_fit)
        z_tilde <- array(0, dim = c(num_groups, num_candidates, num_ballots))
        for (ballot in seq_len(num_ballots)) {
            for (group in seq_len(num_groups)) {
                z_tilde[group, , ballot] <- W_fit[ballot, group] * prob_matrix[group, ]
            }
        }
    }

    if (is.null(z_tilde) || !identical(dim(z_hat), dim(z_tilde))) {
        return(NA_real_)
    }

    het <- 50 * sum(abs(z_hat - z_tilde)) / denominator
    if (!is.finite(het)) {
        return(NA_real_)
    }

    het
}

#' Internal function!
#'
#' Compute the mixture AE metric for a fitted object.
#'
#' @param fit_object An `eim` object returned by `run_em()`.
#' @return A numeric scalar, or `NA_real_` when not computable.
#' @noRd
.run_em_compute_AE_metric <- function(fit_object) {
    W_fit <- .run_em_working_group_matrix(fit_object)
    X_fit <- fit_object$X
    if (is.null(W_fit) || is.null(X_fit)) {
        return(NA_real_)
    }

    prob_forward <- fit_object$prob
    if (is.array(prob_forward) && length(dim(prob_forward)) == 3) {
        num_groups <- dim(prob_forward)[1]
        num_candidates <- dim(prob_forward)[2]
        prob_matrix <- matrix(0, nrow = num_groups, ncol = num_candidates)
        for (group in seq_len(num_groups)) {
            denominator <- sum(W_fit[, group])
            if (is.finite(denominator) && denominator > 0) {
                prob_matrix[group, ] <- colSums(t(prob_forward[group, , ]) * W_fit[, group]) / denominator
            } else {
                prob_matrix[group, ] <- rep(1 / num_candidates, num_candidates)
            }
        }
        prob_forward <- prob_matrix
    }

    if (!is.matrix(prob_forward)) {
        return(NA_real_)
    }

    tau <- sum(abs(as.matrix(X_fit) - as.matrix(W_fit) %*% prob_forward))
    if (!is.finite(tau)) {
        return(NA_real_)
    }

    tau
}

#' Internal function!
#'
#' Restore zero-vote candidates and attach mixture diagnostics.
#'
#' @param fit_object An `eim` object returned by an EM branch.
#' @param zero_vote_state State returned by `.prepare_zero_vote_columns()`.
#' @return The finalized `eim` object.
#' @noRd
.run_em_finalize_fit_object <- function(fit_object, zero_vote_state) {
    fit_object <- .restore_zero_vote_columns(fit_object, zero_vote_state)
    fit_object$HET <- .run_em_compute_HET_metric(fit_object)
    fit_object$AE <- .run_em_compute_AE_metric(fit_object)
    fit_object
}

#' Internal function!
#'
#' Replace invalid empty-ballot outputs by zero conditional probabilities.
#'
#' @param fit_object An `eim` object with `cond_prob`, `expected_outcome`, and `prob`.
#' @param W_matrix Group matrix used by the fit.
#' @return The updated object.
#' @noRd
.run_em_fix_empty_ballot_outputs <- function(fit_object, W_matrix) {
    empty_ballots <- rowSums(W_matrix) == 0
    needs_empty_fix <- any(empty_ballots) ||
        any(!is.finite(fit_object$cond_prob)) ||
        any(!is.finite(fit_object$expected_outcome)) ||
        any(!is.finite(fit_object$prob))

    if (!needs_empty_fix) {
        return(fit_object)
    }

    for (ballot in which(empty_ballots)) {
        fit_object$cond_prob[, , ballot] <- 0
    }
    fit_object$cond_prob[!is.finite(fit_object$cond_prob)] <- 0

    expected_bgc <- sweep(aperm(fit_object$cond_prob, c(3, 1, 2)), c(1, 2), W_matrix, "*")
    fit_object$expected_outcome <- aperm(expected_bgc, c(2, 3, 1))
    dimnames(fit_object$expected_outcome) <- .run_em_dimnames(fit_object, W_matrix)

    numerator <- apply(fit_object$expected_outcome, c(1, 2), sum)
    denominator <- colSums(W_matrix)
    valid_groups <- is.finite(denominator) & denominator > 0
    if (any(valid_groups)) {
        fit_object$prob[valid_groups, ] <- sweep(
            numerator[valid_groups, , drop = FALSE],
            1,
            denominator[valid_groups],
            "/"
        )
        fit_object$prob[valid_groups, ] <- .normalize_prob_rows(fit_object$prob[valid_groups, , drop = FALSE])
    }

    fit_object
}

#' Internal function!
#'
#' Copy parametric EM results back into an `eim` object.
#'
#' @param object An `eim` object.
#' @param resulting_values Output list returned by `EMAlgorithmParametric`.
#' @param W_matrix Group matrix used in the fit.
#' @param V_matrix Attribute matrix used in the fit.
#' @param control List with active `run_em()` controls.
#' @return The updated `eim` object.
#' @noRd
.run_em_assign_parametric_results <- function(object, resulting_values, W_matrix, V_matrix, control) {
    dimnames_3d <- .run_em_dimnames(object, W_matrix)

    object$prob <- resulting_values$prob
    dimnames(object$prob) <- dimnames_3d
    object$cond_prob <- resulting_values$cond_prob
    dimnames(object$cond_prob) <- dimnames_3d
    object$expected_outcome <- resulting_values$expected_outcome
    dimnames(object$expected_outcome) <- dimnames_3d

    object$beta <- resulting_values$beta
    if (!is.null(colnames(W_matrix))) {
        rownames(object$beta) <- colnames(W_matrix)
    }
    if (!is.null(colnames(object$X))) {
        colnames(object$beta) <- colnames(object$X)[-ncol(object$X)]
    }

    object$alpha <- resulting_values$alpha
    if (!is.null(colnames(object$X))) {
        rownames(object$alpha) <- colnames(object$X)[-ncol(object$X)]
    }
    if (!is.null(colnames(V_matrix))) {
        colnames(object$alpha) <- colnames(V_matrix)
    }

    object$iterations <- as.numeric(resulting_values$iter)
    object$logLik <- as.numeric(resulting_values$logLik)
    object$time <- resulting_values$time
    object$maxiter <- control$maxiter
    object$maxtime <- control$maxtime
    object$ll_threshold <- control$ll_threshold
    object$maxnewton <- control$maxnewton
    object$adjust_prob_cond_method <- control$adjust_prob_cond_method
    object$adjust_prob_cond_every <- control$adjust_prob_cond_every
    object$mixture <- as.integer(control$mixture)
    object$K <- as.integer(control$mixture)
    object$row_mixture <- as.integer(control$row_mixture)

    object
}

#' Internal function!
#'
#' Build the recursive call used for the reverse symmetric run.
#'
#' @param base_call Original `run_em()` call.
#' @param object Forward fitted object.
#' @param include_V Whether to include `V` in the recursive call.
#' @param reset_parametric_init Whether to clear parametric initial values.
#' @param transform_initial_prob Whether to transpose a custom non-parametric initial probability.
#' @param all_params Evaluated arguments from the original call.
#' @return A modified call object for the reverse run.
#' @noRd
.run_em_symmetric_inverse_call <- function(base_call,
                                           object,
                                           include_V,
                                           reset_parametric_init,
                                           transform_initial_prob,
                                           all_params) {
    inverse_call <- base_call
    inverse_call$symmetric <- FALSE
    inverse_call$X <- object$W
    inverse_call$W <- object$X
    inverse_call$json_path <- NULL
    inverse_call$object <- NULL
    inverse_call$scale_factor <- 1

    if (include_V) {
        inverse_call$V <- object$V
    }

    if (reset_parametric_init) {
        inverse_call$beta_init <- NULL
        inverse_call$alpha_init <- NULL
    }

    if (transform_initial_prob && !is.null(all_params)) {
        initial_prob <- all_params$initial_prob
        if (is.matrix(initial_prob)) {
            col_tot_X <- colSums(object$X)
            numerator <- sweep(initial_prob, 2, col_tot_X, "*")
            denominator <- rowSums(numerator)
            inverse_call$initial_prob <- sweep(numerator, 1, denominator, "/")
            inverse_call$initial_prob <- t(inverse_call$initial_prob)
        }
    }

    inverse_call
}

#' Internal function!
#'
#' Evaluate a fixed non-parametric probability matrix log-likelihood.
#'
#' @param X_matrix Candidate-vote matrix.
#' @param W_matrix Group-vote matrix.
#' @param prob_matrix Probability matrix to evaluate.
#' @param method_name EM method name.
#' @return Numeric scalar log-likelihood.
#' @noRd
.run_em_loglik_from_prob <- function(X_matrix, W_matrix, prob_matrix, method_name) {
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

#' Internal function!
#'
#' Combine forward and reverse symmetric runs.
#'
#' @param object Forward fitted object.
#' @param control List with active `run_em()` controls.
#' @param W_sym Group matrix used to rebuild symmetric probabilities.
#' @param include_V Whether the reverse call should include `V`.
#' @param reset_parametric_init Whether to clear parametric initial values.
#' @param transform_initial_prob Whether to transpose custom non-parametric initial probabilities.
#' @param prob_from_cond_prob Whether to store the full `cond_prob` array as `prob`.
#' @return The symmetric `eim` object.
#' @noRd
.run_em_symmetric_helper <- function(object,
                                     control,
                                     W_sym,
                                     include_V = FALSE,
                                     reset_parametric_init = FALSE,
                                     transform_initial_prob = FALSE,
                                     prob_from_cond_prob = FALSE) {
    inverse_call <- .run_em_symmetric_inverse_call(
        control$base_call,
        object,
        include_V,
        reset_parametric_init,
        transform_initial_prob,
        control$all_params
    )
    inverse <- eval(inverse_call, control$caller_env)

    object$cond_prob_inv <- inverse$cond_prob
    object$prob_inv <- inverse$prob
    object$expected_outcome_inv <- inverse$expected_outcome
    object$time <- object$time + inverse$time
    object$iterations <- object$iterations + inverse$iterations

    weight_method <- if (!is.null(control$symmetric_weight_method)) {
        as.character(control$symmetric_weight_method)
    } else {
        "average"
    }
    em_method <- if (!is.null(control$method)) as.character(control$method) else "mult"

    weight_original <- 0.5
    weight_reverse <- 0.5
    expected_orig <- object$expected_outcome
    reverse_expected <- aperm(inverse$expected_outcome, c(2, 1, 3))

    single_profile <- (is.null(object$K) || as.integer(object$K) == 1L) &&
        (is.null(object$row_mixture) || as.integer(object$row_mixture) == 1L)
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
        identical(dim(expected_orig), dim(reverse_expected))

    if (can_use_delta) {
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

        LL <- suppressWarnings(as.numeric(object$logLik))
        if (!is.finite(LL)) {
            LL <- .run_em_loglik_from_prob(object$X, W_sym, object$prob, em_method)
        }
        LL_rev <- suppressWarnings(as.numeric(inverse$logLik))
        if (!is.finite(LL_rev)) {
            LL_rev <- .run_em_loglik_from_prob(object$W, object$X, inverse$prob, em_method)
        }

        LL_ind <- .run_em_loglik_from_prob(object$X, W_sym, p_ind, em_method)
        LL_rev_ind <- .run_em_loglik_from_prob(object$W, object$X, p_rev_ind, em_method)

        dLL <- LL - LL_ind
        dLL_rev <- LL_rev - LL_rev_ind
        dLL_pos <- max(0, dLL)
        dLL_rev_pos <- max(0, dLL_rev)
        denominator <- dLL_pos + dLL_rev_pos
        if (is.finite(denominator) && denominator > 0) {
            weight_original <- dLL_pos / denominator
            weight_reverse <- dLL_rev_pos / denominator
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

    object$expected_outcome <- weight_original * expected_orig + weight_reverse * reverse_expected
    dimnames(object$expected_outcome) <- .run_em_dimnames(object, W_sym)

    expected_bgc <- aperm(object$expected_outcome, c(3, 1, 2))
    cond_prob_bgc <- sweep(expected_bgc, c(1, 2), W_sym, "/")
    cond_prob_bgc <- .normalize_cube_last_dim(cond_prob_bgc)
    object$cond_prob <- aperm(cond_prob_bgc, c(2, 3, 1))
    dimnames(object$cond_prob) <- .run_em_dimnames(object, W_sym)

    if (prob_from_cond_prob) {
        object$prob <- object$cond_prob
        dimnames(object$prob) <- .run_em_dimnames(object, W_sym)
    } else {
        numerator <- apply(object$expected_outcome, c(1, 2), sum)
        denominator <- colSums(W_sym)
        object$prob <- sweep(numerator, 1, denominator, "/")
        object$prob <- .normalize_prob_rows(object$prob)
        dimnames(object$prob) <- list(colnames(W_sym), colnames(object$X))
    }

    object
}

#' Internal function!
#'
#' Execute the parametric branch of `run_em()`.
#'
#' @param object An `eim` object.
#' @param control List with active `run_em()` controls.
#' @return The updated `eim` object.
#' @noRd
.run_em_parametric <- function(object, control) {
    if (!identical(control$method, "mult")) {
        stop("run_em: Parametric mode only supports method = \"mult\".")
    }
    if (isTRUE(control$symmetric) && identical(control$symmetric_weight_method, "EM_weight")) {
        stop("run_em: `symmetric_weight_method = \"EM_weight\"` is currently supported only for non-parametric models (`V = NULL`).")
    }
    if (control$mixture > 1L) {
        stop("run_em: `mixture > 1` is currently supported only for non-parametric models (`V = NULL`).")
    }
    if (!is.null(control$HET)) {
        stop("run_em: `HET` is currently supported only for non-parametric models (`V = NULL`).")
    }
    if (!is.null(control$AE)) {
        stop("run_em: `AE` is currently supported only for non-parametric models (`V = NULL`).")
    }
    if (control$row_mixture > 1L) {
        stop("run_em: `row_mixture > 1` is currently supported only for non-parametric models (`V = NULL`).")
    }

    W_matrix <- .run_em_working_group_matrix(object)
    V_matrix <- object$V
    num_candidates <- ncol(object$X)
    num_groups <- ncol(W_matrix)
    num_attributes <- ncol(V_matrix)

    beta <- if (is.null(control$beta_init)) {
        matrix(0, nrow = num_groups, ncol = num_candidates - 1)
    } else {
        control$beta_init
    }
    alpha <- if (is.null(control$alpha_init)) {
        matrix(0, nrow = num_candidates - 1, ncol = num_attributes)
    } else {
        control$alpha_init
    }

    if (!is.matrix(beta) || nrow(beta) != num_groups || ncol(beta) != num_candidates - 1) {
        stop("run_em: 'beta' must be a matrix with dimensions (g x (c-1)).")
    }
    if (!is.matrix(alpha) || nrow(alpha) != num_candidates - 1 || ncol(alpha) != num_attributes) {
        stop("run_em: 'alpha' must be a matrix with dimensions ((c-1) x a).")
    }

    resulting_values <- EMAlgorithmParametric(
        as.matrix(object$X),
        as.matrix(W_matrix),
        as.matrix(V_matrix),
        as.matrix(beta),
        as.matrix(alpha),
        control$maxiter,
        control$maxtime,
        control$ll_threshold,
        control$maxnewton,
        control$verbose,
        control$adjust_prob_cond_method,
        control$adjust_prob_cond_every
    )

    object <- .run_em_assign_parametric_results(object, resulting_values, W_matrix, V_matrix, control)

    if (isTRUE(control$symmetric)) {
        object <- .run_em_symmetric_helper(
            object = object,
            control = control,
            W_sym = .run_em_working_group_matrix(object),
            include_V = TRUE,
            reset_parametric_init = TRUE,
            transform_initial_prob = FALSE,
            prob_from_cond_prob = TRUE
        )
    }

    .run_em_finalize_fit_object(object, control$zero_vote_state)
}

#' Internal function!
#'
#' Attach common non-parametric run metadata.
#'
#' @param object An `eim` object.
#' @param resulting_values Values returned by a non-parametric C routine.
#' @param W_matrix Group matrix used by the fit.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_assign_nonparametric_common <- function(object, resulting_values, W_matrix, control) {
    object$cond_prob <- resulting_values$q
    dimnames(object$cond_prob) <- .run_em_dimnames(object, W_matrix)
    object$expected_outcome <- resulting_values$expected_outcome
    dimnames(object$expected_outcome) <- .run_em_dimnames(object, W_matrix)
    object$prob <- as.matrix(resulting_values$result)
    dimnames(object$prob) <- list(colnames(W_matrix), colnames(object$X))
    object$iterations <- as.numeric(resulting_values$total_iterations)
    object <- .run_em_fix_empty_ballot_outputs(object, W_matrix)

    if (control$compute_ll) {
        log_likelihood <- resulting_values$log_likelihood
        object$logLik <- as.numeric(log_likelihood[length(log_likelihood)])
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
    object$mixture <- as.integer(control$mixture)
    object$K <- as.integer(control$K)
    object$row_mixture <- as.integer(control$row_mixture)
    object$adjust_prob_cond_method <- control$adjust_prob_cond_method
    object$adjust_prob_cond_every <- control$adjust_prob_cond_every

    object
}

#' Internal function!
#'
#' Copy row-mixture EM results back into an `eim` object.
#'
#' @param object An `eim` object.
#' @param resulting_values Output list returned by `EMAlgorithmRowMixture`.
#' @param W_matrix Group matrix used by the fit.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_assign_row_mixture_results <- function(object, resulting_values, W_matrix, control) {
    object <- .run_em_assign_nonparametric_common(object, resulting_values, W_matrix, control)
    object$mixture <- as.integer(control$H_group)
    object$K <- as.integer(control$H_group)
    object$row_mixture <- as.integer(control$H_group)

    object$phi <- as.matrix(resulting_values$phi)
    dimnames(object$phi) <- list(colnames(W_matrix), paste0("H", seq_len(control$H_group)))

    object$responsibilities <- as.matrix(resulting_values$responsibilities)
    assignment_count <- ncol(object$responsibilities)
    dimnames(object$responsibilities) <- list(rownames(object$X), paste0("A", seq_len(assignment_count)))

    object$component_prob <- resulting_values$component_prob
    dimnames(object$component_prob) <- list(
        colnames(W_matrix),
        colnames(object$X),
        paste0("H", seq_len(control$H_group))
    )

    object
}

#' Internal function!
#'
#' Copy mixture EM results back into an `eim` object.
#'
#' @param object An `eim` object.
#' @param resulting_values Output list returned by `EMAlgorithmMixture`.
#' @param W_matrix Group matrix used by the fit.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_assign_mixture_results <- function(object, resulting_values, W_matrix, control) {
    object <- .run_em_assign_nonparametric_common(object, resulting_values, W_matrix, control)
    object$mixture <- as.integer(control$mixture)
    object$K <- as.integer(control$K)
    object$row_mixture <- as.integer(control$H_group)

    expose_mixture_details <- control$K > 1L || isTRUE(control$mixture_explicit)
    if (expose_mixture_details) {
        object$phi <- as.numeric(resulting_values$phi)
        names(object$phi) <- paste0("H", seq_len(control$K))
        object$responsibilities <- as.matrix(resulting_values$responsibilities)
        dimnames(object$responsibilities) <- list(rownames(object$X), paste0("H", seq_len(control$K)))
        object$component_prob <- resulting_values$component_prob
        dimnames(object$component_prob) <- list(
            colnames(W_matrix),
            colnames(object$X),
            paste0("H", seq_len(control$K))
        )
    } else {
        object$phi <- NULL
        object$responsibilities <- NULL
        object$component_prob <- NULL
    }

    object
}

#' Internal function!
#'
#' Execute the row-mixture branch of `run_em()`.
#'
#' @param object An `eim` object.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_row_mixture <- function(object, control) {
    W_matrix <- .run_em_working_group_matrix(object)
    resulting_values <- EMAlgorithmRowMixture(
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
        control$H_group
    )

    object <- .run_em_assign_row_mixture_results(object, resulting_values, W_matrix, control)

    if (isTRUE(control$symmetric)) {
        object <- .run_em_symmetric_helper(
            object = object,
            control = control,
            W_sym = .run_em_working_group_matrix(object),
            include_V = FALSE,
            reset_parametric_init = FALSE,
            transform_initial_prob = TRUE,
            prob_from_cond_prob = FALSE
        )
    }

    .run_em_finalize_fit_object(object, control$zero_vote_state)
}

#' Internal function!
#'
#' Execute the finite-mixture branch of `run_em()`.
#'
#' @param object An `eim` object.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_mixture <- function(object, control) {
    W_matrix <- .run_em_working_group_matrix(object)
    resulting_values <- EMAlgorithmMixture(
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
        control$K,
        control$use_joint_symmetric_em,
        if (control$use_joint_symmetric_em) "EM_weight" else "average"
    )

    object <- .run_em_assign_mixture_results(object, resulting_values, W_matrix, control)

    if (control$use_joint_symmetric_em) {
        object$cond_prob_inv <- resulting_values$q_inv
        dimnames(object$cond_prob_inv) <- list(colnames(object$X), colnames(W_matrix), rownames(object$X))
        object$prob_inv <- as.matrix(resulting_values$prob_inv)
        dimnames(object$prob_inv) <- list(colnames(object$X), colnames(W_matrix))
        object$expected_outcome_inv <- resulting_values$expected_outcome_inv
        dimnames(object$expected_outcome_inv) <- list(colnames(object$X), colnames(W_matrix), rownames(object$X))
        object$symmetric_weight_method <- "EM_weight"
        object$symmetric_weights <- c(original = 0.5, reverse = 0.5)
    } else if (isTRUE(control$symmetric)) {
        object <- .run_em_symmetric_helper(
            object = object,
            control = control,
            W_sym = .run_em_working_group_matrix(object),
            include_V = FALSE,
            reset_parametric_init = FALSE,
            transform_initial_prob = TRUE,
            prob_from_cond_prob = FALSE
        )
    }

    .run_em_finalize_fit_object(object, control$zero_vote_state)
}

#' Internal function!
#'
#' Execute the standard non-parametric branch of `run_em()`.
#'
#' @param object An `eim` object.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_standard_nonparametric <- function(object, control) {
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
        if (is.matrix(control$initial_prob)) control$initial_prob else matrix(-1, nrow = 1, ncol = 1)
    )

    object <- .run_em_assign_nonparametric_common(object, resulting_values, W_matrix, control)

    if (isTRUE(control$symmetric)) {
        object <- .run_em_symmetric_helper(
            object = object,
            control = control,
            W_sym = object$W,
            include_V = FALSE,
            reset_parametric_init = FALSE,
            transform_initial_prob = TRUE,
            prob_from_cond_prob = FALSE
        )
    }

    .run_em_finalize_fit_object(object, control$zero_vote_state)
}

#' Internal function!
#'
#' Search over K or row_mixture using a mixture diagnostic metric.
#'
#' @param control List with active `run_em()` controls.
#' @param metric_name Either `"HET"` or `"AE"`.
#' @param threshold Numeric threshold.
#' @return The selected fitted object.
#' @noRd
.run_em_adaptive_search_metric <- function(control, metric_name, threshold) {
    metric_fun <- switch(metric_name,
        HET = .run_em_compute_HET_metric,
        AE = .run_em_compute_AE_metric,
        stop("run_em: Unsupported adaptive metric.")
    )

    K_max <- 7L
    best_fit <- NULL
    best_metric <- Inf
    last_fit <- NULL
    search_best_only <- isTRUE(all.equal(threshold, 0)) || identical(metric_name, "AE")
    search_over_H <- control$row_mixture > 1L

    if (search_over_H) {
        current <- 1L
        repeat {
            metric_call <- control$base_call
            metric_call$object <- control$full_candidate_object
            metric_call$X <- NULL
            metric_call$W <- NULL
            metric_call$V <- NULL
            metric_call$json_path <- NULL
            metric_call$mixture <- current
            metric_call$row_mixture <- current
            metric_call$H <- NULL
            metric_call$HET <- NULL
            metric_call$AE <- NULL

            fit_metric <- eval(metric_call, control$caller_env)
            last_fit <- fit_metric
            metric_value <- metric_fun(fit_metric)
            fit_metric[[metric_name]] <- metric_value
            fit_metric$K <- as.integer(current)
            fit_metric$mixture <- as.integer(current)
            fit_metric$row_mixture <- as.integer(current)

            if (is.finite(metric_value) && metric_value < best_metric) {
                best_metric <- metric_value
                best_fit <- fit_metric
            }
            if (!search_best_only && is.finite(metric_value) && metric_value < threshold) {
                return(fit_metric)
            }
            if (current >= K_max) {
                if (!is.null(best_fit)) {
                    if (!search_best_only) {
                        warning(sprintf(
                            "run_em: No row_mixture <= %d satisfied the %s threshold. Returning the fit with minimum %s.",
                            K_max,
                            metric_name,
                            metric_name
                        ))
                    }
                    return(best_fit)
                }
                warning(sprintf(
                    "run_em: No finite %s found for row_mixture <= %d. Returning the last fitted model.",
                    metric_name,
                    K_max
                ))
                return(last_fit)
            }

            current <- current + 1L
        }
    }

    current <- 1L
    repeat {
        metric_call <- control$base_call
        metric_call$object <- control$full_candidate_object
        metric_call$X <- NULL
        metric_call$W <- NULL
        metric_call$V <- NULL
        metric_call$json_path <- NULL
        metric_call$mixture <- current
        metric_call$row_mixture <- 1L
        metric_call$H <- NULL
        metric_call$HET <- NULL
        metric_call$AE <- NULL

        fit_metric <- eval(metric_call, control$caller_env)
        last_fit <- fit_metric
        metric_value <- metric_fun(fit_metric)
        fit_metric[[metric_name]] <- metric_value
        fit_metric$K <- as.integer(current)
        fit_metric$mixture <- as.integer(current)
        fit_metric$row_mixture <- 1L

        if (is.finite(metric_value) && metric_value < best_metric) {
            best_metric <- metric_value
            best_fit <- fit_metric
        }
        if (!search_best_only && is.finite(metric_value) && metric_value < threshold) {
            return(fit_metric)
        }
        if (current >= K_max) {
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

        current <- current + 1L
    }
}

#' Internal function!
#'
#' Execute the non-parametric branch of `run_em()`.
#'
#' @param object An `eim` object.
#' @param control List with active `run_em()` controls.
#' @return The updated object.
#' @noRd
.run_em_nonparametric <- function(object, control) {
    if (isTRUE(control$symmetric) &&
        identical(control$symmetric_weight_method, "EM_weight") &&
        control$use_row_mixture) {
        stop("run_em: `symmetric_weight_method = \"EM_weight\"` is currently supported only when `row_mixture = 1`.")
    }

    if (!is.null(control$HET)) {
        return(.run_em_adaptive_search_metric(control, "HET", control$HET))
    }
    if (!is.null(control$AE)) {
        return(.run_em_adaptive_search_metric(control, "AE", control$AE))
    }

    object <- .run_em_apply_nonparametric_defaults(object, control$method, control$all_params)

    if (control$use_row_mixture) {
        return(.run_em_row_mixture(object, control))
    }
    if (control$use_mixture) {
        return(.run_em_mixture(object, control))
    }

    .run_em_standard_nonparametric(object, control)
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
