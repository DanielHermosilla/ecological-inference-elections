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
    is_parametric <- !is.null(object$V)

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

    if (is_parametric) {
        cat("Attribute matrix (V) [b x a]:\n")
        print(object$V[1:min(5, nrow(object$V)), , drop = FALSE]) # nolint
        truncated_V <- (nrow(object$V) > 5)
        if (truncated_V) cat(".\n.\n.\n") else cat("\n")

        if (!is.null(object$beta)) {
            cat("Estimated beta parameters [g x (c-1)]:\n")
            print(object$beta, drop = FALSE)
            cat("\n")
        }
        if (!is.null(object$alpha)) {
            cat("Estimated alpha parameters [(c-1) x a]:\n")
            print(object$alpha, drop = FALSE)
            cat("\n")
        }
        if (!is.null(object$sd_beta)) {
            cat("Standard deviation for beta [g x (c-1)]:\n")
            print(object$sd_beta, drop = FALSE)
            cat("\n")
        }
        if (!is.null(object$sd_alpha)) {
            cat("Standard deviation for alpha [(c-1) x a]:\n")
            print(object$sd_alpha, drop = FALSE)
            cat("\n")
        }
        # if (!is.null(object$prob)) {
        #     prob_dim <- dim(object$prob)
        #     cat("Estimated probabilities [g x c x b]: ", paste(prob_dim, collapse = " x "), "\n", sep = "")
        #     cat("\n")
        # }
    } else {
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
#' @param object An "eim" object.
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
    # Generates the list with the core attribute.
    object_core_attr <- list(
        candidates = ncol(object$X),
        groups = ncol(object$W),
        ballots = nrow(object$X)
    )

    if (!is.null(object$V)) {
        object_core_attr$attributes <- ncol(object$V)
    }

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

    if (!is.null(object$beta)) {
        object_run_em_attr$beta <- object$beta
    }
    if (!is.null(object$alpha)) {
        object_run_em_attr$alpha <- object$alpha
    }

    # Display sd if the bootstrapping method has been called.
    if (!is.null(object$sd)) {
        object_run_em_attr$sd <- object$sd
    }
    if (!is.null(object$avg_prob)) {
        object_run_em_attr$avg_prob <- object$avg_prob
    }
    if (!is.null(object$sd_beta)) {
        object_run_em_attr$sd_beta <- object$sd_beta
    }
    if (!is.null(object$sd_alpha)) {
        object_run_em_attr$sd_alpha <- object$sd_alpha
    }

    final_list <- c(object_core_attr, object_run_em_attr)
    final_list
}

#' Returns the object estimated probability
#'
#' @param object An "eim" object.
#' @param ... Additional arguments that are ignored.
#' @return The probability matrix or array
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
