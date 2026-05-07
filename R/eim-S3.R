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
    if (!is.null(object$avg_prob)) {
        object_run_em_attr$avg_prob <- object$avg_prob
    }
    final_list <- c(object_core_attr, object_run_em_attr)
    final_list
}

#' @title Plot estimated probabilities
#' @description
#'   Plots the estimated probabilities as pie charts using `ggplot2`, one per row of the probability matrix.
#'   Each slice displays its percentage label.
#'
#' @param x An "eim" object.
#' @param title Title for the plot.
#' @param legend_title Title for the legend.
#' @param color_scale A vector of colors or a palette for the candidates.
#' @param min_pct Minimum percentage required to display a label.
#' @param pies_per_row Number of pie charts to display per row. Defaults to `ceiling(sqrt(G))`,
#'   where `G` is the number of groups.
#' @param ... Additional arguments that are ignored.
#'
#' @return Returns a `ggplot2` object representing the pie charts.
#'
#' @examples
#' \donttest{
#' sim <- simulate_election(
#'     num_ballots = 100,
#'     num_candidates = 4,
#'     num_groups = 5,
#'     ballot_voters = rep(40, 100),
#'     seed = 42
#' )
#' fit <- run_em(sim, maxiter = 5)
#'
#' plot(fit, title = "Estimated probabilities", legend_title = "Candidates", min_pct = 7)
#' }
#' @export
plot.eim <- function(x,
                     title = "Estimated probabilities",
                     legend_title = "Candidates",
                     color_scale = NULL,
                     min_pct = 3,
                     pies_per_row = NULL,
                     ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot.eim: package 'ggplot2' is required for plotting.")
    }
    object <- x
    prob <- object$prob

    W_use <- if (!is.null(object$W_agg)) object$W_agg else object$W
    G <- ncol(W_use)

    P <- as.matrix(prob)

    row_names <- rownames(P)
    col_names <- colnames(P)
    if (is.null(row_names)) row_names <- paste0("Row ", seq_len(nrow(P)))
    if (is.null(col_names)) col_names <- paste0("Col ", seq_len(ncol(P)))

    if (is.function(color_scale)) {
        colors <- color_scale(ncol(P))
    } else {
        colors <- color_scale
    }
    if (is.null(colors)) {
        colors <- grDevices::colorRampPalette(c("#4575B4", "#F7F7F7", "#D73027"))(ncol(P))
    } else if (length(colors) < ncol(P)) {
        colors <- grDevices::colorRampPalette(colors)(ncol(P))
    }

    df <- expand.grid(row = row_names, col = col_names, stringsAsFactors = FALSE)
    df$value <- as.vector(P)
    df$row <- factor(df$row, levels = row_names)
    df$col <- factor(df$col, levels = col_names)
    df$label <- ifelse(100 * df$value >= min_pct, sprintf("%.1f%%", 100 * df$value), "")

    if (is.null(pies_per_row)) {
        pies_per_row <- ceiling(sqrt(G))
    }
    pies_per_row <- max(1L, as.integer(pies_per_row))

    plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = 0.5, y = value, fill = col)) +
        ggplot2::geom_col(width = 1, color = NA) +
        ggplot2::coord_polar(theta = "y", clip = "off") +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::facet_wrap(~row, strip.position = "top", ncol = pies_per_row) +
        ggplot2::geom_text(
            ggplot2::aes(label = label),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 3
        ) +
        ggplot2::scale_fill_manual(values = colors, name = legend_title) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold"),
            legend.position = "right"
        )

    print(plot_obj)
    invisible(plot_obj)
}

#' Returns the object estimated probability.
#'
#' @param object An "eim" object.
#' @param ... Additional arguments that are ignored.
#' @return The global probability matrix
#' @noRd
#' @export
as.matrix.eim <- function(x, ...) {
    object <- x
    if (is.null(object$prob)) {
        stop(paste0(
            "Probability matrix not available. Run 'run_em()'."
        ))
    }
    prob <- object$prob

    P <- as.matrix(prob)
    return(P)
}

#' @title Extract log-likelihood
#' @description
#'   Return the log-likelihood of the last EM iteration, or compute it for a fixed probability matrix.
#'
#' @param object An `eim` object
#' @param prob Optional probability matrix `(g x c)` used to compute a fixed-parameter log-likelihood.
#' @param method Optional E-step method. Defaults to the fitted method in `object`, or `"mult"` if unavailable.
#' @param ... Additional optional parameters for fixed-probability evaluation:
#'   `mcmc_stepsize`, `mcmc_samples`, `mvncdf_method`, `mvncdf_error`, `mvncdf_samples`,
#'   `miniter`, `adjust_prob_cond_method`, and `adjust_prob_cond_every`.
#'
#' @return A numeric value with the log-likelihood.
#' @noRd
#' @export
logLik.eim <- function(object, prob = NULL, method = NULL, ...) {
    if (!inherits(object, "eim")) {
        stop("The object must be initialized with the `eim()` function.")
    }

    if (is.null(prob)) {
        if (is.null(object$logLik)) {
            stop("The `run_em()` method must be called for getting the log-likelihood.")
        }
        return(tail(object$logLik, 1))
    }

    method_alias <- function(m) {
        m <- trimws(tolower(as.character(m)))
        if (m == "multinomial") {
            return("mult")
        }
        m
    }

    method_use <- if (is.null(method)) {
        if (!is.null(object$method)) object$method else "mult"
    } else {
        method
    }
    method_use <- method_alias(method_use)
    valid_methods <- c("mult", "mcmc", "mvn_cdf", "mvn_pdf", "exact")
    if (!(method_use %in% valid_methods)) {
        stop("Invalid 'method'. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    if (is.null(object$X) || is.null(object$W)) {
        stop("The object must include 'X' and 'W' matrices.")
    }

    W_use <- if (is.null(object$W_agg)) object$W else object$W_agg
    prob_matrix <- as.matrix(prob)
    expected_dims <- c(ncol(W_use), ncol(object$X))
    if (!all(dim(prob_matrix) == expected_dims)) {
        stop(
            "Invalid 'prob' dimensions. Expected (", expected_dims[1], " x ", expected_dims[2], ")."
        )
    }

    dots <- list(...)
    step_size <- if (!is.null(dots$mcmc_stepsize)) dots$mcmc_stepsize else if (!is.null(object$mcmc_stepsize)) object$mcmc_stepsize else 3000
    samples <- if (!is.null(dots$mcmc_samples)) dots$mcmc_samples else if (!is.null(object$mcmc_samples)) object$mcmc_samples else 1000
    monte_method <- if (!is.null(dots$mvncdf_method)) dots$mvncdf_method else if (!is.null(object$mvncdf_method)) object$mvncdf_method else "genz"
    monte_error <- if (!is.null(dots$mvncdf_error)) dots$mvncdf_error else if (!is.null(object$mvncdf_error)) object$mvncdf_error else 1e-3
    monte_iter <- if (!is.null(dots$mvncdf_samples)) dots$mvncdf_samples else if (!is.null(object$mvncdf_samples)) object$mvncdf_samples else 5000
    miniter <- if (!is.null(dots$miniter)) dots$miniter else if (!is.null(object$miniter)) object$miniter else 0
    lp_method <- if (!is.null(dots$adjust_prob_cond_method)) dots$adjust_prob_cond_method else if (!is.null(object$adjust_prob_cond_method)) object$adjust_prob_cond_method else ""
    project_every <- if (!is.null(dots$adjust_prob_cond_every)) dots$adjust_prob_cond_every else if (!is.null(object$adjust_prob_cond_every)) object$adjust_prob_cond_every else FALSE

    ll <- EMLogLikFromProb(
        t(object$X),
        W_use,
        prob_matrix,
        method_use,
        as.integer(step_size),
        as.integer(samples),
        as.character(monte_method),
        as.numeric(monte_error),
        as.integer(monte_iter),
        as.integer(miniter),
        as.character(lp_method),
        as.logical(project_every)
    )
    as.numeric(ll)
}
