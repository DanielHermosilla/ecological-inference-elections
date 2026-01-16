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

#' @title Plot estimated probabilities
#' @description
#'   Plots the estimated probabilities as pie charts using `ggplot2`, one per row of the probability matrix.
#'   Each slice displays its percentage label. For the parametric case, it does a weighted average over groups to retrieve the global probabilities.
#'
#' @param x An "eim" object.
#' @param title Title for the plot.
#' @param legend_title Title for the legend.
#' @param color_scale A vector of colors or a palette for the candidates.
#' @param min_pct Minimum percentage required to display a label.
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
#'     ballot_voters = 40,
#'     parametric = TRUE,
#'     num_attributes = 2,
#'     num_districts = 2,
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
                     ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("plot.eim: package 'ggplot2' is required for plotting.")
    }
    object <- x
    prob <- object$prob

    # Case it's parametric
    if (is.array(prob) && length(dim(prob)) == 3) {
        W <- as.matrix(object$W)
        G <- dim(prob)[1]
        C <- dim(prob)[2]
        P <- matrix(0, nrow = G, ncol = C)
        for (g in seq_len(G)) {
            P[g, ] <- colSums(t(prob[g, , ]) * W[, g]) / sum(W[, g])
        }
        dimnames(P) <- dimnames(prob)[1:2]
    } else {
        P <- as.matrix(prob)
    }

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

    plot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = 0.5, y = value, fill = col)) +
        ggplot2::geom_col(width = 1, color = NA) +
        ggplot2::coord_polar(theta = "y", clip = "off") +
        ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        ggplot2::facet_wrap(~row, strip.position = "top") +
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

#' Returns the object estimated probability. In case it is a parametric model, it returns the global probabilities weighted by W.
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

    # Return the global probabilities in case of parametric model, weighted by W
    if (is.array(prob) && length(dim(prob)) == 3) {
        W <- as.matrix(object$W)
        G <- dim(prob)[1]
        C <- dim(prob)[2]
        P <- matrix(0, nrow = G, ncol = C)
        for (g in seq_len(G)) {
            P[g, ] <- colSums(t(prob[g, , ]) * W[, g]) / sum(W[, g])
        }
        dimnames(P) <- dimnames(prob)[1:2]
    } else {
        P <- as.matrix(prob)
    }
    return(P)
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
