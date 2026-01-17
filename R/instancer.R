#' Simulate an Election
#'
#' @description
#' This function simulates an election by creating matrices representing candidate votes `(X)` and voters' demographic group `(W)` across a specified number of ballot-boxes. It either (i) receives as input or (ii) generates a probability matrix `(prob)`, indicating how likely each demographic group is to vote for each candidate.
#' It supports both non-parametric and parametric simulations; set `num_covariates` and `num_districts` greater than zero to generate `V`, `real_alpha`, and `real_beta`.
#'
#' By default, the number of voters per ballot box `(ballot_voters)` is set to a vector of 100 with
#' length `num_ballots`. You can optionally override this by providing a custom vector.
#'
#' Optional parameters are available to control the distribution of votes:
#' \itemize{
#'   \item \strong{`group_proportions`}: A vector of length `num_groups` specifying
#'         the overall proportion of each demographic group. Entries must sum to one and be non-negative.
#'   \item \strong{`prob`}: A user-supplied probability matrix of dimension
#'         (`num_groups` \eqn{\times} `num_candidates`). If provided, this matrix is used directly. Otherwise, voting probabilities for each group are drawn from a Dirichlet distribution.
#' }
#'
#' @param num_ballots
#'   Number of ballot boxes (`b`).
#'
#' @param num_candidates
#'   Number of candidates (`c`).
#'
#' @param num_groups
#'   Number of demographic groups (`g`).
#'
#' @param ballot_voters
#'   A vector of length `num_ballots` representing the number of voters per ballot
#'   box. Defaults to `rep(100, num_ballots)`.
#'
#' @param lambda
#' A numeric value between 0 and 1 that represents the fraction of voters that are randomly assigned to ballot-boxes. The remaining voters are assigned sequentially according to their demographic group.
#'   - `lambda = 0`: The assignment of voters to ballot-boxes is fully sequential in terms of their demographic group. This leads to a **high heterogeneity** of the voters' groups across ballot-boxes.
#'   - `lambda = 1`: The assignment of voters to ballot-boxes is fully random. This leads to a **low heterogeneity** of the voters' group across ballot-boxes.
#'
#'  Default value is set to `0.5`. See **Shuffling Mechanish** for more details.
#'
#' @param seed
#'   If provided, overrides the current global seed. Defaults to `NULL`.
#'
#' @param group_proportions
#'   Optional. A vector of length `num_groups` that indicates the fraction of voters that belong to each group. Default is that all groups are of the same size.
#'
#' @param prob
#'   Optional. A user-supplied probability matrix of dimension `(g x c)`.
#'   If provided, this matrix is used as the underlying voting probability distribution. If not supplied, each row is sampled from a Dirichlet distribution with each parameter set to one.
#'
#' @param num_covariates
#'   Number of covariates (`a`) used to build the parametric covariates matrix `V`.
#'
#' @param num_districts
#'   Number of districts used to assign ballot boxes when `parametric = TRUE`.
#'
#' @section Shuffling Mechanism:
#' Without loss of generality, consider an order relation of groups and ballot-boxes. The shuffling step is controlled by the `lambda` parameter and operates as follows:
#'
#' 1. **Initial Assignment**: Voters are assigned to each ballot-box sequentially according to their demographic
#' group. More specifically, the first ballot-boxes receive voters from the first group. Then, the next ballot-boxes receive voters from the second group. This continues until all voters have been assigned. Note that most ballot-boxes will contain voters from a single group (as long as the number of ballot-boxes exceeds the number of groups).
#'
#' 2. **Shuffling**: A fraction `lambda` of voters who have already been assigned is selected at random. Then, the ballot-box assignment of this sample is shuffled. Hence, different `lambda` values are interpreted as follows:
#'
#' 	- `lambda = 0` means no one is shuffled (the initial assignment remains).
#' 	- `lambda = 1` means all individuals are shuffled.
#' 	- Intermediate values like `lambda = 0.5` shuffle half the voters.
#'
#' Using a high level of `lambda` (greater than 0.7) is not recommended, as this could make identification of the voting probabilities difficult. This is because higher values of lambda induce similar ballot-boxes in terms of voters' group.
#'
#' @return An eim object. For the non-parametric case it contains:
#' \describe{
#'   \item{\code{X}}{A \code{(b x c)} matrix with candidates' votes for each ballot box.}
#'   \item{\code{W}}{A \code{(b x g)} matrix with voters' groups for each ballot-box.}
#'   \item{\code{real_prob}}{A \code{(g x c)} matrix with the probability that a voter from each group votes for each candidate. If prob is provided, it would equal such probability.}
#' 	 \item{\code{outcome}}{A \code{(b x g x c)} array with the number of votes for each candidate in each ballot box, broken down by group.}
#' }
#'
#' When `parametric = TRUE`, it returns:
#' \describe{
#'   \item{\code{X}}{A \code{(b x c)} matrix with candidates' votes for each ballot box.}
#'   \item{\code{W}}{A \code{(b x g)} matrix with voters' groups for each ballot-box.}
#'   \item{\code{V}}{A \code{(b x a)} matrix with ballot-box attributes.}
#'   \item{\code{real_prob}}{A \code{(g x c x b)} array with ballot-box probabilities.}
#'   \item{\code{real_alpha}}{A \code{((c-1) x a)} matrix of true attribute parameters.}
#'   \item{\code{real_beta}}{A \code{(g x (c-1))} matrix of true group parameters.}
#' }
#'
#' @references
#' The algorithm is fully explained in ['Thraves, C. Ubilla, P. and Hermosilla, D.: *"A Fast Ecological Inference Algorithm for the R×C Case"*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834).
#'
#' @examples
#' # Example 1: Default usage with 200 ballot boxes, each having 100 voters
#' result1 <- simulate_election(
#'     num_ballots = 200,
#'     num_candidates = 3,
#'     num_groups = 5
#' )
#'
#' # Example 2: Using a custom ballot_voters vector
#' result2 <- simulate_election(
#'     num_ballots = 340,
#'     num_candidates = 3,
#'     num_groups = 7,
#'     ballot_voters = rep(200, 340)
#' )
#'
#' # Example 3: Supplying group_proportions
#' result3 <- simulate_election(
#'     num_ballots = 93,
#'     num_candidates = 3,
#'     num_groups = 4,
#'     group_proportions = c(0.3, 0.5, 0.1, 0.1)
#' )
#'
#' # Example 4: Providing a user-defined prob matrix
#' custom_prob <- matrix(c(
#'     0.9,  0.1,
#'     0.4,  0.6,
#'     0.25, 0.75,
#'     0.32, 0.68,
#'     0.2,  0.8
#' ), nrow = 5, byrow = TRUE)
#'
#' result4 <- simulate_election(
#'     num_ballots = 200,
#'     num_candidates = 2,
#'     num_groups = 5,
#'     lambda = 0.3,
#'     prob = custom_prob
#' )
#'
#' result4$real_prob == custom_prob
#' # The attribute of the output real_prob matches the input custom_prob.
#' @export
simulate_election <- function(num_ballots,
                              num_candidates,
                              num_groups,
                              ballot_voters = rep(100, num_ballots),
                              lambda = 0.5,
                              seed = NULL,
                              group_proportions = rep(1 / num_groups, num_groups),
                              prob = NULL,
                              num_covariates = 0,
                              num_districts = 0) {
    # If user provides a seed, override the current global seed
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # If the user gives num_covariates and num_districts > 0, set parametric to TRUE
    parametric <- ifelse(num_covariates > 0 && num_districts > 0, TRUE, FALSE)
    if (num_covariates > 0 && num_districts == 0 || num_districts > 0 && num_covariates == 0) {
        stop("`num_covariates` and `num_districts` must both be greater than 0 to enable parametric mode.")
    }

    if (parametric) {
        if (length(ballot_voters) != 1 && length(ballot_voters) != num_ballots) {
            stop("`ballot_voters` must be length 1 or length `num_ballots` in parametric mode.")
        }
        if (length(ballot_voters) == num_ballots && length(unique(ballot_voters)) != 1) {
            stop("`ballot_voters` must be constant across ballot boxes in parametric mode.")
        }
    } else {
        # Validate length of ballot_voters
        if (length(ballot_voters) != num_ballots) {
            stop("`ballot_voters` must be a vector of length `num_ballots`.")
        }
    }

    # Validate group proportions
    if (length(group_proportions) != num_groups) {
        stop("group_proportions must have length equal to num_groups.")
    }
    tol <- .Machine$double.eps^0.5
    if (abs(sum(group_proportions) - 1) > tol) {
        stop("group_proportions must sum to 1.")
    }
    if (lambda < 0 || lambda > 1) {
        stop("Lambda must be a decimal between 0 and 1.")
    }

    if (parametric) {
        if (!is.null(prob)) {
            warning("simulate_election: 'prob' is ignored in parametric mode.")
        }
        num_voters <- ballot_voters[1]

        alpha <- matrix(
            rnorm((num_candidates - 1) * num_covariates),
            nrow = num_candidates - 1,
            ncol = num_covariates
        )
        beta <- matrix(
            rnorm(num_groups * (num_candidates - 1)),
            nrow = num_groups,
            ncol = num_candidates - 1
        )

        random_one_in_each_row <- function(n, m) {
            mat <- matrix(0L, nrow = n, ncol = m)
            idx <- sample.int(m, n, replace = TRUE)
            mat[cbind(seq_len(n), idx)] <- 1L
            d <- min(n, m)
            mat[seq_len(d), seq_len(d)] <- diag(1L, d)
            mat
        }

        e_bd <- random_one_in_each_row(num_ballots, num_districts)

        v_da <- matrix(
            runif(num_districts * num_covariates, min = 0, max = 1),
            nrow = num_districts,
            ncol = num_covariates
        )
        v_ba <- e_bd %*% v_da

        p_dgc <- array(0, dim = c(num_districts, num_groups, num_candidates))
        for (d in seq_len(num_districts)) {
            for (g in seq_len(num_groups)) {
                u <- alpha %*% v_da[d, ] + beta[g, ]
                u <- c(u, 0)
                p_dgc[d, g, ] <- exp(u) / sum(exp(u))
            }
        }

        W <- matrix(0L, nrow = num_ballots, ncol = num_groups)
        X <- matrix(0L, nrow = num_ballots, ncol = num_candidates)

        for (d in seq_len(num_districts)) {
            boxes <- which(e_bd[, d] == 1L)
            B_d <- length(boxes)
            I_d <- num_voters * B_d

            base_counts <- floor(I_d * group_proportions)
            remainder <- I_d - sum(base_counts)
            if (remainder > 0) {
                fractional_parts <- I_d * group_proportions - base_counts
                indices <- order(fractional_parts, decreasing = TRUE)
                base_counts[indices[seq_len(remainder)]] <- base_counts[indices[seq_len(remainder)]] + 1
            }

            omega0 <- rep(seq_len(num_groups), times = base_counts)
            if (length(omega0) > I_d) {
                omega0 <- omega0[seq_len(I_d)]
            }
            if (length(omega0) < I_d) {
                omega0 <- c(
                    omega0,
                    sample(seq_len(num_groups), I_d - length(omega0), replace = TRUE)
                )
            }
            omega <- omega0

            n_mix <- round(lambda * I_d)
            if (n_mix > 0) {
                idx <- sample.int(I_d, n_mix)
                vals_to_mix <- omega0[idx]
                omega[idx] <- sample(vals_to_mix, length(vals_to_mix))
            }

            for (j in seq_along(boxes)) {
                b <- boxes[j]
                start <- (j - 1) * num_voters + 1
                end <- j * num_voters
                G_b <- omega[start:end]
                votes_b <- vapply(G_b, function(g) {
                    sample.int(num_candidates, 1, prob = p_dgc[d, g, ])
                }, integer(1))
                W[b, ] <- tabulate(G_b, nbins = num_groups)
                X[b, ] <- tabulate(votes_b, nbins = num_candidates)
            }
        }

        real_prob <- array(0, dim = c(num_groups, num_candidates, num_ballots))
        for (b in seq_len(num_ballots)) {
            d <- which(e_bd[b, ] == 1L)
            real_prob[, , b] <- p_dgc[d, , ]
        }

        toReturn <- list(
            W = W,
            X = X,
            V = v_ba,
            real_prob = real_prob,
            real_alpha = alpha,
            real_beta = beta
        )

        class(toReturn) <- "eim"
        return(toReturn)
    }

    # Build W (b x g), the distribution of groups in each ballot box
    W <- matrix(0, nrow = num_ballots, ncol = num_groups)

    # Get the total amount of voters
    n <- sum(ballot_voters)

    # Group assignment according group proportions
    counts <- floor(n * group_proportions)
    # Case when it's not an exact distribution
    remainder <- n - sum(counts)
    if (remainder > 0) {
        # Arbitrarily distribute the "fractional" person to the highest fractional parts
        fractional_parts <- n * group_proportions - counts
        indices <- order(fractional_parts, decreasing = TRUE)
        counts[indices[seq_len(remainder)]] <- counts[indices[seq_len(remainder)]] + 1
    }
    # Assign each of the n people to a group (uniformly) in a vector 'omega0'
    omega0 <- rep(seq_len(num_groups), times = counts)
    # omega0 <- ceiling(seq_len(n) * num_groups / n)
    # Copy for shuffling
    omega <- omega0

    # Shuffle lambda * n individuals among groups
    size_shuffle <- floor(lambda * n)
    if (size_shuffle > 0) {
        v <- sample.int(n, size_shuffle, replace = FALSE)
        v_sorted <- sort(v)
        for (k in seq_len(size_shuffle)) {
            # Shuffle individuals
            omega[v[k]] <- omega0[v_sorted[k]]
        }
    }

    # Get the cumulative sum for later creating the partition between ballot boxes
    # This will allow defining a range of voters between values of a[]. For
    # example: a[2] - a[1] = voters of ballot box 1.
    a <- c(0, cumsum(ballot_voters))

    # Fill in W by counting how many from each group are in each box
    for (b_idx in seq_len(num_ballots)) {
        start_idx <- a[b_idx] + 1
        end_idx <- a[b_idx + 1]
        sub_omega <- omega[start_idx:end_idx]
        counts_bg <- tabulate(sub_omega, nbins = num_groups)
        W[b_idx, ] <- counts_bg
    }

    # Build the probality matrix (g x c)
    # 1) If 'prob' is provided, use it directly
    # 2) Otherwise, sample each row from a Dirichlet with alpha=1
    if (!is.null(prob)) {
        # Validate dimension
        if (!all(dim(prob) == c(num_groups, num_candidates))) {
            stop("`prob` must be a matrix with dimensions (num_groups x num_candidates).")
        }
        p <- prob
    } else {
        # Sample from Dirichlet(alpha=1)
        rdirichlet_1c <- function() {
            x <- rgamma(num_candidates, shape = 1)
            x / sum(x)
        }
        p <- matrix(0, nrow = num_groups, ncol = num_candidates)
        for (g_idx in seq_len(num_groups)) {
            p[g_idx, ] <- rdirichlet_1c()
        }
    }

    # Build X (b x c) by aggregating multinomial draws per group
    z_bgc <- array(
        0L,
        dim = c(num_groups, num_candidates, num_ballots),
        dimnames = list(
            group = seq_len(num_groups),
            candidate = seq_len(num_candidates),
            `ballot box` = seq_len(num_ballots)
        )
    )
    X <- matrix(0, nrow = num_ballots, ncol = num_candidates)
    for (b_idx in seq_len(num_ballots)) {
        for (g_idx in seq_len(num_groups)) {
            w_bg <- W[b_idx, g_idx]
            if (w_bg > 0) {
                # Draw from Multinomial(w_bg, p[g_idx, ])
                z_draw <- rmultinom(1, size = w_bg, prob = p[g_idx, ])
                z_bgc[g_idx, , b_idx] <- as.integer(z_draw)
                # Acumulate in X
                X[b_idx, ] <- X[b_idx, ] + z_draw
            }
        }
    }

    toReturn <- list(
        W = W,
        X = X,
        real_prob = p,
        outcome = z_bgc
    )

    # Return it as an eim object
    class(toReturn) <- "eim"

    return(toReturn)
}
