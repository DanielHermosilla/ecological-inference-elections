#' Simulate an Election
#'
#' @description
#' This function simulates an election by creating matrices representing candidate votes
#' (`X`) and demographic group votes (`W`) across a specified number of ballot boxes.
#' It either receives (or accepts) a probability matrix (`prob`) indicating how likely
#' each demographic group is to vote for each candidate.
#'
#' By default, the number of voters per ballot box (`ballot_voters`) is set to a vector
#' of 50s with length `num_ballots`. You can optionally override this by providing a
#' custom vector.
#'
#' Several optional parameters are available to control the distribution of votes:
#' \itemize{
#'   \item \strong{`group_proportions`}: A vector of length `num_groups` specifying
#'         the overall proportion of each demographic group. Entries must sum to one and be non-negative.
#'   \item \strong{`prob`}: A user-supplied probability matrix of dimension
#'         (`num_groups` \eqn{\times} `num_candidates`). If provided, this matrix is used directly
#'         instead of drawing from a Dirichlet distribution.
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
#'   box. Defaults to `rep(50, num_ballots)`.
#'
#' @param lambda
#' 	A numeric value between `0` and `1` that controls the **shuffling proportion**
#'  of individuals among demographic groups (`num_groups`). This parameter determines
#'  how much **random reassignment** occurs among groups, affecting group heterogeneity.
#'
#'   - **`lambda = 0`**: No shuffling occurs. Individuals remain in their original
#'     assigned groups, leading to **maximum heterogeneity**.
#'   - **`lambda = 1`**: Complete shuffling. Individuals are **fully reassigned**,
#'     leading to **homogeneous groups**.
#'   - **Intermediate values** (e.g., `lambda = 0.5`): A fraction of individuals is **randomly selected and reassigned**, balancing between heterogeneity and homogeneity.
#'
#'  Default value is set to `0.5`.
#'
#' @param seed
#'   If provided, overrides the current global seed. Defaults to `NULL`.
#'
#' @param group_proportions
#'   Optional. A vector specifying the overall proportion of each group among all of the ballot boxes. Defaults to an uniform distribution.
#'
#'
#' @param prob
#'   Optional. A user-supplied probability matrix of dimension `(g x c)`.
#'   If provided, this matrix is used directly instead of
#'   drawing from a Dirichlet distribution with each parameter set to 1 (yielding an uninformative prior).
#'
#' @section Shuffling Mechanism:
#' The shuffling step is controlled by the `lambda` (\eqn{\lambda}) parameter and operates as follows:
#'
#' 1. **Initial Group Assignment**: Each of the individuals (`I`) is evenly assigned to one of the `num_groups`:
#' 		\deqn{\omega_{i}^{0} = \lceil i \cdot \lvert \text{num\_groups}\rvert \cdot I^{-1} \rceil}
#'    This ensures all groups initially contain an equal count of individuals. If `group_proportions` isn't uniform, it assigns \eqn{\omega_{i}^{0}} according the given parameters.
#'
#' 2. **Determining the Shuffling Fraction**: The fraction of individuals to shuffle is given by \eqn{\lambda \cdot I}. Hence, different `lambda` values are interpreted as follows:
#'
#' 	- `lambda = 0` means no one is shuffled (groups remain as is).
#' 	- `lambda = 1` means all individuals are shuffled.
#' 	- Intermediate values like `lambda = 0.5` shuffle half the population.
#'
#' 3. **Random Selection of Individuals**: From the total population, \eqn{\lambda\cdot I} individuals are randomly sampled without replacement into a set `v`. This set dictates which individuals
#'    will have their group assignments swapped.
#'
#' 4. **Sorting and Swapping**: The selected indices in `v` are sorted into a new vector `v_sorted`. Then, each individual `v[i]` obtains the group assignment of `v_sorted[i]`. This approach preserves the total number of people in each group while introducing a new arrangement of individuals.
#'
#' @return A list with three components:
#' \describe{
#'   \item{\code{W}}{A \code{(b x g)} matrix of demographic votes per ballot box.}
#'   \item{\code{X}}{A \code{(b x c)} matrix of candidate votes per ballot box.}
#'   \item{\code{prob}}{A \code{(g x c)} matrix of probabilities that each group votes for each candidate. If `prob` is provided, it would equal such probability.}
#' }
#'
#' @references
#' The algorithm is fully explained in ['Thraves, C. and Ubilla, P.: *"Fast Ecological Inference Algorithm for the R×C Case"*](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834).
#'
#' @examples
#' # Example 1: Default usage with 2 ballot boxes, each having 50 voters
#' result1 <- simulate_election(
#'     num_ballots = 55,
#'     num_candidates = 3,
#'     num_groups = 5
#' )
#'
#' # Example 2: Using a custom ballot_voters vector
#' result2 <- simulate_election(
#'     num_ballots = 80,
#'     num_candidates = 2,
#'     num_groups = 7,
#'     ballot_voters = c(100, 50, 75)
#' )
#'
#' # Example 3: Supplying group_proportions
#' result3 <- simulate_election(
#'     num_ballots = 40,
#'     num_candidates = 7,
#'     num_groups = 4,
#'     group_proportions = c(0.3, 0.5, 0.1, 0.1)
#' )
#'
#' # Example 4: Providing a user-defined prob matrix
#' custom_prob <- matrix(c(
#'     0.9, 0.1,
#'     0.4, 0.6
#' ), nrow = 2, byrow = TRUE)
#'
#' result4 <- simulate_election(
#'     num_ballots = 2,
#'     num_candidates = 2,
#'     num_groups = 2,
#'     prob = custom_prob
#' )
#'
#' result4$prob == custom_prob # TRUE
#' @export
simulate_election <- function(num_ballots,
                              num_candidates,
                              num_groups,
                              ballot_voters = rep(50, num_ballots),
                              lambda = 0.5,
                              seed = NULL,
                              group_proportions = rep(1 / num_groups, num_groups),
                              prob = NULL) {
    # If user provides a seed, override the current global seed
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Validate length of ballot_voters
    if (length(ballot_voters) != num_ballots) {
        stop("`ballot_voters` must be a vector of length `num_ballots`.")
    }

    # Validate group proportions
    if (length(group_proportions) != num_groups) {
        stop("group_proportions must have length equal to num_groups.")
    }
    tol <- .Machine$double.eps^0.5
    if (abs(sum(group_proportions) - 1) > tol) {
        stop("group_proportions must sum to 1.")
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
    X <- matrix(0, nrow = num_ballots, ncol = num_candidates)
    for (b_idx in seq_len(num_ballots)) {
        for (g_idx in seq_len(num_groups)) {
            w_bg <- W[b_idx, g_idx]
            if (w_bg > 0) {
                # Draw from Multinomial(w_bg, p[g_idx, ])
                z_bgc <- rmultinom(1, size = w_bg, prob = p[g_idx, ])
                X[b_idx, ] <- X[b_idx, ] + z_bgc
            }
        }
    }

    # Return the list
    list(
        W = W,
        X = X,
        prob = p
    )
}
