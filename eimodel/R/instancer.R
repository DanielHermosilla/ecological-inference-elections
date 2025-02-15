#' Simulates an election by creating the candidate and group matrices and their results.
#'
#' @description
#' Given the number of ballots, groups, candidates, and votes per ballot, it will simulate an election.
#' Additionally, it generates the resulting matrix (of dimension g x c) that represents the unobserved probabilities
#' that a demographic group votes for a given candidate. These probabilities are drawn from a Dirichlet distribution
#' and a \code{lambda} value, which represents the heterogeneity of the groups.
#'
#' @param num_ballots \emph{(integer(1))} The number of ballot boxes \emph{("b")}.
#'
#' @param num_candidates \emph{(integer(1))} The number of candidates \emph{("c")}.
#'
#' @param num_groups \emph{(integer(1))} The number of demographic groups \emph{("g")}.
#'
#' @param ballot_voters \emph{(integer(num_ballots))} A vector of size \code{num_ballots} with the number of votes per ballot box.
#'
#' @param lambda \emph{(numeric(1))} A value between 0 and 1 representing the group heterogeneity. Values near 0 are more heterogeneous,
#' but less realistic. The same goes the other way around. (default: 0.5)
#'
#' @param seed \emph{(numeric(1)} Optional. If provided, it overrides the current global seed. (default: \code{NULL})
#'
#' @note A reminder of the Dirichlet distribution is that alpha = 1 represents a distribution that is
#' uniform around the mean. On the other hand, alpha < 1 tends to produce a sparser probability vector, with many
#' values close to zero and one or a few larger values. The inverse case occurs when alpha > 1, tightening
#' the concentration around the mean.
#'
#' @seealso [random_samples()]
#'
#' @return A list with components:
#' \item{X}{A matrix (b x c) with candidate votes per ballot box.}
#' \item{W}{A matrix (b x g) with demographic votes per ballot box.}
#' \item{real_p}{A matrix (g x c) with the estimated (unobserved) probabilities that a demographic group votes for a given candidate.}
#'
#' @examples
#' result <- simulate_elections(num_ballots = 10, num_candidates = 5, num_groups = 3, ballot_voters = c(100, 10))
#' result$X # Candidate matrix (b x c)
#' result$W # Group matrix (b x g)
#' result$real_p # Probability matrix (g x c)
#'
#' @export
simulate_elections <- function(num_ballots, # Number of ballot boxes
                               num_candidates, # Number of candidates
                               num_groups, # Number of demographic groups
                               ballot_voters, # Vector of length "ballot boxes" with the number of voters per box
                               lambda = 0.5, # Proportion of individuals to "shuffle"
                               seed = NULL) {
  # If user provides a seed, override the current global seed:
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(ballot_voters) != num_ballots) {
    stop("`ballot_voters` must be a vector of length `num_ballots`.")
  }

  # Form the total population size:
  n <- sum(ballot_voters)

  # Split the population roughly evenly among the G groups.
  omega0 <- ceiling((1:n) * num_groups / n)

  # Get a copy of the vector for the shuffle
  omega <- omega0

  # Sample λ·I individuals v (without replacement) and shuffle them:
  size_shuffle <- floor(lambda * n)
  if (size_shuffle > 0) {
    v <- sample.int(n, size_shuffle, replace = FALSE) # the λ·I indices
    v_ <- sort(v) # sorted
    for (k in seq_len(size_shuffle)) {
      omega[v[k]] <- omega0[v_[k]]
    }
  }

  # Get the cumulative sum for later creating the partition between ballot boxes
  # This will allow defining a range of voters between values of a[]. For
  # example: a[2] - a[1] = voters of ballot box 1.
  a <- c(0, cumsum(ballot_voters))

  # Build W
  w <- matrix(0, nrow = num_ballots, ncol = num_groups)
  for (b_idx in seq_len(num_ballots)) {
    # Indices of the population that fall into box b
    start_idx <- a[b_idx] + 1
    end_idx <- a[b_idx + 1]
    # Tally the groups in this segment:
    sub_omega <- omega[start_idx:end_idx]
    # Count how many are in each group:
    # (Note that it's shuffled according to lambda)
    counts_bg <- tabulate(sub_omega, nbins = num_groups)
    w[b_idx, ] <- counts_bg
  }

  # Define the probability matrix
  # small helper for Dirichlet(1C) draws:
  rdirichlet_1c <- function() {
    # shape = rep(1, C)
    x <- rgamma(num_candidates, shape = 1)
    x / sum(x)
  }

  p <- matrix(0, nrow = num_groups, ncol = num_candidates)
  for (g_idx in seq_len(num_groups)) {
    p[g_idx, ] <- rdirichlet_1c()
  }

  # For each (b,g), draw z_{bg} (C-dimensional) from Multinomial(w_{bg}, p_g).
  # Then x_{bc} = sum_{g} z_{bgc}.
  x <- matrix(0, nrow = num_ballots, ncol = num_candidates)

  for (b_idx in seq_len(num_ballots)) {
    for (g_idx in seq_len(num_groups)) {
      w_bg <- w[b_idx, g_idx] # number of people of group g in box b
      if (w_bg > 0) {
        # draw from multinomial
        z_bgc <- as.vector(rmultinom(1, size = w_bg, prob = p[g_idx, ]))
        x[b_idx, ] <- x[b_idx, ] + z_bgc
      }
    }
  }

  list(
    W = w,
    X = x,
    real_p = p
  )
}

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
#' @seealso [simulate_elections()]
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
#' @export
random_samples <- function(ballots_range, # Arguments must be vectors of size 2
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

  choosen_values <- list(
    ballots = num_ballots,
    candidates = num_candidates,
    groups = num_groups,
    total_votes = total_votes
  )

  result <- simulate_elections(
    num_ballots = num_ballots,
    num_candidates = num_candidates,
    num_groups = num_groups,
    ballot_voters = total_votes,
    seed = seed
  )

  appended_list <- c(result, choosen_values)
  appended_list
}
