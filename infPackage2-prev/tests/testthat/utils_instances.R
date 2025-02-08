library(testthat)

# Helper function for simulating random voting instances.
simulate_votation <- function(num_ballots, candidate_alpha, num_demographics,
                              demographic_alpha = rep(1, num_demographics),
                              total_votes = 100) {
    # If total_votes is a single number, make it a vector (one total per ballot box)
    if (length(total_votes) == 1) {
        total_votes <- rep(total_votes, num_ballots)
    } else if (length(total_votes) != num_ballots) {
        stop("total_votes must be either a single number or a vector of length num_ballots")
    }

    num_candidates <- length(candidate_alpha)

    # Preallocate matrices:
    # X will be (ballot boxes x candidates) -- remains unchanged.
    X <- matrix(0, nrow = num_ballots, ncol = num_candidates)
    # W will now be (ballot boxes x demographic groups)
    W <- matrix(0, nrow = num_ballots, ncol = num_demographics)

    # Function to generate a Dirichlet draw.
    # It draws 'n' samples, each of length = length(alpha)
    rdirichlet <- function(n, alpha) {
        k <- length(alpha)
        samples <- matrix(rgamma(n * k, shape = alpha, rate = 1), ncol = k, byrow = TRUE)
        samples / rowSums(samples)
    }

    # Loop over ballot boxes
    for (i in 1:num_ballots) {
        # --- Candidate votes (X) ---
        # Draw candidate vote shares for ballot box i:
        p <- rdirichlet(1, candidate_alpha)
        # Allocate the total_votes[i] among candidates:
        votes_candidates <- as.vector(rmultinom(1, size = total_votes[i], prob = p))
        X[i, ] <- votes_candidates

        # --- Demographic votes (W) ---
        # Draw demographic vote shares for ballot box i:
        q <- rdirichlet(1, demographic_alpha)
        # Allocate the same total votes among demographic groups:
        votes_demographics <- as.vector(rmultinom(1, size = total_votes[i], prob = q))
        W[i, ] <- votes_demographics
    }

    # At this point, for each ballot box i:
    #   sum(X[i, ]) == total_votes[i]  and  sum(W[i, ]) == total_votes[i].
    # That is, the candidate vote totals (by ballot box) equal the demographic vote totals.

    return(list(X = X, W = W))
}

random_samples <- function(ballots_range, # Arguments must be vectors of size 2
                           candidates_range,
                           demographic_range,
                           voting_range) {
    # Randomly draw the number of ballots, candidates and dirichlet parameters
    # Randomly choose a ballot box between 5 and 20
    num_ballots <- sample(ballots_range[1]:ballots_range[2], 1)
    # Randomly choose a candidate between 1 and 5
    num_candidates <- sample(candidates_range[1]:candidates_range[2], 1)
    # Randomly choose dirichlet parameters for each candidate between (0.1, 3)
    candidate_alpha <- runif(num_candidates, 0.1, 3)
    # Randomly choose demographic groups between 1 and 5
    num_demographics <- sample(demographic_range[1]:demographic_range[2], 1)
    # Randomly choose dirichlet parameters for each demographic group between (0.5, 3)
    demographic_alpha <- runif(num_demographics, 0.1, 5)
    # Randomly choose the total amount of votes per ballot box
    total_votes <- sample(voting_range[1]:voting_range[2], 1)

    choosen_values <- list(
        ballots = num_ballots,
        candidates = num_candidates,
        candidates_param = candidate_alpha,
        groups = num_demographics,
        groups_param = demographic_alpha,
        total_votes = total_votes
    )

    result <- simulate_votation(
        num_ballots = num_ballots,
        candidate_alpha = candidate_alpha,
        num_demographics = num_demographics,
        demographic_alpha = demographic_alpha,
        total_votes = total_votes
    )

    appended_list <- c(result, choosen_values)
    appended_list
}

# ========= Sampling simulation test ========= #
tests <- 10 # Number of random tests to run
test_that("The sampling function works", {
    for (i in 1:tests) {
        # Randomly draw parameters for the simulation:

        result <- random_samples(
            c(5, 50), # Ballots range
            c(1, 5), # Candidates range
            c(1, 5), # Group range
            c(50, 150) # Total votes range
        )

        # Check that for each ballot box, candidate votes (rows of X) match demographic votes (rows of W)
        expect_equal(rowSums(result$X), rowSums(result$W))
        # Check dimensions for candidates and demographic groups
        expect_equal(ncol(result$X), result$candidates)
        expect_equal(ncol(result$W), result$groups)
        # Check that the number of ballot boxes is consistent in both matrices
        expect_equal(nrow(result$X), result$ballots)
        expect_equal(nrow(result$W), result$ballots)
    }
})
