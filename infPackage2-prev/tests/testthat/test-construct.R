library(testthat)
library(infPackage2)
# Instance maker;
# simulate_votation <- function(num_ballots, candidate_alpha, num_demographics,
#                              demographic_alpha = rep(1, num_demographics),
#                              total_votes = 100)
source("utils_instances.R")
# ========= Constructor function test ========= #
tests <- 50
test_that("The constructor works", {
    # The "border case" would be to handle the case where there's no votes for a candidate
    for (i in 1:tests) {
        # Randomly draw the number of ballots, candidates and dirichlet parameters
        # Randomly choose a ballot box between 5 and 50
        num_ballots <- sample(5:50, 1)
        # Randomly choose a candidate between 1 and 5
        num_candidates <- sample(1:5, 1)
        # Randomly choose dirichlet parameters for each candidate between (0.1, 3)
        candidate_alpha <- runif(num_candidates, 0.1, 3)
        # Randomly choose demographic groups between 1 and 5
        num_demographics <- sample(2:5, 1)
        # Randomly choose dirichlet parameters for each demographic group between (0.5, 3)
        demographic_alpha <- runif(num_demographics, 0.1, 3)
        # Randomly choose the total amount of votes per ballot box
        total_votes <- sample(50:150, 1)
        # random_samples <- function(ballots_range,  # Arguments must be vectors of size 2
        # 						   candidates_range,
        # 						   demographic_range,
        # 						   voting_range){

        result <- random_samples(
            c(6, 200), # Ballots range
            c(1, 5), # Candidates range
            c(1, 5), # Demographic group range
            c(200, 250) # Votes per ballot range
        )

        # Construct the model
        constructing_model <- eim$new(result$X, result$W)
        # Check if it's equal to the X matrix
        expect_equal(constructing_model$X, result$X)
        # Check if it's equal to the W matrix
        expect_equal(constructing_model$W, result$W)
        # Check that it didn't create other attributes
        expect_equal(constructing_model$method, NULL)
        expect_equal(constructing_model$probability, NULL)
        expect_equal(constructing_model$logLikelihood, NULL)
        expect_equal(constructing_model$total_iterations, NULL)
        expect_equal(constructing_model$total_time, NULL)
        expect_equal(constructing_model$total_time, NULL)
        expect_error(constructing_model$samples)
        expect_error(constructing_model$step_size)
        expect_error(constructing_model$multivariate_method)
        expect_error(constructing_model$multivariate_error)
        expect_error(constructing_model$multivariate_iterations)
        # We should also expect an error if there's an incongruency
        # between the rows dimension between X and W

        # We won't get an error though if it's a square matrix, so we
        # impose that the ballots range > candidates, groups.
        wrong_matrix <- t(result$X)
        expect_error(wrong_matrix, result$W)
    }
})
# ========= Precomputation function test ========= #
# For avoiding long precomputations, we'll use a small sample => reduce ballots boxes
tests <- 10
test_that("Precomputation works", {
    for (i in 1:tests) {
        result <- random_samples(
            c(5, 20), # Ballots range
            c(1, 5), # Candidates range
            c(1, 5), # Demographic group range
            c(200, 250) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)

        candidates <- result$candidates
    }
})

