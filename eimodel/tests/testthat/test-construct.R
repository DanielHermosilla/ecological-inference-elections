library(testthat)
library(infPackage2)
# Instance maker;
# simulate_votation <- function(num_ballots, candidate_alpha, num_demographics,
#                              demographic_alpha = rep(1, num_demographics),
#                              total_votes = 105)
source("utils_instances.R")
# ========= Constructor function test ========= #
tests <- 50
test_that("The constructor works with random matrices", {
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
    result <- random_samples(
      c(6, 200), # Ballots range
      c(2, 5), # Candidates range
      c(2, 5), # Demographic group range
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
    rm(constructing_model)
  }
  gc()
})

test_that("The constructor doesn't work with invalid files paths", {
  expect_error(eim$new("invalid.json"))
  expect_error(eim$new("NA."))
  expect_error(eim$new("n"))
})

test_that("The constructor doesn't work with matrices with different row size.", {
  # We will create matrices such that their rows aren't equal.
  for (i in 1:tests) {
    matrix1 <- matrix(0, nrow = sample(2:10, 1), ncol = 4)
    matrix2 <- matrix(0, nrow = sample(2:10, 1), ncol = 3)

    # We don't want to check the case where their rows are equal
    if (nrow(matrix1) == nrow(matrix2)) next

    expect_error(eim$new(matrix1, matrix2))
    expect_error(eim$new(matrix2, matrix1))
  }
  gc()
})

test_that("The constructor doesn't work with matrices of dimension < 1", {
  # We will create matrices such that their rows aren't equal.
  matrix1 <- matrix(0, nrow = 2, ncol = 1)
  valid <- matrix(0, nrow = 2, ncol = 2)

  expect_error(eim$new(matrix1, valid))
  expect_error(eim$new(valid, matrix1))
})

# ========= Precomputation function test ========= #
# For avoiding long precomputations, we'll use a small sample => reduce ballots boxes
tests <- 10
test_that("Precomputation works", {
  for (i in 1:tests) {
    result <- random_samples(
      c(45, 60), # Ballots range
      c(2, 5), # Candidates range
      c(2, 5), # Demographic group range
      c(100, 150) # Votes per ballot range
    )

    # From the test before, we can assume the constructing function work
    model <- eim$new(result$X, result$W)
    step <- sample(500:1000, 1)
    samples <- sample(100:200, 1)
    candidates <- result$candidates
    groups <- result$groups
    expect_invisible(model$precompute("Exact"))
    expect_invisible(model$precompute("Hit and Run", step_size = step, samples = samples))
    expect_error(model$precompute("Hit and Run"))
    expect_error(model$precompute("Hit and Run", step_size = 3.14, samples = samples))
    expect_error(model$precompute("Hit and Run", step_size = step, samples = 3.14))
    expect_error(model$precompute("Multinomial"))

    # We should also expect that the other attributes aren't updated
    expect_equal(model$method, NULL)
    expect_equal(model$probability, NULL)
    expect_error(model$samples)
    expect_error(model$step_size)
    rm(model)
    gc()
  }
})
