library(testthat)

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
