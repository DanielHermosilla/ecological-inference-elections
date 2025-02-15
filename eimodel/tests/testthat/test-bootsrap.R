# ========= Bootstrap function test ========= #
# For avoiding long precomputations, we'll use a small sample => reduce ballots boxes
tests <- 50
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
    model$bootstrap(bootstrap_iterations = 50, ballot_boxes = 20)

    # We should expect a visible matrix on the standard deviation attribute
    expect_type(model$standard_deviation, "matrix")
    expect_equal(nrow(model$standard_deviation), ncol(model$W))
    expect_equal(ncol(model$standard_deviation), ncol(model$X))

    # Unusual cases
    expect_error(model$bootstrap())
    expect_error(model$bootstrap(bootstrap_iterations = 8, ballot_boxes = nrow(result$X)))
    expect_invisible(model$bootstrap(main_method = "Hit and Run", step_size = 100, samples = 100))
    expect_invisible(model$bootstrap(bootstrap_iterations = 4, ballot_boxes = 5, main_method = "MVN CDF", method = "Genz"))

    # We should also expect that the other attributes aren't updated
    expect_equal(model$method, NULL)
    expect_equal(model$probability, NULL)
    expect_error(model$samples)
    expect_error(model$step_size)
    rm(model)
    gc()
  }
})
