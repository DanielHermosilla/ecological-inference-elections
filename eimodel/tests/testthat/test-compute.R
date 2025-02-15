# ========= Computation function test ========= #
# For avoiding long computations, we'll use a small sample => reduce ballots boxes
# This may take a while...
tests <- 10
test_that("Computation works for the exact method", {
    difference <- 0
    for (i in 1:tests) {
        result <- random_samples(
            c(5, 20), # Ballots range
            c(2, 3), # Candidates range
            c(2, 3), # Demographic group range
            c(20, 100) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)

        candidates <- result$candidates
        groups <- result$groups

        # Exact method
        message("\nDEBUG:\tRunning test with ", candidates, " candidates, ", groups, " groups, and method: EXACT")
        expect_invisible(model$compute("Exact"))
        expect_equal(model$method, "Exact")
        expect_gt(model$total_time, 0)
        expect_gt(model$total_iterations, 0)
        expect_true(is.integer(model$total_iterations))
        expect_true(is.matrix(model$probability))
        expect_true((is.numeric(model$logLikelihood) && is.atomic(model$logLikelihood)))
        difference <- difference + mean(abs(model$probability - result$ground_truth))

        rm(model)
        gc()
    }
    expect_true(difference / tests < 0.1)
})

test_that("Computation works for the Hit and Run method", {
    difference <- 0
    for (i in 1:tests) {
        result <- random_samples(
            c(45, 60), # Ballots range
            c(2, 4), # Candidates range
            c(2, 4), # Demographic group range
            c(80, 120) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)
        step <- sample(2000:4000, 1)
        samples <- sample(1000:1500, 1)

        candidates <- result$candidates
        groups <- result$groups
        # Hit and Run method
        message("\nDEBUG:\tRunning test with ", candidates, " candidates, ", groups, " groups, and method: HIT AND RUN")
        expect_invisible(model$compute("Hit and Run", step_size = step, samples = samples))
        expect_equal(model$method, "Hit and Run")
        expect_equal(model$samples, samples)
        expect_equal(model$step_size, step)
        expect_gte(model$total_time, 0)
        expect_gte(model$total_iterations, 0)
        expect_true(is.integer(model$total_iterations))
        expect_true(is.matrix(model$probability))
        expect_true((is.numeric(model$logLikelihood) && is.atomic(model$logLikelihood)))
        difference <- difference + mean(abs(model$probability - result$ground_truth))


        rm(model)
        gc()
    }
    expect_true(difference / tests < 0.1)
})
test_that("Computation works for the MVN CDF method", {
    difference <- 0
    for (i in 1:tests) {
        result <- random_samples(
            c(45, 60), # Ballots range
            c(2, 4), # Candidates range
            c(2, 4), # Demographic group range
            c(100, 120) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)
        step <- sample(2000:4000, 1)
        samples <- sample(900:1100, 1)

        candidates <- result$candidates
        groups <- result$groups
        # MVN CDF method
        message("\nDEBUG:\tRunning test with ", candidates, " candidates, ", groups, " groups, and method: MVN CDF")
        expect_invisible(model$compute("MVN CDF"))
        expect_equal(model$method, "MVN CDF")
        expect_equal(model$multivariate_method, "Genz2")
        expect_equal(model$multivariate_error, 1e-5)
        expect_equal(model$multivariate_iterations, 5000)
        expect_gte(model$total_time, 0)
        expect_gte(model$total_iterations, 0)
        expect_true(is.integer(model$total_iterations))
        expect_true(is.matrix(model$probability))
        expect_true((is.numeric(model$logLikelihood) && is.atomic(model$logLikelihood)))
        difference <- difference + mean(abs(model$probability - result$ground_truth))


        rm(model)
        gc()
    }
    expect_true(difference / tests < 0.1)
})

test_that("Computation works for the MVN PDF method", {
    difference <- 0
    for (i in 1:tests) {
        result <- random_samples(
            c(100, 120), # Ballots range
            c(2, 4), # Candidates range
            c(2, 4), # Demographic group range
            c(80, 120) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)
        step <- sample(500:1000, 1)
        samples <- sample(100:200, 1)

        candidates <- result$candidates
        groups <- result$groups
        # MVN PDF method
        message("\nDEBUG:\tRunning test with ", candidates, " candidates, ", groups, " groups, and method: MVN PDF")
        expect_invisible(model$compute("MVN PDF"))
        expect_equal(model$method, "MVN PDF")
        expect_gte(model$total_time, 0)
        expect_gte(model$total_iterations, 0)
        expect_true(is.integer(model$total_iterations))
        expect_true(is.matrix(model$probability))
        expect_true((is.numeric(model$logLikelihood) && is.atomic(model$logLikelihood)))
        difference <- difference + mean(abs(model$probability - result$ground_truth))

        rm(model)
        gc()
    }
    expect_true(difference / tests < 0.1)
})

test_that("Computation works for the Multinomial method", {
    difference <- 0
    for (i in 1:tests) {
        result <- random_samples(
            c(100, 200), # Ballots range
            c(2, 10), # Candidates range
            c(2, 10), # Demographic group range
            c(80, 120) # Votes per ballot range
        )

        # From the test before, we can assume the constructing function work
        model <- eim$new(result$X, result$W)
        step <- sample(500:1000, 1)
        samples <- sample(100:200, 1)

        candidates <- result$candidates
        groups <- result$groups
        # Multinomial method
        message("\nDEBUG:\tRunning test with ", candidates, " candidates, ", groups, " groups, and method: Multinomial")
        expect_invisible(model$compute("Multinomial"))
        expect_equal(model$method, "Multinomial")
        expect_gte(model$total_time, 0)
        expect_gte(model$total_iterations, 0)
        expect_true(is.integer(model$total_iterations))
        expect_true(is.matrix(model$probability))
        expect_true((is.numeric(model$logLikelihood) && is.atomic(model$logLikelihood)))
        difference <- difference + mean(abs(model$probability - result$ground_truth))

        rm(model)
        gc()
    }
    expect_true(difference / tests < 0.1)
})
