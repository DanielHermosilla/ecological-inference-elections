# ========= Computation function test ========= #
# For avoiding long computations, we'll use a small sample => reduce ballots boxes
# This may take a while...
tests <- 5
results <- vector("list", tests)
for (i in 1:tests) {
    results[[i]] <- fastei:::.random_samples(
        c(15, 25), # Ballots range
        c(2, 4), # Candidates range
        c(2, 4), # Demographic group range
        c(30, 40), # Votes per ballot range
        seed = i
    )
}
# For precision tests, we'll use big samples, except for the exact and hnr (whom time of computation would be very long)
precisionTests <- 50
precisionsResults <- vector("list", precisionTests)
for (i in 1:precisionTests) {
    precisionsResults[[i]] <- fastei:::.random_samples(
        c(200, 400), # Ballots range
        c(2, 10), # Candidates range
        c(2, 10), # Demographic group range
        c(200, 400), # Votes per ballot range
        seed = i
    )
}
available_methods <- c("exact", "mult", "hnr", "mvn_cdf", "mvn_pdf")
method_division <- tests %/% 5

test_that("run_em:\tProbability output is coherent", {
    for (k in 1:5) {
        for (i in 1:tests) {
            result <- run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = 10
            )
            expect_is(result$prob, "matrix")
            for (r in 1:ncol(results[[i]]$X)) { # nolint
                expect_equal(rowSums(result$prob), rep(1, nrow(result$prob)))
            }
            expect_equal(ncol(result$prob), ncol(results[[i]]$X))
            expect_equal(nrow(result$prob), ncol(results[[i]]$W))
            rm(result)
        }
    }
})

test_that("run_em:\tLog likelihood and iterations output is coherent", {
    for (k in 1:5) {
        for (i in 1:tests) {
            result <- run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = 10
            )
            expect_type(result$logLik, "double")
            expect_equal(length(result$logLik), result$iterations)
            rm(result)
        }
    }
})

test_that("run_em:\tStatus and messages output is coherent", {
    for (k in 1:5) {
        for (i in 1:tests) {
            result <- run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = 7, maxtime = 10
            )
            if (result$status == 0) {
                expect_equal(result$message, "Convergence achieved")
            } else if (result$status == 1) {
                expect_equal(result$message, "Log-likelihood decrease")
            } else if (result$status == 2) {
                expect_equal(result$message, "Maximum time reached")
            } else {
                expect_equal(result$message, "Maximum iterations reached")
            }
            rm(result)
        }
    }
})

test_that("run_em:\tComputation works for the exact method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- run_em(model, method = "exact")
        expect_equal(result$method, "exact")
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat(sprintf("\nPrecision result for exact is: %.6f\n", difference / tests))
    expect_lt(difference / tests, 0.1) # Less than 0.1 of error
})

test_that("run_em:\tComputation works for the Hit and Run method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:5) {
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- run_em(model, method = "hnr", maxtime = 10)
        expect_equal(result$method, "hnr")
        expect_equal(result$samples, 1000)
        expect_equal(result$step_size, 3000)
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat(sprintf("\nPrecision result for Hit and Run is: %.6f\n", difference / 5))
    expect_lt(difference / 5, 0.1) # Less than 0.1 of error
})

test_that("run_em:\tComputation works for the Hit and Run method with distinct step_size and sample size", {
    # Test with distinct parameters
    for (i in 1:5) {
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- run_em(model, method = "hnr", step_size = 1000, maxtime = 10)
        expect_equal(result$step_size, 1000)
        expect_equal(result$samples, 1000)

        result <- run_em(model, method = "hnr", samples = 800, maxtime = 10)
        expect_equal(result$samples, 800)
        expect_equal(result$step_size, 3000)

        result <- run_em(model, method = "hnr", step_size = 123, samples = 321, maxtime = 10)
        expect_equal(result$samples, 321)
        expect_equal(result$step_size, 123)
        rm(model)
        rm(result)
        gc()
    }
})

test_that("run_em:\tComputation works for the Multinomial method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:precisionTests) {
        model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
        result <- run_em(model, method = "mult", maxtime = 10)
        expect_equal(result$method, "mult")
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(result$prob - precisionsResults[[i]]$prob))
        rm(model)
        rm(result)
        gc()
    }
    cat(sprintf("\nPrecision result for Multinomial is: %.6f\n", difference / tests))
    expect_lt(difference / precisionTests, 0.1) # Less than 0.1 of error
})

test_that("run_em:\tComputation works for the MVN CDF method with 'Genz2' method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:precisionTests) {
        model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
        result <- run_em(model, method = "mvn_cdf", maxtime = 10)
        expect_equal(result$method, "mvn_cdf")
        expect_equal(result$mc_method, "genz2")
        expect_equal(result$mc_error, 1e-6)
        expect_equal(result$mc_samples, 5000)
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(result$prob - precisionsResults[[i]]$prob))
        rm(model)
        rm(result)
        gc()
    }
    cat(sprintf("\nPrecision result for MVN CDF is: %.6f\n", difference / precisionTests))
    expect_lt(difference / precisionTests, 0.1) # Less than 0.1 of error
})

test_that("run_em:\tComputation works for the MVN CDF method with 'Genz' method", {
    # Test with the "genz" method
    difference <- 0
    for (i in 1:precisionTests) {
        model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
        result <- run_em(model, method = "mvn_cdf", mc_method = "genz", maxtime = 10)
        expect_equal(result$method, "mvn_cdf")
        expect_equal(result$mc_method, "genz")
        expect_equal(result$mc_error, 1e-6)
        expect_equal(result$mc_samples, 5000)
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(result$prob - precisionsResults[[i]]$prob))
        rm(model)
        rm(result)
        gc()
    }
    cat(sprintf("\nPrecision result for MVN with Genz2 is: %.6f\n", difference / precisionTests))
    expect_lt(difference / precisionTests, 0.1) # Less than 0.1 of error
})

test_that("run_em:\tComputation works for the MVN CDF method with distinct iterations and sampling methods", {
    # Test with distinct parameters
    for (i in 1:tests) {
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- run_em(model, method = "mvn_cdf", mc_method = "genz2", mc_samples = 3000, maxtime = 10)
        expect_equal(result$mc_samples, 3000)
        expect_equal(result$mc_error, 1e-6)
        expect_is(result$prob, "matrix")

        result <- run_em(model, method = "mvn_cdf", mc_method = "genz2", mc_error = 1e-8, maxtime = 10)
        expect_equal(result$mc_samples, 5000)
        expect_equal(result$mc_error, 1e-8)
        expect_is(result$prob, "matrix")

        result <- run_em(model, method = "mvn_cdf", mc_method = "genz", mc_error = 1e-8, mc_samples = 10000, maxtime = 10)
        expect_equal(result$mc_error, 1e-8)
        expect_equal(result$mc_samples, 10000)
        expect_is(result$prob, "matrix")
        rm(model)
        rm(result)
        gc()
    }
})

test_that("run_em:\tComputation works for the MVN PDF method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:precisionTests) {
        model <- eim(precisionsResults[[i]]$X, precisionsResults[[i]]$W)
        result <- run_em(model, method = "mvn_pdf", maxtime = 10)
        expect_equal(result$method, "mvn_pdf")
        expect_is(result$prob, "matrix")
        difference <- difference + mean(abs(model$prob - precisionsResults[[i]]$prob))
        rm(model)
        rm(result)
        gc()
    }
    cat(sprintf("\nPrecision result for MVN PDF is: %.6f\n", difference / precisionTests))
    expect_lt(difference / precisionTests, 0.1) # Less than 0.1 of error
})

test_that("run_em: Verbose works correctly.", {
    for (k in 1:5) {
        for (i in 1:tests) {
            expect_output(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], verbose = TRUE, maxtime = 3
            ))
        }
    }
})

test_that("run_em: Stop threshold correctly.", {
    for (k in 1:5) {
        for (i in 1:tests) {
            x <- runif(1, min = .Machine$double.eps, max = 1)
            expect_no_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], stop_threshold = x, maxtime = 10
            ))
        }
    }
})

test_that("run_em: Errors when an invalid method is provided", {
    for (i in 1:3) {
        expect_error(run_em(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "Hit and Run", maxtime = 10
        ))
        expect_error(run_em(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "HNR", maxtime = 10
        ))
        expect_error(run_em(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "EXakt", maxtime = 10
        ))
    }
})

test_that("run_em: Errors when an object and matrices are provided", {
    for (i in 1:tests) {
        expect_error(run_em(
            X = results[[i]]$X, W = results[[i]]$W,
            object = eim(X = results[[i]]$X, W = results[[i]]$W), maxtime = 10
        ))
        dummyobj <- eim(X = results[[i]]$X, W = results[[i]]$W)
        expect_error(run_em(
            X = results[[i]]$X, W = results[[i]]$W,
            object = dummyobj
        ))
        expect_error(run_em(X = results[[i]]$X, object = dummyobj))
        expect_error(run_em(W = results[[i]]$W, object = dummyobj))
    }
})

test_that("run_em: Errors when an invalid json is provided", {
    expect_error(run_em(json_path = "invalid.json"))
})

test_that("run_em: Errors when an invalid iteration parameter is provided", {
    for (k in 1:5) {
        for (i in 1:tests) {
            x <- runif(1, min = .Machine$double.eps, max = 1)
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = x, maxtime = 10
            ))
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = 0, maxtime = 10
            ))
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = -3, maxtime = 10
            ))
        }
    }
})

test_that("run_em: Errors when an invalid time parameter is provided", {
    for (k in 1:5) {
        for (i in 1:tests) {
            x <- runif(1, min = .Machine$double.eps, max = 1)
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = -3
            ))
        }
    }
})

test_that("run_em: Errors when an invalid stop threshold is provided", {
    for (k in 1:5) {
        for (i in 1:tests) {
            x <- runif(1, min = .Machine$double.eps, max = 1)
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], stop_threshold = -3
            ))
        }
    }
})

test_that("run_em: Errors when an invalid verbose is provided", {
    for (k in 1:5) {
        for (i in 1:tests) {
            x <- runif(1, min = .Machine$double.eps, max = 1)
            expect_error(run_em(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], verbose = 0.1
            ))
        }
    }
})

test_that("run_em: Multiple samples runs without any error", {
    for (i in 1:1000) {
        sample <- fastei:::.random_samples(
            c(400, 500), # Ballots range
            c(2, 10), # Candidates range
            c(2, 10), # Demographic group range
            c(200, 500), # Votes per ballot range
            seed = i + tests
        )
        expect_no_error(run_em(
            X = sample$X, W = sample$W, maxtime = 60
        ))
    }
})
