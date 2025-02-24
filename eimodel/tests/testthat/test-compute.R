# ========= Computation function test ========= #
# For avoiding long computations, we'll use a small sample => reduce ballots boxes
# This may take a while...
tests <- 10
results <- vector("list", tests)
for (i in 1:tests) {
    cat(sprintf("Generating sample %d...\n", i))
    results[[i]] <- eimodel:::.random_samples(
        c(30, 50), # Ballots range
        c(2, 4), # Candidates range
        c(2, 4), # Demographic group range
        c(30, 60), # Votes per ballot range
        seed = i
    )
}

available_methods <- c("exact", "mult", "hnr", "mvn_cdf", "mvn_pdf")
method_division <- tests %/% 5

test_that("Compute:\tProbability output is coherent", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'prob output is coherent' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            result <- compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = 10
            )
            cat("\nAbout to run expect_type(result$prob, 'double')\n")
            expect_type(result$prob, "double")
            for (r in 1:ncol(results[[i]]$X)) { # nolint
                cat(sprintf("\nAbout to run expect_equal(rowSums(result$prob), rep(1, nrow(result$prob))) for candidate %d\n", r))
                expect_equal(rowSums(result$prob), rep(1, nrow(result$prob)))
            }
            cat("\nAbout to run expect_equal(ncol(result$prob), ncol(results[[i]]$X))\n")
            expect_equal(ncol(result$prob), ncol(results[[i]]$X))
            cat("\nAbout to run expect_equal(nrow(result$prob), ncol(results[[i]]$W))\n")
            expect_equal(nrow(result$prob), ncol(results[[i]]$W))
            rm(result)
        }
    }
})

test_that("Compute:\tLog likelihood and iterations output is coherent", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Log likelihood and iterations' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            result <- compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = 10
            )
            cat("\nAbout to run expect_type(result$logLik, 'double')\n")
            expect_type(result$logLik, "double")
            cat("\nAbout to run expect_equal(length(result$logLik), result$iterations)\n")
            expect_equal(length(result$logLik), result$iterations)
            rm(result)
        }
    }
})

test_that("Compute:\tStatus and messages output is coherent", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Status and messages' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            result <- compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = 7, maxtime = 10
            )
            if (result$status == 0) {
                cat("\nAbout to run expect_equal(result$message, 'Convergence achieved')\n")
                expect_equal(result$message, "Convergence achieved")
            } else if (result$status == 1) {
                cat("\nAbout to run expect_equal(result$message, 'Log-likelihood decrease')\n")
                expect_equal(result$message, "Log-likelihood decrease")
            } else if (result$status == 2) {
                cat("\nAbout to run expect_equal(result$message, 'Maximum time reached')\n")
                expect_equal(result$message, "Maximum time reached")
            } else {
                cat("\nAbout to run expect_equal(result$message, 'Maximum iterations reached')\n")
                expect_equal(result$message, "Maximum iterations reached")
            }
            rm(result)
        }
    }
})

test_that("Compute:\tComputation works for the exact method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'exact method' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "exact")
        cat("\nAbout to run expect_equal(result$method, 'exact')\n")
        expect_equal(result$method, "exact")
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tComputation works for the Hit and Run method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'Hit and Run method' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "hnr", maxtime = 10)
        cat("\nAbout to run expect_equal(result$method, 'hnr')\n")
        expect_equal(result$method, "hnr")
        cat("\nAbout to run expect_equal(result$samples, 1000)\n")
        expect_equal(result$samples, 1000)
        cat("\nAbout to run expect_equal(result$step_size, 3000)\n")
        expect_equal(result$step_size, 3000)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tComputation works for the Hit and Run method with distinct step_size and sample size", {
    # Test with distinct parameters
    for (i in 1:tests) {
        cat(sprintf("Running 'Hit and Run distinct parameters' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "hnr", step_size = 1000, maxtime = 10)
        cat("\nAbout to run expect_equal(result$step_size, 1000)\n")
        expect_equal(result$step_size, 1000)
        cat("\nAbout to run expect_equal(result$samples, 1000)\n")
        expect_equal(result$samples, 1000)

        result <- compute(model, method = "hnr", samples = 800, maxtime = 10)
        cat("\nAbout to run expect_equal(result$samples, 800)\n")
        expect_equal(result$samples, 800)
        cat("\nAbout to run expect_equal(result$step_size, 3000)\n")
        expect_equal(result$step_size, 3000)

        result <- compute(model, method = "hnr", step_size = 123, samples = 321, maxtime = 10)
        cat("\nAbout to run expect_equal(result$samples, 321)\n")
        expect_equal(result$samples, 321)
        cat("\nAbout to run expect_equal(result$step_size, 321)\n")
        expect_equal(result$step_size, 321)
        rm(model)
    }
})

test_that("Compute:\tComputation works for the Multinomial method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'Multinomial method' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "mult", maxtime = 10)
        cat("\nAbout to run expect_equal(result$method, 'mult')\n")
        expect_equal(result$method, "mult")
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tComputation works for the MVN CDF method with 'Genz2' method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'MVN CDF Genz2' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "mvn_cdf", maxtime = 10)
        cat("\nAbout to run expect_equal(result$method, 'mvn_cdf')\n")
        expect_equal(result$method, "mvn_cdf")
        cat("\nAbout to run expect_equal(result$mc_method, 'genz2')\n")
        expect_equal(result$mc_method, "genz2")
        cat("\nAbout to run expect_equal(result$mc_error, 1e-6)\n")
        expect_equal(result$mc_error, 1e-6)
        cat("\nAbout to run expect_equal(result$mc_samples, 5000)\n")
        expect_equal(result$mc_samples, 5000)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tComputation works for the MVN CDF method with 'Genz' method", {
    # Test with the "genz" method
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'MVN CDF Genz' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "mvn_cdf", mc_method = "genz", maxtime = 10)
        cat("\nAbout to run expect_equal(result$method, 'mvn_cdf')\n")
        expect_equal(result$method, "mvn_cdf")
        cat("\nAbout to run expect_equal(result$mc_method, 'genz')\n")
        expect_equal(result$mc_method, "genz")
        cat("\nAbout to run expect_equal(result$mc_error, 1e-6)\n")
        expect_equal(result$mc_error, 1e-6)
        cat("\nAbout to run expect_equal(result$mc_samples, 5000)\n")
        expect_equal(result$mc_samples, 5000)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(result$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tComputation works for the MVN CDF method with distinct iterations and sampling methods", {
    # Test with distinct parameters
    for (i in 1:tests) {
        cat(sprintf("Running 'MVN CDF distinct parameters' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "mvn_cdf", mc_method = "genz2", mc_samples = 3000, maxtime = 10)
        cat("\nAbout to run expect_equal(result$mc_samples, 3000)\n")
        expect_equal(result$mc_samples, 3000)
        cat("\nAbout to run expect_equal(result$mc_error, 1e-6)\n")
        expect_equal(result$mc_error, 1e-6)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")

        result <- compute(model, method = "mvn_cdf", mc_method = "genz2", mc_error = 1e-8, maxtime = 10)
        cat("\nAbout to run expect_equal(result$mc_samples, 5000)\n")
        expect_equal(result$mc_samples, 5000)
        cat("\nAbout to run expect_equal(result$mc_error, 1e-8)\n")
        expect_equal(result$mc_error, 1e-8)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")

        result <- compute(model, method = "mvn_cdf", mc_method = "genz", mc_error = 1e-8, mc_samples = 10000, maxtime = 10)
        cat("\nAbout to run expect_equal(result$mc_error, 1e-8)\n")
        expect_equal(result$mc_error, 1e-8)
        cat("\nAbout to run expect_equal(result$mc_samples, 10000)\n")
        expect_equal(result$mc_samples, 10000)
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        rm(model)
        gc()
    }
})

test_that("Compute:\tComputation works for the MVN PDF method", {
    # From the test before, we can assume the constructing function works
    difference <- 0
    for (i in 1:tests) {
        cat(sprintf("Running 'MVN PDF method' test, iteration %d\n", i))
        model <- eim(results[[i]]$X, results[[i]]$W)
        result <- compute(model, method = "mvn_pdf", maxtime = 10)
        cat("\nAbout to run expect_equal(result$method, 'mvn_pdf')\n")
        expect_equal(result$method, "mvn_pdf")
        cat("\nAbout to run expect_type(model$prob, 'double')\n")
        expect_type(model$prob, "double")
        difference <- difference + mean(abs(model$prob - results[[i]]$prob))
        rm(model)
        gc()
    }
    cat("\nAbout to run expect_true(difference / tests < 0.1)\n")
    expect_true(difference / tests < 0.1) # Less than 0.1 of error
})

test_that("Compute:\tVerbose works correctly.", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Verbose' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            cat("\nAbout to run expect_output for verbose computation\n")
            expect_output(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], verbose = TRUE, maxtime = 10
            ))
        }
    }
})

test_that("Compute:\tStop threshold correctly.", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Stop threshold' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            x <- runif(1, min = .Machine$double.eps, max = 1)
            cat("\nAbout to run expect_no_error for stop_threshold\n")
            expect_no_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], stop_threshold = x, maxtime = 10
            ))
        }
    }
})

test_that("Compute:\tErrors when an invalid method is provided", {
    for (i in 1:method_division) {
        cat(sprintf("Running 'Invalid method' test, iteration %d\n", i))
        cat("\nAbout to run expect_error for invalid method 'Hit and Run'\n")
        expect_error(compute(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "Hit and Run", maxtime = 10
        ))
        cat("\nAbout to run expect_error for invalid method 'HNR'\n")
        expect_error(compute(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "HNR", maxtime = 10
        ))
        cat("\nAbout to run expect_error for invalid method 'EXakt'\n")
        expect_error(compute(
            X = results[[i]]$X, W = results[[i]]$W,
            method = "EXakt", maxtime = 10
        ))
    }
})

test_that("Compute:\tErrors when an object and matrices are provided", {
    for (i in 1:method_division) {
        cat(sprintf("Running 'Object and matrices provided' test, iteration %d\n", i))
        cat("\nAbout to run expect_error for providing object and matrices together\n")
        expect_error(compute(
            X = results[[i]]$X, W = results[[i]]$W,
            object = eim(X = results[[i]]$X, W = results[[i]]$W), maxtime = 10
        ))
        dummyobj <- eim(X = results[[i]]$X, W = results[[i]]$W)
        cat("\nAbout to run expect_error for providing object in second call\n")
        expect_error(compute(
            X = results[[i]]$X, W = results[[i]]$W,
            object = dummyobj
        ))
        cat("\nAbout to run expect_error for providing object with only X\n")
        expect_error(compute(X = results[[i]]$X, object = dummyobj))
        cat("\nAbout to run expect_error for providing object with only W\n")
        expect_error(compute(W = results[[i]]$W, object = dummyobj))
    }
})

test_that("Compute:\tErrors when an invalid json is provided", {
    cat("Running 'Invalid json' test\n")
    cat("\nAbout to run expect_error for invalid json path\n")
    expect_error(compute(json_path = "invalid.json"))
})

test_that("Compute:\tErrors when an invalid iteration parameter is provided", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Invalid iteration parameter' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            x <- runif(1, min = .Machine$double.eps, max = 1)
            cat("\nAbout to run expect_error for maxiter = x\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = x, maxtime = 10
            ))
            cat("\nAbout to run expect_error for maxiter = 0\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = 0, maxtime = 10
            ))
            cat("\nAbout to run expect_error for maxiter = -3\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxiter = -3, maxtime = 10
            ))
        }
    }
})

test_that("Compute:\tErrors when an invalid time parameter is provided", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Invalid time parameter' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            x <- runif(1, min = .Machine$double.eps, max = 1)
            cat("\nAbout to run expect_error for maxtime = -3\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], maxtime = -3
            ))
        }
    }
})

test_that("Compute:\tErrors when an invalid stop threshold is provided", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Invalid stop threshold' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            x <- runif(1, min = .Machine$double.eps, max = 1)
            cat("\nAbout to run expect_error for stop_threshold = -3\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], stop_threshold = -3
            ))
        }
    }
})

test_that("Compute:\tErrors when an invalid verbose is provided", {
    for (k in 1:5) {
        for (i in 1:method_division) {
            cat(sprintf(
                "Running 'Invalid verbose' test: method '%s', iteration %d\n",
                available_methods[k], i
            ))
            x <- runif(1, min = .Machine$double.eps, max = 1)
            cat("\nAbout to run expect_error for verbose = 0.1\n")
            expect_error(compute(
                X = results[[i]]$X, W = results[[i]]$W,
                method = available_methods[k], verbose = 0.1
            ))
        }
    }
})
