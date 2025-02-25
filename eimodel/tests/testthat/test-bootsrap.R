# ========= Bootstrap function test ========= #
# For avoiding long precomputations, we'll use a small sample => reduce ballots boxes
tests <- 50
results <- vector("list", tests)
for (i in 1:tests) {
    results[[i]] <- eimodel:::.random_samples(
        c(30, 50), # Ballots range
        c(2, 4), # Candidates range
        c(2, 4), # Demographic group range
        c(30, 60), # Votes per ballot range
        seed = i
    )
}

test_that("Bootstrap accepts matrices", {
    for (i in 1:tests) {
        test <- bootstrap(nboot = 10, X = results[[i]]$X, W = results[[i]]$W)
        expect_type(test$sd, "double")
        expect_equal(test$X, results[[i]]$X)
        expect_equal(test$W, results[[i]]$W)
        expect_s3_class(test, "eim")
    }
})
test_that("Bootstrap accepts 'eim' objects", {
    for (i in 1:tests) {
        test <- eim(results[[i]]$X, results[[i]]$W)
        test <- bootstrap(10, test)
        expect_type(test$sd, "double")
        expect_s3_class(test, "eim")
    }
})

test_that("Bootstrap doesn't work when providing multiple sources", {
    for (i in 1:tests) {
        test <- eim(results[[i]]$X, results[[i]]$W)
        expect_error(bootstrap(nboot = 10, test, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = 10, test, W = results[[i]]$W))
        expect_error(bootstrap(nboot = 10, test, X = results[[i]]$X))
    }
})

test_that("Bootstrap doesn't work when providing a single matrix", {
    for (i in 1:tests) {
        expect_error(bootstrap(nboot = 10, X = results[[i]]$X))
    }
})


test_that("Bootstrap doesn't work when providing a invalid json path", {
    for (i in 1:tests) {
        expect_error(bootstrap(json_path = "invalid.json"))
        expect_error(bootstrap(json_path = "invalid"))
        expect_error(bootstrap(json_path = results[[i]]))
    }
})


test_that("Bootstrap doesn't work when providing an invalid iteration parameter", {
    for (i in 1:tests) {
        expect_error(bootstrap(nboot = 1.1, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = 0.21353, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = 1e-9, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = -2, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = -2.1, X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = "a", X = results[[i]]$X, W = results[[i]]$W))
        expect_error(bootstrap(nboot = NULL, X = results[[i]]$X, W = results[[i]]$W))
    }
})

test_that("Standard deviation matrix has dimensional coherence", {
    for (i in 1:tests) {
        test <- eim(results[[i]]$X, results[[i]]$W)
        test <- bootstrap(10, test)
        # Have the same group dimension
        expect_equal(nrow(test$sd), ncol(test$W))
        # Have the same candidate dimension
        expect_equal(ncol(test$sd), ncol(test$X))
    }
})

test_that("Bootstrap computes the probability matrix if required", {
    half <- tests %/% 2
    for (i in 1:half) {
        test <- eim(results[[i]]$X, results[[i]]$W)
        test <- bootstrap(10, test, get_p = TRUE)
        expect_type(test$prob, "double")
        expect_equal(test$method, "mult")
        expect_gt(test$iterations, 0)
        expect_type(logLik(test), "double")
    }
    for (i in 1:half) {
        test <- eim(results[[i]]$X, results[[i]]$W)
        test <- bootstrap(nboot = 30, object = test, get_p = FALSE)
        expect_null(test$prob)
        expect_null(test$method)
        expect_null(test$iterations)
        expect_null(test$logLik)
    }
})
