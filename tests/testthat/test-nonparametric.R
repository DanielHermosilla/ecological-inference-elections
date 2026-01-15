test_that("simulate_election returns non-parametric data", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        seed = 123
    )

    expect_s3_class(sim, "eim")
    expect_equal(dim(sim$X), c(10, 3))
    expect_equal(dim(sim$W), c(10, 2))
    expect_true(is.matrix(sim$real_prob))
    expect_equal(dim(sim$real_prob), c(2, 3))
})

test_that("run_em returns non-parametric probability matrix", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        seed = 124
    )

    model <- eim(X = sim$X, W = sim$W)
    fit <- run_em(
        object = model,
        method = "mult",
        maxiter = 5,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(fit, "eim")
    expect_true(is.matrix(fit$prob))
    expect_equal(dim(fit$prob), c(2, 3))
    expect_true(is.array(fit$cond_prob))
    expect_equal(dim(fit$cond_prob), c(2, 3, 10))
})

test_that("bootstrap returns non-parametric standard deviations", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        seed = 125
    )

    model <- eim(X = sim$X, W = sim$W)
    boot <- bootstrap(
        object = model,
        nboot = 3,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(boot, "eim")
    expect_true(is.matrix(boot$sd))
    expect_equal(dim(boot$sd), c(2, 3))
})
