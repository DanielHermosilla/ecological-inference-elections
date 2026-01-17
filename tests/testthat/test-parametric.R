test_that("simulate_election returns parametric data", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 50,
        num_covariates = 2,
        num_districts = 2,
        seed = 200
    )

    expect_s3_class(sim, "eim")
    expect_true(is.matrix(sim$V))
    expect_equal(dim(sim$V), c(8, 2))
    expect_true(is.matrix(sim$real_alpha))
    expect_true(is.matrix(sim$real_beta))
    expect_equal(dim(sim$real_alpha), c(2, 2))
    expect_equal(dim(sim$real_beta), c(2, 2))
})

test_that("run_em works in parametric mode", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 201
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    fit <- run_em(
        object = model,
        method = "mult",
        beta_init = sim$real_beta,
        alpha_init = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    expect_s3_class(fit, "eim")
    expect_true(is.array(fit$prob))
    expect_equal(dim(fit$prob), c(2, 3, 8))
    expect_equal(dim(fit$beta), c(2, 2))
    expect_equal(dim(fit$alpha), c(2, 2))
})

test_that("bootstrap works in parametric mode", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 202
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    boot <- bootstrap(
        object = model,
        nboot = 3,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    expect_s3_class(boot, "eim")
    expect_true(is.matrix(boot$sd_beta))
    expect_true(is.matrix(boot$sd_alpha))
    expect_equal(dim(boot$sd_beta), c(2, 2))
    expect_equal(dim(boot$sd_alpha), c(2, 2))
})

test_that("parametric run_em rejects non-mult methods", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        num_covariates = 2,
        num_districts = 2,
        seed = 203
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    expect_error(
        run_em(object = model, method = "mvn_cdf"),
        "Parametric mode only supports method"
    )
})
