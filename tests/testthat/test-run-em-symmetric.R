test_that("run_em symmetric mode returns inverse probabilities", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 140
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_true(is.matrix(fit$prob_inv))
    expect_true(is.array(fit$cond_prob_inv))
    expect_equal(dim(fit$prob_inv), c(ncol(sim$X), ncol(sim$W)))
    expect_equal(dim(fit$cond_prob_inv), c(ncol(sim$X), ncol(sim$W), nrow(sim$X)))

    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)

    sums_inv <- apply(fit$cond_prob_inv, c(1, 3), sum)
    expect_true(all(abs(sums_inv - 1) < 1e-6))
})

test_that("run_em symmetric works in parametric mode", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 141
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        symmetric = TRUE,
        maxiter = 4,
        maxtime = 2
    )

    expect_true(is.array(fit$prob))
    expect_true(is.array(fit$cond_prob))
    expect_true(is.array(fit$prob_inv))
    expect_true(is.array(fit$cond_prob_inv))

    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$prob_inv), c(ncol(sim$X), ncol(sim$W), nrow(sim$X)))
    expect_equal(dim(fit$cond_prob_inv), c(ncol(sim$X), ncol(sim$W), nrow(sim$X)))

    expect_prob_array(fit$prob)
    expect_prob_array(fit$cond_prob)

    sums_inv <- apply(fit$cond_prob_inv, c(1, 3), sum)
    expect_true(all(abs(sums_inv - 1) < 1e-6))
})
