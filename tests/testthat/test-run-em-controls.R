test_that("run_em respects iteration controls", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        lambda = 0.4,
        seed = 120
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 3,
        miniter = 2,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_true(fit$iterations >= 2)
    expect_true(fit$iterations <= 3)
})

test_that("run_em rejects mismatched totals when disabled", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(25, 4),
        seed = 121
    )

    sim$X[1, 1] <- sim$X[1, 1] + 1

    expect_error(
        run_em(X = sim$X, W = sim$W, allow_mismatch = FALSE),
        "Row-wise mismatch"
    )
})
