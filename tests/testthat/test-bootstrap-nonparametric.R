test_that("bootstrap returns standard deviation and average probability", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 8),
        seed = 150
    )

    boot <- bootstrap(
        X = sim$X,
        W = sim$W,
        nboot = 3,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(boot, "eim")
    expect_equal(dim(boot$sd), c(2, 3))
    expect_equal(dim(boot$avg_prob), c(2, 3))
    expect_true(all(boot$sd >= 0))
    expect_prob_matrix(boot$avg_prob)
})
