test_that("waldtest returns p-values matrix", {
    sim1 <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 190
    )
    sim2 <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 191
    )

    result <- waldtest(
        X1 = sim1$X,
        W1 = sim1$W,
        X2 = sim2$X,
        W2 = sim2$W,
        nboot = 2,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_true(is.matrix(result$pvals))
    expect_equal(dim(result$pvals), c(2, 3))
    expect_s3_class(result$eim1, "eim")
    expect_s3_class(result$eim2, "eim")
})
