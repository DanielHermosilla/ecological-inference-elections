test_that("parametric bootstrap returns beta and alpha deviations", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        parametric = TRUE,
        num_attributes = 2,
        num_districts = 2,
        seed = 160
    )

    boot <- bootstrap(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        nboot = 3,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    expect_s3_class(boot, "eim")
    expect_equal(dim(boot$sd_beta), c(2, 2))
    expect_equal(dim(boot$sd_alpha), c(2, 2))
    expect_true(all(boot$sd_beta >= 0))
    expect_true(all(boot$sd_alpha >= 0))
    expect_equal(boot$nboot, 3)
})
