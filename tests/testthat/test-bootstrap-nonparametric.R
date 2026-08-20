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
    expect_false(boot$bootstrap_symmetric)
    expect_null(boot$bootstrap_symmetric_weight_method)
})

test_that("non-parametric bootstrap runs the joint symmetric EM in each replicate", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 8),
        seed = 151
    )

    args <- list(
        X = sim$X,
        W = sim$W,
        nboot = 4,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE,
        adjust_prob_cond_method = "",
        seed = 901
    )

    regular <- do.call(bootstrap, args)
    symmetric <- do.call(bootstrap, c(args, list(
        symmetric = TRUE,
        symmetric_weight_method = "joint"
    )))

    expect_true(symmetric$bootstrap_symmetric)
    expect_equal(symmetric$bootstrap_symmetric_weight_method, "joint")
    expect_prob_matrix(symmetric$avg_prob)
    expect_false(isTRUE(all.equal(symmetric$avg_prob, regular$avg_prob, tolerance = 1e-12)))
})

test_that("symmetric bootstrap rejects unsupported weighting methods", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = rep(20, 4),
        seed = 152
    )

    expect_error(
        bootstrap(
            X = sim$X,
            W = sim$W,
            nboot = 2,
            symmetric = TRUE,
            symmetric_weight_method = "average"
        ),
        "supports only symmetric_weight_method = \"joint\"",
        fixed = TRUE
    )
})
