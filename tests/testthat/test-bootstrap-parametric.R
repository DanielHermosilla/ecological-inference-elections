test_that("parametric bootstrap returns beta and alpha deviations", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
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
    expect_false(boot$bootstrap_symmetric)
    expect_null(boot$bootstrap_symmetric_weight_method)
})

test_that("parametric bootstrap runs the joint symmetric EM in each replicate", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 8),
        num_covariates = 2,
        num_districts = 2,
        seed = 161
    )

    args <- list(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        nboot = 4,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        maxnewton = 1,
        adjust_prob_cond_method = "",
        seed = 900
    )

    regular <- do.call(bootstrap, args)
    symmetric <- do.call(bootstrap, c(args, list(
        symmetric = TRUE,
        symmetric_weight_method = "joint"
    )))

    expect_true(symmetric$bootstrap_symmetric)
    expect_equal(symmetric$bootstrap_symmetric_weight_method, "joint")
    expect_equal(dim(symmetric$sd_beta), dim(regular$sd_beta))
    expect_equal(dim(symmetric$sd_alpha), dim(regular$sd_alpha))
    expect_null(symmetric$sd_beta_inv)
    expect_null(symmetric$sd_alpha_inv)
    expect_false(isTRUE(all.equal(
        c(symmetric$sd_beta, symmetric$sd_alpha),
        c(regular$sd_beta, regular$sd_alpha),
        tolerance = 1e-12
    )))
})

test_that("non-joint parametric bootstrap returns deviations for both directions", {
    sim <- simulate_election(
        num_ballots = 7,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 7),
        num_covariates = 2,
        num_districts = 2,
        seed = 163
    )

    for (weight_method in c("average", "delta_ll", "mae_inverse")) {
        boot <- bootstrap(
            X = sim$X,
            W = sim$W,
            V = sim$V,
            nboot = 3,
            method = "mult",
            maxiter = 2,
            maxtime = 2,
            maxnewton = 1,
            adjust_prob_cond_method = "",
            symmetric = TRUE,
            symmetric_weight_method = weight_method,
            seed = 903
        )

        expect_true(boot$bootstrap_symmetric, info = weight_method)
        expect_equal(boot$bootstrap_symmetric_weight_method, weight_method)
        expect_equal(dim(boot$sd_beta), c(2, 2), info = weight_method)
        expect_equal(dim(boot$sd_alpha), c(2, 2), info = weight_method)
        expect_equal(dim(boot$sd_beta_inv), c(3, 1), info = weight_method)
        expect_equal(dim(boot$sd_alpha_inv), c(1, 2), info = weight_method)
        expect_true(all(is.finite(boot$sd_beta_inv)), info = weight_method)
        expect_true(all(is.finite(boot$sd_alpha_inv)), info = weight_method)
    }
})

test_that("parametric symmetric bootstrap supports the exact E-step", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = rep(8, 4),
        num_covariates = 1,
        num_districts = 2,
        seed = 162
    )

    boot <- bootstrap(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        nboot = 2,
        method = "exact",
        maxiter = 1,
        maxtime = 5,
        maxnewton = 1,
        adjust_prob_cond_method = "",
        symmetric = TRUE,
        symmetric_weight_method = "joint",
        seed = 902
    )

    expect_true(boot$bootstrap_symmetric)
    expect_equal(boot$bootstrap_symmetric_weight_method, "joint")
    expect_true(all(is.finite(boot$sd_beta)))
    expect_true(all(is.finite(boot$sd_alpha)))
})
