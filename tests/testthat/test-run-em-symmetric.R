test_that("run_em symmetric uses joint by default", {
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
        maxiter = 6,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_equal(fit$symmetric_weight_method, "joint")
    expect_equal(fit$symmetric_weights, c(original = 0.5, reverse = 0.5))
    expect_null(fit$prob_inv)
    expect_null(fit$cond_prob_inv)
    expect_null(fit$expected_outcome_inv)
    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
})

test_that("run_em symmetric average preserves inverse outputs", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 146
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "average",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_equal(fit$symmetric_weight_method, "average")
    expect_equal(fit$symmetric_weights, c(original = 0.5, reverse = 0.5))
    expect_true(is.matrix(fit$prob_inv))
    expect_true(is.array(fit$cond_prob_inv))
    expect_equal(dim(fit$prob_inv), c(ncol(sim$X), ncol(sim$W)))
    expect_equal(dim(fit$cond_prob_inv), c(ncol(sim$X), ncol(sim$W), nrow(sim$X)))
    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)

    sums_inv <- apply(fit$cond_prob_inv, c(1, 3), sum)
    expect_true(all(abs(sums_inv - 1) < 1e-6))
})

test_that("run_em symmetric supports delta_ll weighting", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 142
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "delta_ll",
        maxiter = 4,
        maxtime = 2,
        compute_ll = TRUE
    )

    expect_true(is.matrix(fit$prob_inv))
    expect_true(is.array(fit$cond_prob_inv))
    expect_equal(fit$symmetric_weight_method, "delta_ll")
    expect_true(is.numeric(fit$symmetric_weights))
    expect_equal(names(fit$symmetric_weights), c("original", "reverse"))
    expect_equal(sum(fit$symmetric_weights), 1, tolerance = 1e-8)

    expect_true(is.numeric(fit$LL_ind))
    expect_true(is.numeric(fit$LL_rev_ind))
    expect_true(is.numeric(fit$dLL))
    expect_true(is.numeric(fit$dLL_rev))
    expect_true(is.numeric(fit$nu))
    expect_true(is.numeric(fit$nu_rev))
})

test_that("run_em symmetric supports mae_inverse weighting", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 143
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "mae_inverse",
        maxiter = 4,
        maxtime = 2,
        compute_ll = TRUE
    )

    expect_equal(fit$symmetric_weight_method, "mae_inverse")
    expect_true(is.numeric(fit$symmetric_weights))
    expect_equal(names(fit$symmetric_weights), c("original", "reverse"))
    expect_equal(sum(fit$symmetric_weights), 1, tolerance = 1e-8)
    expect_true(is.numeric(fit$err_forward))
    expect_true(is.numeric(fit$err_inverse))
})

test_that("run_em symmetric supports joint", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 144
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "joint",
        maxiter = 6,
        maxtime = 2,
        compute_ll = TRUE
    )

    expect_equal(fit$symmetric_weight_method, "joint")
    expect_equal(fit$symmetric_weights, c(original = 0.5, reverse = 0.5))
    expect_null(fit$prob_inv)
    expect_null(fit$cond_prob_inv)
    expect_null(fit$expected_outcome_inv)
    expect_false(isTRUE(fit$adjust_prob_cond_every))
    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
})

test_that("run_em joint supports lp with adjust_prob_cond_every TRUE/FALSE", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 145
    )

    fit_false <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "joint",
        adjust_prob_cond_method = "lp",
        adjust_prob_cond_every = FALSE,
        maxiter = 6,
        maxtime = 2,
        compute_ll = FALSE
    )

    fit_true <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        symmetric = TRUE,
        symmetric_weight_method = "joint",
        adjust_prob_cond_method = "lp",
        adjust_prob_cond_every = TRUE,
        maxiter = 6,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_false(isTRUE(fit_false$adjust_prob_cond_every))
    expect_true(isTRUE(fit_true$adjust_prob_cond_every))
    expect_prob_matrix(fit_false$prob)
    expect_prob_matrix(fit_true$prob)
    expect_prob_array(fit_false$cond_prob)
    expect_prob_array(fit_true$cond_prob)
})
