test_that("run_em restores zero-vote candidate columns in standard outputs", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(50, 10),
        lambda = 0.3,
        seed = 310
    )

    rownames(sim$X) <- paste0("B", seq_len(nrow(sim$X)))
    rownames(sim$W) <- rownames(sim$X)
    colnames(sim$W) <- paste0("G", seq_len(ncol(sim$W)))

    X_zero <- cbind(
        C1 = sim$X[, 1],
        C0 = 0,
        C2 = sim$X[, 2],
        C3 = sim$X[, 3]
    )
    rownames(X_zero) <- rownames(sim$X)

    fit <- run_em(
        X = X_zero,
        W = sim$W,
        method = "mult",
        maxiter = 5,
        miniter = 1,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_equal(fit$X, X_zero)
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(X_zero)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(X_zero), nrow(X_zero)))
    expect_equal(dim(fit$expected_outcome), c(ncol(sim$W), ncol(X_zero), nrow(X_zero)))

    expect_equal(unname(fit$prob[, "C0"]), rep(0, ncol(sim$W)))
    expect_equal(unname(fit$cond_prob[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))
    expect_equal(unname(fit$expected_outcome[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))

    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
    expect_equal(dimnames(fit$prob)[[2]], colnames(X_zero))
    expect_dimnames_match(fit$cond_prob, colnames(sim$W), colnames(X_zero), rownames(X_zero))

    expected_by_group <- apply(fit$expected_outcome, c(1, 3), sum)
    expect_equal(expected_by_group, t(sim$W), tolerance = 1e-6)
})

test_that("run_em restores zero-vote candidate columns in mixture and symmetric outputs", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 6),
        lambda = 0.25,
        seed = 311
    )

    rownames(sim$X) <- paste0("B", seq_len(nrow(sim$X)))
    rownames(sim$W) <- rownames(sim$X)
    colnames(sim$W) <- paste0("G", seq_len(ncol(sim$W)))

    X_zero <- cbind(
        C1 = sim$X[, 1],
        C0 = 0,
        C2 = sim$X[, 2],
        C3 = sim$X[, 3]
    )
    rownames(X_zero) <- rownames(sim$X)

    fit <- run_em(
        X = X_zero,
        W = sim$W,
        method = "mult",
        mixture = 2,
        symmetric = TRUE,
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_equal(dim(fit$component_prob), c(ncol(sim$W), ncol(X_zero), 2))
    expect_equal(dim(fit$prob_inv), c(ncol(X_zero), ncol(sim$W)))
    expect_equal(dim(fit$cond_prob_inv), c(ncol(X_zero), ncol(sim$W), nrow(X_zero)))

    expect_equal(unname(fit$component_prob[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = 2))
    expect_equal(unname(fit$prob_inv["C0", ]), rep(0, ncol(sim$W)))
    expect_equal(unname(fit$cond_prob_inv["C0", , ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))

    expect_equal(dimnames(fit$prob_inv)[[1]], colnames(X_zero))
    expect_equal(dimnames(fit$cond_prob_inv)[[1]], colnames(X_zero))
})

test_that("run_em restores zero-vote candidate columns in parametric outputs", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 313
    )

    rownames(sim$X) <- paste0("B", seq_len(nrow(sim$X)))
    rownames(sim$W) <- rownames(sim$X)
    colnames(sim$W) <- paste0("G", seq_len(ncol(sim$W)))
    colnames(sim$V) <- paste0("A", seq_len(ncol(sim$V)))

    X_zero <- cbind(
        C1 = sim$X[, 1],
        C0 = 0,
        C2 = sim$X[, 2],
        C3 = sim$X[, 3]
    )
    rownames(X_zero) <- rownames(sim$X)

    fit <- run_em(
        X = X_zero,
        W = sim$W,
        V = sim$V,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        maxnewton = 1
    )

    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(X_zero), nrow(X_zero)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(X_zero), nrow(X_zero)))
    expect_equal(dim(fit$expected_outcome), c(ncol(sim$W), ncol(X_zero), nrow(X_zero)))
    expect_equal(dim(fit$beta), c(ncol(sim$W), ncol(X_zero) - 1))
    expect_equal(dim(fit$alpha), c(ncol(X_zero) - 1, ncol(sim$V)))

    expect_equal(unname(fit$prob[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))
    expect_equal(unname(fit$cond_prob[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))
    expect_equal(unname(fit$expected_outcome[, "C0", ]), matrix(0, nrow = ncol(sim$W), ncol = nrow(X_zero)))

    expect_equal(colnames(fit$beta), colnames(X_zero)[-ncol(X_zero)])
    expect_equal(rownames(fit$alpha), colnames(X_zero)[-ncol(X_zero)])
})
