test_that("run_em returns parametric probabilities and parameters", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        parametric = TRUE,
        num_covariates = 2,
        num_districts = 2,
        seed = 170
    )

    colnames(sim$X) <- paste0("C", seq_len(ncol(sim$X)))
    colnames(sim$W) <- paste0("G", seq_len(ncol(sim$W)))
    rownames(sim$X) <- paste0("B", seq_len(nrow(sim$X)))
    rownames(sim$W) <- rownames(sim$X)
    colnames(sim$V) <- paste0("A", seq_len(ncol(sim$V)))

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 8,
        maxtime = 3,
        maxnewton = 1
    )

    expect_s3_class(fit, "eim")
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$expected_outcome), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$beta), c(ncol(sim$W), ncol(sim$X) - 1))
    expect_equal(dim(fit$alpha), c(ncol(sim$X) - 1, ncol(sim$V)))
    expect_equal(fit$maxnewton, 1)

    expect_prob_array(fit$prob)
    expect_prob_array(fit$cond_prob)

    expect_equal(dimnames(fit$prob)[[1]], colnames(sim$W))
    expect_equal(dimnames(fit$prob)[[2]], colnames(sim$X))
    expect_equal(dimnames(fit$prob)[[3]], rownames(sim$X))
})

test_that("parametric run_em approximates simulated probabilities", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 50,
        parametric = TRUE,
        num_covariates = 2,
        num_districts = 2,
        seed = 171
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 10,
        maxtime = 4,
        maxnewton = 1
    )

    G <- ncol(sim$W)
    C <- ncol(sim$X)
    B <- nrow(sim$X)
    true_prob <- sim$real_prob
    error_mae <- mae(true_prob, fit$prob)
    expect_lt(error_mae, 0.5)
})
