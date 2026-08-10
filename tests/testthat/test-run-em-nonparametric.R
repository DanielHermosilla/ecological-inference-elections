test_that("run_em returns consistent non-parametric outputs", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(50, 10),
        lambda = 0.3,
        seed = 110
    )

    colnames(sim$X) <- paste0("C", seq_len(ncol(sim$X)))
    colnames(sim$W) <- paste0("G", seq_len(ncol(sim$W)))
    rownames(sim$X) <- paste0("B", seq_len(nrow(sim$X)))
    rownames(sim$W) <- rownames(sim$X)

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 5,
        miniter = 1,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(fit, "eim")
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$expected_outcome), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))

    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)

    expect_dimnames_match(fit$cond_prob, colnames(sim$W), colnames(sim$X), rownames(sim$X))
    expect_equal(dimnames(fit$prob)[[1]], colnames(sim$W))
    expect_equal(dimnames(fit$prob)[[2]], colnames(sim$X))

    expected_by_group <- apply(fit$expected_outcome, c(1, 3), sum)
    expect_equal(expected_by_group, t(sim$W), tolerance = 1e-6)
})

test_that("run_em returns finite outputs for ballot boxes with zero totals", {
    base_X <- matrix(c(4, 3, 3, 3, 2, 5, 2, 5, 3), nrow = 3, byrow = TRUE)
    base_W <- matrix(c(6, 4, 7, 3, 5, 5), nrow = 3, byrow = TRUE)

    scenarios <- list(
        both_zero = list(X = base_X, W = base_W),
        W_zero = list(X = base_X, W = base_W),
        X_zero = list(X = base_X, W = base_W)
    )
    scenarios$both_zero$X[2, ] <- 0
    scenarios$both_zero$W[2, ] <- 0
    scenarios$W_zero$W[2, ] <- 0
    scenarios$X_zero$X[2, ] <- 0

    for (method in c("mult", "mvn_pdf")) {
        for (scenario_name in names(scenarios)) {
            scenario <- scenarios[[scenario_name]]
            fit <- run_em(
                X = scenario$X,
                W = scenario$W,
                method = method,
                maxiter = 2,
                miniter = 1,
                maxtime = 2,
                compute_ll = FALSE
            )

            info <- paste("method:", method, "scenario:", scenario_name)
            expect_true(all(is.finite(fit$cond_prob)), info = info)
            expect_true(all(is.finite(fit$expected_outcome)), info = info)
            expect_true(all(is.finite(fit$prob)), info = info)
            expect_equal(
                fit$cond_prob[, , 2],
                matrix(0, nrow = ncol(scenario$W), ncol = ncol(scenario$X)),
                info = info
            )
            expect_equal(
                fit$expected_outcome[, , 2],
                matrix(0, nrow = ncol(scenario$W), ncol = ncol(scenario$X)),
                info = info
            )
            expect_prob_matrix(fit$prob)
        }
    }
})

test_that("run_em approximates simulated probabilities", {
    sim <- simulate_election(
        num_ballots = 20,
        num_candidates = 3,
        num_groups = 3,
        ballot_voters = rep(80, 20),
        lambda = 0.2,
        seed = 111
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 30,
        miniter = 5,
        maxtime = 4,
        compute_ll = FALSE
    )

    error_mae <- mae(sim$real_prob, fit$prob)
    expect_lt(error_mae, 0.4)
})
