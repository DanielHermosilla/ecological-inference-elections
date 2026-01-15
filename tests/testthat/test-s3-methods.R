test_that("S3 methods behave as expected", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(20, 4),
        seed = 210
    )

    model <- eim(X = sim$X, W = sim$W)
    expect_error(as.matrix(model), "Probability matrix not available")

    fit <- run_em(
        object = model,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_true(is.matrix(as.matrix(fit)))
    summary_out <- summary(fit)
    expect_true(all(c("candidates", "groups", "ballots") %in% names(summary_out)))
    expect_error(logLik(fit), "run_em", fixed = TRUE)

    fit_ll <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = TRUE
    )
    expect_true(is.numeric(logLik(fit_ll)))
})
