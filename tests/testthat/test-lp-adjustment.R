test_that("LP adjustment matches candidate totals", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(60, 8),
        lambda = 0.3,
        seed = 130
    )

    fit_lp <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE,
        adjust_prob_cond_method = "lp",
        adjust_prob_cond_every = TRUE
    )

    fit_project <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE,
        adjust_prob_cond_method = "project_lp",
        adjust_prob_cond_every = TRUE
    )

    xhat_lp <- expected_votes_from_q(sim$W, fit_lp$cond_prob)
    xhat_project <- expected_votes_from_q(sim$W, fit_project$cond_prob)

    expect_equal(xhat_lp, sim$X, tolerance = 1e-4)
    expect_equal(xhat_project, sim$X, tolerance = 1e-4)
}) 
