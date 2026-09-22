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

    fit_kl <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE,
        adjust_prob_cond_method = "kl",
        adjust_prob_cond_every = TRUE
    )

    xhat_lp <- expected_votes_from_q(sim$W, fit_lp$cond_prob)
    xhat_project <- expected_votes_from_q(sim$W, fit_project$cond_prob)
    xhat_kl <- expected_votes_from_q(sim$W, fit_kl$cond_prob)

    expect_equal(xhat_lp, sim$X, tolerance = 1e-4)
    expect_equal(xhat_project, sim$X, tolerance = 1e-4)
    expect_equal(xhat_kl, sim$X, tolerance = 1e-7)
}) 

test_that("symmetric LP adjustment is stable with large ballot counts", {
    set.seed(1)
    B <- 20L
    G <- 4L
    C <- 4L
    ballot_voters <- pmax(1000L, round(exp(rnorm(B, log(17000), 1))))
    W <- t(vapply(
        ballot_voters,
        function(n) as.vector(rmultinom(1, n, rep(1 / G, G))),
        numeric(G)
    ))

    transition <- matrix(0.1, nrow = G, ncol = C)
    diag(transition) <- 0.7
    X <- matrix(0, nrow = B, ncol = C)
    for (b in seq_len(B)) {
        for (g in seq_len(G)) {
            X[b, ] <- X[b, ] + as.vector(rmultinom(1, W[b, g], transition[g, ]))
        }
    }

    fit <- run_em(
        X = X,
        W = W,
        method = "mult",
        symmetric = TRUE,
        adjust_prob_cond_method = "lp",
        maxiter = 3,
        initial_prob = "group_proportional",
        compute_ll = FALSE,
        verbose = FALSE
    )

    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
    expect_lt(max(abs(expected_votes_from_q(W, fit$cond_prob) - X)), 1e-6)
})
