test_that("run_em aggregates groups when group_agg is provided", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 4,
        ballot_voters = rep(30, 6),
        seed = 180
    )

    group_agg <- c(2, 4)
    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        group_agg = group_agg,
        maxiter = 4,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_equal(fit$group_agg, group_agg)
    expect_equal(dim(fit$W_agg), c(6, 2))
    expect_equal(dim(fit$prob), c(2, 3))
})

test_that("get_agg_proxy returns aggregated results for non-parametric models", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 3,
        ballot_voters = rep(30, 6),
        seed = 181
    )

    proxy <- get_agg_proxy(
        X = sim$X,
        W = sim$W,
        nboot = 2,
        sd_threshold = 0.2,
        feasible = FALSE,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(proxy, "eim")
    expect_true(!is.null(proxy$W_agg))
    expect_true(!is.null(proxy$group_agg))
    expect_true(is.matrix(proxy$sd))
})

test_that("get_agg_opt handles aggregation outputs", {
    sim <- simulate_election(
        num_ballots = 5,
        num_candidates = 3,
        num_groups = 3,
        ballot_voters = rep(30, 5),
        seed = 182
    )

    opt <- get_agg_opt(
        X = sim$X,
        W = sim$W,
        nboot = 2,
        sd_threshold = 0.2,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    expect_s3_class(opt, "eim")
    if (!is.null(opt$group_agg)) {
        expect_true(!is.null(opt$W_agg))
        expect_true(is.matrix(opt$sd))
        expect_equal(ncol(opt$W_agg), length(opt$group_agg))
    }
})
