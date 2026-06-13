test_that("run_em supports finite mixtures with mixture > 1", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(40, 8),
        lambda = 0.3,
        seed = 901
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        mixture = 2,
        maxiter = 4,
        miniter = 1,
        maxtime = 5,
        compute_ll = TRUE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_equal(length(fit$phi), 2)
    expect_equal(sum(fit$phi), 1, tolerance = 1e-6)
    expect_equal(dim(fit$responsibilities), c(nrow(sim$X), 2))
    expect_equal(dim(fit$component_prob), c(ncol(sim$W), ncol(sim$X), 2))
    expect_equal(rowSums(fit$responsibilities), rep(1, nrow(sim$X)), tolerance = 1e-6)

    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
    expect_prob_array(fit$component_prob)
})

test_that("run_em accepts component-specific initial probabilities for finite mixtures", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = rep(40, 8),
        lambda = 0.3,
        seed = 911
    )

    initial_prob <- array(0, dim = c(2, 2, 2))
    initial_prob[, , 1] <- matrix(c(
        0.8, 0.2,
        0.7, 0.3
    ), nrow = 2, byrow = TRUE)
    initial_prob[, , 2] <- matrix(c(
        0.2, 0.8,
        0.3, 0.7
    ), nrow = 2, byrow = TRUE)

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        mixture = 2,
        initial_prob = initial_prob,
        maxiter = 3,
        miniter = 1,
        maxtime = 5,
        compute_ll = TRUE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_equal(dim(fit$component_prob), c(ncol(sim$W), ncol(sim$X), 2))
    expect_equal(dim(fit$responsibilities), c(nrow(sim$X), 2))
    expect_prob_array(fit$component_prob)
})

test_that("run_em supports saddlepoint method with finite mixtures", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 904
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "saddlepoint",
        mixture = 2,
        maxiter = 3,
        miniter = 1,
        maxtime = 5,
        compute_ll = TRUE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_equal(length(fit$phi), 2)
    expect_equal(sum(fit$phi), 1, tolerance = 1e-6)
    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
})

test_that("run_em validates mixture argument", {
    sim <- simulate_election(
        num_ballots = 5,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 5),
        seed = 902
    )

    expect_error(
        run_em(X = sim$X, W = sim$W, mixture = 0),
        "Invalid 'mixture'"
    )
    expect_error(
        run_em(X = sim$X, W = sim$W, mixture = 1.5),
        "Invalid 'mixture'"
    )
})

test_that("run_em allows row_mixture > G without truncation", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 911
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        row_mixture = 5,
        maxiter = 3,
        miniter = 1,
        maxtime = 5,
        compute_ll = FALSE
    )

    expect_equal(fit$row_mixture, 5L)
    expect_equal(fit$K, 5L)
    expect_equal(fit$mixture, 5L)
    expect_equal(dim(fit$phi), c(ncol(sim$W), 5))
})

test_that("run_em prioritizes row-level backend when row_mixture > 1 and mixture > 1", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 912
    )

    expect_warning(
        fit <- run_em(
            X = sim$X,
            W = sim$W,
            method = "mult",
            mixture = 3,
            row_mixture = 2,
            maxiter = 3,
            miniter = 1,
            maxtime = 5,
            compute_ll = FALSE
        ),
        "implies `K = row_mixture`"
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$row_mixture, 2L)
    expect_equal(fit$K, 2L)
    expect_equal(fit$mixture, 2L)
    expect_equal(dim(fit$phi), c(ncol(sim$W), 2))
    expect_equal(ncol(fit$responsibilities), 2^ncol(sim$W))
})

test_that("run_em supports symmetric finite mixtures", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(35, 6),
        seed = 903
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        mixture = 2,
        symmetric = TRUE,
        maxiter = 4,
        miniter = 1,
        maxtime = 5,
        compute_ll = FALSE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_false(is.null(fit$cond_prob_inv))
    expect_false(is.null(fit$prob_inv))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X)))
    expect_prob_matrix(fit$prob)
    expect_prob_array(fit$cond_prob)
})

test_that("run_em supports symmetric finite mixtures with EM_weight joint EM", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(35, 6),
        seed = 904
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        mixture = 2,
        symmetric = TRUE,
        symmetric_weight_method = "EM_weight",
        maxiter = 4,
        miniter = 1,
        maxtime = 5,
        compute_ll = FALSE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_equal(fit$symmetric_weight_method, "EM_weight")
    expect_equal(fit$symmetric_weights, c(original = 0.5, reverse = 0.5))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$cond_prob_inv), c(ncol(sim$X), ncol(sim$W), nrow(sim$X)))
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X)))
    expect_equal(dim(fit$prob_inv), c(ncol(sim$X), ncol(sim$W)))
    expect_prob_matrix(fit$prob)
    expect_prob_matrix(fit$prob_inv)
    expect_prob_array(fit$cond_prob)
    expect_prob_array(fit$cond_prob_inv)
    expect_true(all(abs(rowSums(fit$responsibilities) - 1) < 1e-6))
})

test_that("run_em symmetric finite mixtures ignore log-likelihood threshold", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(35, 6),
        seed = 904
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        mixture = 2,
        symmetric = TRUE,
        symmetric_weight_method = "EM_weight",
        maxiter = 6,
        miniter = 1,
        maxtime = 5,
        param_threshold = 0,
        ll_threshold = Inf,
        compute_ll = TRUE
    )

    expect_equal(fit$status, 2L)
    expect_equal(fit$iterations, 6)
    expect_true(is.numeric(fit$logLik))
    expect_equal(length(fit$logLik), 1)
})

test_that("run_em rejects EM_weight when row_mixture > 1", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(35, 6),
        seed = 905
    )

    expect_error(
        run_em(
            X = sim$X,
            W = sim$W,
            method = "mult",
            row_mixture = 2,
            symmetric = TRUE,
            symmetric_weight_method = "EM_weight",
            maxiter = 4,
            miniter = 1,
            maxtime = 5
        ),
        "row_mixture = 1"
    )
})

test_that("finite mixture multistart generates same-sign and crossed initializations", {
    X <- matrix(c(
        30, 20,
        25, 25,
        20, 30
    ), nrow = 3, byrow = TRUE)
    W <- matrix(c(
        20, 30,
        25, 25,
        30, 20
    ), nrow = 3, byrow = TRUE)

    samples <- fastei:::.run_em_mixture_initial_prob_samples(
        X_matrix = X,
        W_matrix = W,
        K = 2,
        S = 2,
        amplitude = 1,
        noise_sd = 0
    )

    expect_equal(dim(samples), c(2L, 2L, 2L, 2L))
    for (s in seq_len(2)) {
        for (k in seq_len(2)) {
            expect_equal(rowSums(samples[, , k, s]), rep(1, 2), tolerance = 1e-8)
        }
    }

    profile_delta <- samples[, 1, 1, ] - samples[, 1, 2, ]
    expect_true(any(profile_delta[1, ] > 0 & profile_delta[2, ] > 0))
    expect_true(any(profile_delta[1, ] > 0 & profile_delta[2, ] < 0))
})

test_that("EMLogLikFromProbMixture evaluates a fixed component probability array", {
    X <- matrix(c(
        30, 20,
        25, 25,
        20, 30
    ), nrow = 3, byrow = TRUE)
    W <- matrix(c(
        20, 30,
        25, 25,
        30, 20
    ), nrow = 3, byrow = TRUE)
    initial_prob <- fastei:::.run_em_mixture_initial_prob_samples(
        X_matrix = X,
        W_matrix = W,
        K = 2,
        S = 1,
        amplitude = 1,
        noise_sd = 0
    )
    initial_prob <- initial_prob[, , , 1]

    ll <- EMLogLikFromProbMixture(
        t(X),
        W,
        initial_prob,
        "mult",
        as.integer(3000),
        as.integer(1000),
        "genz",
        as.numeric(1e-3),
        as.integer(5000),
        as.integer(0),
        "project_lp",
        FALSE
    )

    expect_type(ll, "double")
    expect_length(ll, 1)
    expect_true(is.finite(ll))
})

test_that("run_em finite mixture multistart returns the best sampled initialization", {
    X <- matrix(c(
        32, 18,
        28, 22,
        22, 28,
        18, 32
    ), nrow = 4, byrow = TRUE)
    W <- matrix(c(
        35, 15,
        30, 20,
        20, 30,
        15, 35
    ), nrow = 4, byrow = TRUE)

    output <- capture.output({
        fit <- run_em(
            X = X,
            W = W,
            method = "mult",
            mixture = 2,
            initial_prob = matrix(1, nrow = 1, ncol = 1),
            S = 2,
            maxiter = 3,
            miniter = 1,
            maxtime = 5,
            compute_ll = TRUE,
            verbose = TRUE,
            seed = 123
        )
    })

    expect_match(paste(output, collapse = "\n"), "Sampling between 2 samples for an initial probability")
    expect_match(paste(output, collapse = "\n"), "Staying with the biggest log-likelihood sample")
    expect_equal(fit$mixture, 2L)
    expect_equal(dim(fit$initial_prob), c(2L, 2L, 2L))
    expect_equal(length(fit$initial_prob_multistart_logLik), 2L)
    expect_true(all(is.finite(fit$initial_prob_multistart_logLik)))
    expect_true(fit$initial_prob_multistart_best %in% seq_len(2))
})
