test_that("simulate_election returns parametric data", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 50,
        num_covariates = 2,
        num_districts = 2,
        seed = 200
    )

    expect_s3_class(sim, "eim")
    expect_true(is.matrix(sim$V))
    expect_equal(dim(sim$V), c(8, 2))
    expect_true(is.matrix(sim$real_alpha))
    expect_true(is.matrix(sim$real_beta))
    expect_equal(dim(sim$real_alpha), c(2, 2))
    expect_equal(dim(sim$real_beta), c(2, 2))
})

test_that("run_em works in parametric mode", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 201
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    fit <- run_em(
        object = model,
        method = "mult",
        beta_init = sim$real_beta,
        alpha_init = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    expect_s3_class(fit, "eim")
    expect_true(is.array(fit$prob))
    expect_equal(dim(fit$prob), c(2, 3, 8))
    expect_equal(dim(fit$beta), c(2, 2))
    expect_equal(dim(fit$alpha), c(2, 2))
})

test_that("bootstrap works in parametric mode", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 202
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    boot <- bootstrap(
        object = model,
        nboot = 3,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    expect_s3_class(boot, "eim")
    expect_true(is.matrix(boot$sd_beta))
    expect_true(is.matrix(boot$sd_alpha))
    expect_equal(dim(boot$sd_beta), c(2, 2))
    expect_equal(dim(boot$sd_alpha), c(2, 2))
})

test_that("parametric run_em supports mvn methods", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        num_covariates = 2,
        num_districts = 2,
        seed = 203
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    fit_pdf <- run_em(
        object = model,
        method = "mvn_pdf",
        maxiter = 2,
        maxtime = 2,
        maxnewton = 1
    )
    fit_cdf <- run_em(
        object = model,
        method = "mvn_cdf",
        maxiter = 2,
        maxtime = 2,
        maxnewton = 1,
        mvncdf_samples = 50,
        mvncdf_error = 1e-2
    )

    expect_equal(dim(fit_pdf$prob), c(2, 3, 6))
    expect_equal(dim(fit_cdf$prob), c(2, 3, 6))
    expect_equal(dim(fit_pdf$beta), c(2, 2))
    expect_equal(dim(fit_cdf$alpha), c(2, 2))
})

test_that("parametric run_em rejects unsupported methods", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        num_covariates = 2,
        num_districts = 2,
        seed = 204
    )

    model <- eim(X = sim$X, W = sim$W, V = sim$V)
    expect_error(
        run_em(object = model, method = "exact"),
        "Parametric mode only supports"
    )
})

test_that("run_em supports parametric finite mixtures", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 205
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        mixture = 2,
        beta_init = matrix(0, nrow = ncol(sim$V), ncol = 1),
        maxiter = 2,
        maxtime = 3,
        maxnewton = 1,
        compute_ll = TRUE
    )

    expect_s3_class(fit, "eim")
    expect_equal(fit$mixture, 2L)
    expect_equal(fit$K, 2L)
    expect_equal(dim(fit$phi), c(nrow(sim$X), 2))
    expect_equal(rowSums(fit$phi), rep(1, nrow(sim$X)), tolerance = 1e-6)
    expect_equal(dim(fit$responsibilities), c(nrow(sim$X), 2))
    expect_equal(rowSums(fit$responsibilities), rep(1, nrow(sim$X)), tolerance = 1e-6)
    expect_equal(dim(fit$prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$cond_prob), c(ncol(sim$W), ncol(sim$X), nrow(sim$X)))
    expect_equal(dim(fit$component_prob), c(ncol(sim$W), ncol(sim$X), 2))
    expect_equal(dim(fit$beta), c(ncol(sim$V), 1))
    expect_null(fit$alpha)
})

test_that("parametric matrix mixture initializes components from group proportional base", {
    object <- list(
        X = matrix(
            c(
                80, 20,
                20, 80
            ),
            nrow = 2,
            byrow = TRUE
        )
    )
    W_matrix <- matrix(
        c(
            90, 10,
            10, 90
        ),
        nrow = 2,
        byrow = TRUE
    )
    control <- list(K = 2, initial_prob = "group_proportional")

    init <- .run_em_matrix_mixture_component_prob_array(object, W_matrix, control)

    group_weights <- sweep(W_matrix, 1, rowSums(W_matrix), "/")
    base <- crossprod(group_weights, object$X)
    base <- sweep(base, 1, rowSums(base), "/")

    expect_equal(dim(init), c(2, 2, 2))
    for (k in seq_len(2)) {
        jitter <- matrix(
            1 + 0.05 * sin(seq_len(4) + k),
            nrow = 2,
            ncol = 2
        )
        expected <- base * jitter
        expected <- sweep(expected, 1, rowSums(expected), "/")
        expect_equal(init[, , k], expected, tolerance = 1e-12, ignore_attr = TRUE)
    }
    expect_gt(abs(init[1, 1, 1] - init[2, 1, 1]), 0.3)
})

test_that("parametric matrix mixture updates membership beta", {
    sim <- simulate_election(
        num_ballots = 80,
        num_candidates = 2,
        num_groups = 5,
        ballot_voters = 80,
        lambda = 0.5,
        num_covariates = 2,
        num_districts = 2,
        mixture = 3,
        seed = 12345
    )

    beta_init <- matrix(0, nrow = ncol(sim$V), ncol = 2)
    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mvn_pdf",
        mixture = 3,
        beta_init = beta_init,
        maxiter = 2,
        maxtime = 10,
        maxnewton = 10,
        compute_ll = TRUE,
        seed = 12345
    )

    expect_equal(dim(fit$beta), dim(beta_init))
    expect_gt(max(abs(fit$beta - beta_init)), 1e-4)
    expect_false(isTRUE(all.equal(fit$phi, matrix(1 / 3, nrow = nrow(sim$X), ncol = 3))))
})

test_that("parametric matrix mixtures support mvn methods and component probability initial values", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 2,
        seed = 206
    )

    beta_init <- matrix(0, nrow = ncol(sim$V), ncol = 1)
    initial_prob <- array(0, dim = c(ncol(sim$W), ncol(sim$X), 2))
    initial_prob[, , 1] <- apply(sim$real_prob, c(1, 2), mean)
    initial_prob[, , 2] <- initial_prob[, , 1] * matrix(c(1.05, 0.95, 1.02, 0.98, 1.01, 0.99), nrow = ncol(sim$W))
    initial_prob[, , 2] <- sweep(initial_prob[, , 2], 1, rowSums(initial_prob[, , 2]), "/")

    fit_pdf <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mvn_pdf",
        mixture = 2,
        beta_init = beta_init,
        initial_prob = initial_prob,
        maxiter = 2,
        maxtime = 3,
        maxnewton = 1
    )
    fit_cdf <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mvn_cdf",
        mixture = 2,
        beta_init = beta_init,
        initial_prob = initial_prob,
        maxiter = 2,
        maxtime = 3,
        maxnewton = 1,
        mvncdf_samples = 50,
        mvncdf_error = 1e-2
    )

    expect_equal(dim(fit_pdf$component_prob), dim(initial_prob))
    expect_equal(dim(fit_cdf$component_prob), dim(initial_prob))
    expect_null(fit_pdf$alpha)
    expect_equal(dim(fit_cdf$beta), dim(beta_init))
    expect_equal(rowSums(fit_pdf$phi), rep(1, nrow(sim$X)), tolerance = 1e-6)
})

test_that("parametric matrix mixtures respect minimum EM iterations", {
    sim <- simulate_election(
        num_ballots = 30,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 25,
        mixture = 2,
        num_covariates = 2,
        num_districts = 3,
        seed = 207
    )

    fit <- suppressWarnings(run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        mixture = 2,
        maxiter = 6,
        miniter = 5,
        param_threshold = 1e9,
        maxnewton = 1
    ))

    expect_gte(fit$iterations, 5)
})

test_that("parametric matrix mixtures support multistart initial probabilities", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = 30,
        mixture = 2,
        num_covariates = 2,
        num_districts = 2,
        seed = 208
    )

    output <- capture.output(
        fit <- run_em(
            X = sim$X,
            W = sim$W,
            V = sim$V,
            method = "mult",
            mixture = 2,
            S = 2,
            maxiter = 2,
            maxnewton = 1,
            verbose = TRUE
        )
    )

    expect_match(paste(output, collapse = "\n"), "Sampling between 2 samples for an initial probability")
    expect_match(paste(output, collapse = "\n"), "Staying with the biggest log-likelihood sample")
    expect_equal(dim(fit$initial_prob), c(2L, 2L, 2L))
    expect_equal(length(fit$initial_prob_multistart_logLik), 2L)
    expect_true(all(is.finite(fit$initial_prob_multistart_logLik)))
    expect_true(fit$initial_prob_multistart_best %in% seq_len(2))
    expect_equal(dim(fit$component_prob), c(2L, 2L, 2L))

    selected_ll <- EMLogLikFromProbParametricMixture(
        as.matrix(sim$X),
        sim$W,
        sim$V,
        matrix(0, nrow = ncol(sim$V), ncol = 1),
        fit$initial_prob,
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
    expect_equal(max(fit$initial_prob_multistart_logLik), selected_ll, tolerance = 1e-8)
})

test_that("parametric run_em rejects row mixtures", {
    sim <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        num_covariates = 2,
        num_districts = 2,
        seed = 207
    )

    expect_error(
        run_em(X = sim$X, W = sim$W, V = sim$V, row_mixture = 2),
        "row_mixture > 1"
    )
})
