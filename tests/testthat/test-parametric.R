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
        method = "mvn_pdf",
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
        ballot_voters = rep(30, 6),
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
    expect_equal(dim(fit_pdf$cond_prob), c(2, 3, 6))
    expect_equal(dim(fit_cdf$cond_prob), c(2, 3, 6))
    expect_true(all(is.finite(fit_pdf$prob)))
    expect_true(all(is.finite(fit_cdf$prob)))
    expect_true(all(is.finite(fit_pdf$cond_prob)))
    expect_true(all(is.finite(fit_cdf$cond_prob)))
    expect_equal(unname(apply(fit_pdf$cond_prob, c(1, 3), sum)), matrix(1, 2, 6), tolerance = 1e-6)
    expect_equal(unname(apply(fit_cdf$cond_prob, c(1, 3), sum)), matrix(1, 2, 6), tolerance = 1e-6)
    expect_equal(dim(fit_pdf$beta), c(2, 2))
    expect_equal(dim(fit_cdf$alpha), c(2, 2))

    expect_error(
        run_em(object = model, method = "mcmc"),
        "Parametric mode only supports"
    )
})

test_that("parametric run_em supports the exact E-step", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = rep(8, 4),
        num_covariates = 1,
        num_districts = 2,
        seed = 204
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "exact",
        maxiter = 2,
        maxtime = 5,
        maxnewton = 1,
        adjust_prob_cond_method = ""
    )

    expect_equal(dim(fit$prob), c(2, 2, 4))
    expect_equal(dim(fit$cond_prob), c(2, 2, 4))
    expect_true(all(is.finite(fit$prob)))
    expect_true(all(is.finite(fit$cond_prob)))
    expect_equal(unname(apply(fit$cond_prob, c(1, 3), sum)), matrix(1, 2, 4), tolerance = 1e-6)

    W_mismatch <- sim$W
    W_mismatch[1, 1] <- W_mismatch[1, 1] + 1
    fit_mismatch <- suppressMessages(run_em(
        X = sim$X,
        W = W_mismatch,
        V = sim$V,
        method = "exact",
        maxiter = 1,
        maxtime = 5,
        maxnewton = 1,
        allow_mismatch = TRUE,
        adjust_prob_cond_method = ""
    ))
    expect_equal(rowSums(fit_mismatch$W), rowSums(fit_mismatch$X))
})

test_that("joint parametric EM supports mvn and exact E-steps", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 2,
        num_groups = 2,
        ballot_voters = rep(8, 4),
        num_covariates = 1,
        num_districts = 2,
        seed = 205
    )

    for (method in c("mvn_pdf", "mvn_cdf", "exact")) {
        fit <- run_em(
            X = sim$X,
            W = sim$W,
            V = sim$V,
            method = method,
            symmetric = TRUE,
            symmetric_weight_method = "joint",
            maxiter = 1,
            maxtime = 5,
            maxnewton = 1,
            mvncdf_samples = 30,
            mvncdf_error = 1e-2,
            adjust_prob_cond_method = ""
        )

        expect_equal(dim(fit$prob), c(2, 2, 4), info = method)
        expect_true(all(is.finite(fit$prob)), info = method)
        expect_true(all(is.finite(fit$cond_prob)), info = method)
    }
})
