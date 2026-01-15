test_that("save_eim writes JSON, RDS, and CSV for non-parametric models", {
    sim <- simulate_election(
        num_ballots = 5,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(20, 5),
        seed = 200
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        method = "mult",
        maxiter = 3,
        maxtime = 2,
        compute_ll = FALSE
    )

    out_rds <- tempfile(fileext = ".rds")
    out_json <- tempfile(fileext = ".json")
    out_csv <- tempfile(fileext = ".csv")

    save_eim(fit, out_rds)
    save_eim(fit, out_json)
    save_eim(fit, out_csv)

    expect_true(file.exists(out_rds))
    expect_true(file.exists(out_json))
    expect_true(file.exists(out_csv))

    json_data <- jsonlite::fromJSON(out_json)
    expect_true(all(c("X", "W", "prob") %in% names(json_data)))

    reloaded <- eim(json_path = out_json)
    expect_s3_class(reloaded, "eim")
    expect_equal(dim(reloaded$X), dim(fit$X))
    expect_equal(dim(reloaded$W), dim(fit$W))

    csv_data <- read.csv(out_csv)
    expect_equal(nrow(csv_data), nrow(fit$prob))
    expect_equal(ncol(csv_data), ncol(fit$prob) + 1)
})

test_that("save_eim rejects CSV export for parametric probabilities", {
    sim <- simulate_election(
        num_ballots = 4,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 20,
        parametric = TRUE,
        num_attributes = 2,
        num_districts = 2,
        seed = 201
    )

    fit <- run_em(
        X = sim$X,
        W = sim$W,
        V = sim$V,
        method = "mult",
        beta = sim$real_beta,
        alpha = sim$real_alpha,
        maxiter = 3,
        maxtime = 2,
        maxnewton = 1
    )

    out_csv <- tempfile(fileext = ".csv")
    expect_error(
        save_eim(fit, out_csv),
        "CSV export is not supported"
    )
})
