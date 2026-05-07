test_that("save_eim writes JSON, RDS, and CSV", {
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
