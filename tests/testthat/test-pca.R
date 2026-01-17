test_that("PCA reduces V with a fixed number of components", {
    sim <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 3,
        num_districts = 2,
        seed = 220
    )

    rownames(sim$V) <- paste0("B", seq_len(nrow(sim$V)))

    reduced <- PCA(sim, components = 2)

    expect_s3_class(reduced, "eim")
    expect_equal(dim(reduced$V), c(nrow(sim$V), 2))
    expect_equal(colnames(reduced$V), c("PCA 1", "PCA 2"))
    expect_equal(rownames(reduced$V), rownames(sim$V))
})

test_that("PCA selects components using sd_threshold", {
    sim <- simulate_election(
        num_ballots = 10,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 4,
        num_districts = 2,
        seed = 221
    )

    threshold <- 0.85
    pca <- stats::prcomp(sim$V, center = TRUE, scale. = TRUE)
    var_ratio <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    expected_k <- which(var_ratio >= threshold)[1]

    reduced <- PCA(sim, sd_threshold = threshold)

    expect_equal(ncol(reduced$V), expected_k)
    expect_equal(colnames(reduced$V), paste("PCA", seq_len(expected_k)))
})

test_that("PCA validates parametric inputs", {
    sim_np <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = rep(30, 6),
        seed = 222
    )

    expect_error(
        PCA(sim_np, components = 2),
        "parametric object"
    )

    sim_p <- simulate_election(
        num_ballots = 6,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 30,
        num_covariates = 3,
        num_districts = 2,
        seed = 223
    )

    expect_error(
        PCA(sim_p, components = 2, sd_threshold = 0.9),
        "either"
    )

    expect_error(
        PCA(sim_p, components = 1.5),
        "components"
    )

    expect_error(
        PCA(sim_p, sd_threshold = 1.2),
        "sd_threshold"
    )
})
