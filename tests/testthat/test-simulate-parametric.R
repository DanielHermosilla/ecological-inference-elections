test_that("simulate_election creates consistent parametric data", {
    num_ballots <- 10
    num_candidates <- 3
    num_groups <- 2
    num_covariates <- 2
    num_districts <- 2
    voters <- 50

    sim <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = voters,
        num_covariates = num_covariates,
        num_districts = num_districts,
        seed = 202
    )

    expect_s3_class(sim, "eim")
    expect_equal(dim(sim$X), c(num_ballots, num_candidates))
    expect_equal(dim(sim$W), c(num_ballots, num_groups))
    expect_equal(dim(sim$V), c(num_ballots, num_covariates))
    expect_equal(dim(sim$real_alpha), c(num_candidates - 1, num_covariates))
    expect_equal(dim(sim$real_beta), c(num_groups, num_candidates - 1))

    expect_true(is.array(sim$real_prob))
    expect_equal(dim(sim$real_prob), c(num_groups, num_candidates, num_ballots))
    expect_prob_array(sim$real_prob)

    expect_equal(rowSums(sim$W), rep(voters, num_ballots))
    expect_equal(rowSums(sim$X), rep(voters, num_ballots))
})

test_that("simulate_election creates matrix-mixture parametric data with covariates", {
    num_ballots <- 12
    num_candidates <- 3
    num_groups <- 2
    num_covariates <- 2
    num_districts <- 3
    mixture <- 2
    voters <- 40

    sim <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = voters,
        num_covariates = num_covariates,
        num_districts = num_districts,
        mixture = mixture,
        seed = 303
    )

    expect_s3_class(sim, "eim")
    expect_equal(dim(sim$X), c(num_ballots, num_candidates))
    expect_equal(dim(sim$W), c(num_ballots, num_groups))
    expect_equal(dim(sim$V), c(num_ballots, num_covariates))
    expect_equal(dim(sim$real_beta), c(num_covariates, mixture - 1))
    expect_equal(dim(sim$real_phi), c(num_ballots, mixture))
    expect_equal(dim(sim$component_prob), c(num_groups, num_candidates, mixture))
    expect_equal(dim(sim$real_prob), c(num_groups, num_candidates, num_ballots))
    expect_equal(dim(sim$latent_matrix), c(num_ballots, mixture))
    expect_equal(dim(sim$outcome), c(num_groups, num_candidates, num_ballots))

    expect_equal(rowSums(sim$real_phi), rep(1, num_ballots), tolerance = 1e-6)
    expect_equal(rowSums(sim$latent_matrix), rep(1, num_ballots))
    expect_true(all(sim$latent_profile %in% seq_len(mixture)))
    expect_prob_array(sim$component_prob)
    expect_prob_array(sim$real_prob)
    expect_equal(rowSums(sim$W), rep(voters, num_ballots))
    expect_equal(rowSums(sim$X), rep(voters, num_ballots))
    expect_equal(unname(apply(sim$outcome, 3, sum)), rep(voters, num_ballots))

    for (b in seq_len(num_ballots)) {
        expect_equal(
            sim$real_prob[, , b],
            sim$component_prob[, , sim$latent_profile[b]],
            tolerance = 1e-12,
            ignore_attr = TRUE
        )
    }
})

test_that("simulate_election supports finite L deviations around latent profiles", {
    sim <- simulate_election(
        num_ballots = 12,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 40,
        num_covariates = 2,
        num_districts = 3,
        mixture = 2,
        L = 5,
        seed = 304
    )

    expect_prob_array(sim$real_prob)
    differs_from_profile <- vapply(seq_len(nrow(sim$X)), function(b) {
        !isTRUE(all.equal(
            sim$real_prob[, , b],
            sim$component_prob[, , sim$latent_profile[b]],
            tolerance = 1e-8,
            check.attributes = FALSE
        ))
    }, logical(1))
    expect_true(any(differs_from_profile))
})

test_that("simulate_election accepts matrix and array component probabilities", {
    prob_matrix <- matrix(
        c(
            2, 1, 1,
            1, 2, 1
        ),
        nrow = 2,
        byrow = TRUE
    )
    sim_matrix <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 30,
        num_covariates = 2,
        num_districts = 2,
        mixture = 2,
        prob = prob_matrix,
        seed = 305
    )

    expected_matrix <- sweep(prob_matrix, 1, rowSums(prob_matrix), "/")
    for (k in seq_len(2)) {
        expect_equal(
            sim_matrix$component_prob[, , k],
            expected_matrix,
            tolerance = 1e-12,
            ignore_attr = TRUE
        )
    }

    prob_array <- array(0, dim = c(2, 3, 2))
    prob_array[, , 1] <- matrix(c(2, 1, 1, 1, 2, 1), nrow = 2, byrow = TRUE)
    prob_array[, , 2] <- matrix(c(1, 1, 2, 2, 1, 1), nrow = 2, byrow = TRUE)
    sim_array <- simulate_election(
        num_ballots = 8,
        num_candidates = 3,
        num_groups = 2,
        ballot_voters = 30,
        num_covariates = 2,
        num_districts = 2,
        mixture = 2,
        prob = prob_array,
        seed = 306
    )

    expect_prob_array(sim_array$component_prob)
    for (k in seq_len(2)) {
        expected_slice <- sweep(prob_array[, , k], 1, rowSums(prob_array[, , k]), "/")
        expect_equal(
            sim_array$component_prob[, , k],
            expected_slice,
            tolerance = 1e-12,
            ignore_attr = TRUE
        )
    }
})

test_that("simulate_election validates matrix-mixture parametric inputs", {
    expect_s3_class(
        simulate_election(
            num_ballots = 6,
            num_candidates = 3,
            num_groups = 2,
            mixture = 2
        ),
        "eim"
    )

    expect_error(
        simulate_election(
            num_ballots = 6,
            num_candidates = 3,
            num_groups = 2,
            num_covariates = 2,
            num_districts = 2,
            mixture = 2,
            L = 0
        ),
        "`L` must be"
    )

    bad_prob <- array(1 / 3, dim = c(2, 3, 3))
    expect_error(
        simulate_election(
            num_ballots = 6,
            num_candidates = 3,
            num_groups = 2,
            num_covariates = 2,
            num_districts = 2,
            mixture = 2,
            prob = bad_prob
        ),
        "array with dimensions"
    )
})
