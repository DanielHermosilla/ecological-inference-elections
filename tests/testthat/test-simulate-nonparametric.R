test_that("simulate_election creates consistent non-parametric data", {
    num_ballots <- 12
    num_candidates <- 3
    num_groups <- 2
    voters <- rep(40, num_ballots)

    sim <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = voters,
        lambda = 0.2,
        seed = 101
    )

    expect_s3_class(sim, "eim")
    expect_equal(dim(sim$X), c(num_ballots, num_candidates))
    expect_equal(dim(sim$W), c(num_ballots, num_groups))
    expect_equal(dim(sim$real_prob), c(num_groups, num_candidates))
    expect_equal(dim(sim$outcome), c(num_groups, num_candidates, num_ballots))

    expect_equal(rowSums(sim$W), voters)
    expect_equal(rowSums(sim$X), voters)

    outcome_by_group <- apply(sim$outcome, c(1, 3), sum)
    outcome_by_candidate <- apply(sim$outcome, c(2, 3), sum)
    expect_equal(outcome_by_group, t(sim$W), ignore_attr = TRUE)
    expect_equal(outcome_by_candidate, t(sim$X), ignore_attr = TRUE)

    expect_prob_matrix(sim$real_prob)
})

test_that("simulate_election creates non-parametric finite mixtures", {
    num_ballots <- 12
    num_candidates <- 3
    num_groups <- 2
    mixture <- 2
    voters <- rep(40, num_ballots)

    sim <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = voters,
        lambda = 0.2,
        seed = 102,
        mixture = mixture
    )

    expect_s3_class(sim, "eim")
    expect_equal(dim(sim$X), c(num_ballots, num_candidates))
    expect_equal(dim(sim$W), c(num_ballots, num_groups))
    expect_equal(dim(sim$component_prob), c(num_groups, num_candidates, mixture))
    expect_equal(dim(sim$real_prob), c(num_groups, num_candidates, num_ballots))
    expect_equal(dim(sim$V), c(num_ballots, mixture))
    expect_equal(dim(sim$real_phi), c(num_ballots, mixture))
    expect_equal(dim(sim$latent_matrix), c(num_ballots, mixture))
    expect_equal(dim(sim$outcome), c(num_groups, num_candidates, num_ballots))

    expect_null(sim$real_beta)
    expect_equal(rowSums(sim$W), voters)
    expect_equal(rowSums(sim$X), voters)
    expect_equal(rowSums(sim$real_phi), rep(1, num_ballots), tolerance = 1e-6)
    expect_equal(sim$real_phi, sim$latent_matrix, ignore_attr = TRUE)
    expect_equal(sim$V, sim$latent_matrix, ignore_attr = TRUE)
    expect_equal(rowSums(sim$latent_matrix), rep(1, num_ballots))
    expect_true(all(sim$real_phi %in% c(0, 1)))
    expect_true(all(sim$latent_profile %in% seq_len(mixture)))
    expect_prob_array(sim$component_prob)
    expect_prob_array(sim$real_prob)

    for (b in seq_len(num_ballots)) {
        expect_equal(
            sim$real_prob[, , b],
            sim$component_prob[, , sim$latent_profile[b]],
            tolerance = 1e-12,
            ignore_attr = TRUE
        )
    }
})
