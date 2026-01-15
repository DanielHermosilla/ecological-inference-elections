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
