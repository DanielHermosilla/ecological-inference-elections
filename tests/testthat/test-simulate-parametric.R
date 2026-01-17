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
