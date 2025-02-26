# ========= Internals function test ========= #
# Note: .validate_compute will be tested on the compute file.
# ========= .validate_json_eim ========= #

test_that(".validate_json_eim:\tRejects incorrect json paths.", {
    expect_error(fastei:::.validate_json_eim("invalid/path.json"))
    expect_error(fastei:::.validate_json_eim("invalid/file"))
    expect_error(fastei:::.validate_json_eim("~/"))
    expect_error(fastei:::.validate_json_eim(daniel))
})

# ========= .validate_eim ========= #

test_that(".validate_eim:\tValid inputs return TRUE", {
    X <- matrix(1:12, nrow = 3, ncol = 4)
    W <- matrix(13:24, nrow = 3, ncol = 4)
    expect_true(fastei:::.validate_eim(X, W))
})

test_that(".validate_eim:\tError if X is NULL", {
    W <- matrix(13:24, nrow = 3, ncol = 4)
    expect_error(
        fastei:::.validate_eim(NULL, W),
        "Either provide X and W matrices, or a valid JSON path containing them."
    )
})

test_that(".validate_eim:\tError if W is NULL", {
    X <- matrix(1:12, nrow = 3, ncol = 4)
    expect_error(
        fastei:::.validate_eim(X, NULL),
        "Either provide X and W matrices, or a valid JSON path containing them."
    )
})

test_that(".validate_eim:\tError on row mismatch when mismatch = FALSE", {
    X <- matrix(1:8, nrow = 2, ncol = 4)
    W <- matrix(9:20, nrow = 3, ncol = 4)
    expect_error(
        fastei:::.validate_eim(X, W),
        "Mismatch in the number of ballot boxes: 'X' has 2 rows, but 'W' has 3 rows."
    )
})

test_that(".validate_eim:\tPasses on row mismatch when mismatch = TRUE", {
    X <- matrix(1:8, nrow = 2, ncol = 4)
    W <- matrix(9:20, nrow = 3, ncol = 4)
    expect_true(fastei:::.validate_eim(X, W, mismatch = TRUE))
})

test_that(".validate_eim:\tError if candidate matrix X has fewer than 2 columns", {
    X <- matrix(1:3, nrow = 3, ncol = 1)
    W <- matrix(1:6, nrow = 3, ncol = 2)
    expect_error(
        fastei:::.validate_eim(X, W),
        "Candidate matrix 'X' must have at least 2 columns."
    )
})

test_that(".validate_eim:\tError if group matrix W has fewer than 2 columns", {
    X <- matrix(1:6, nrow = 3, ncol = 2)
    W <- matrix(7:9, nrow = 3, ncol = 1)
    expect_error(
        fastei:::.validate_eim(X, W),
        "Group matrix 'W' must have at least 2 columns."
    )
})

test_that(".validate_eim:\tError if X contains NA", {
    X <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 3, ncol = 2)
    W <- matrix(7:12, nrow = 3, ncol = 2)
    expect_error(
        fastei:::.validate_eim(X, W),
        "Matrices 'X' and 'W' cannot contain missing values"
    )
})

test_that(".validate_eim:\tError if W contains NA", {
    X <- matrix(1:6, nrow = 3, ncol = 2)
    W <- matrix(c(7, 8, NA, 10, 11, 12), nrow = 3, ncol = 2)
    expect_error(
        fastei:::.validate_eim(X, W),
        "Matrices 'X' and 'W' cannot contain missing values"
    )
})

test_that(".validate_eim:\tData frames are coerced to matrices", {
    # When X and W are data frames, as.matrix should convert them properly.
    X <- data.frame(a = 1:4, b = 5:8)
    W <- data.frame(c = 9:12, d = 13:16)
    expect_true(fastei:::.validate_eim(X, W))
})

# ========= .random_samples ========= #

test_that(".random_samples:\tThe sampling function has dimensional coherence", {
    tests <- 50
    for (i in 1:tests) {
        # Randomly draw parameters for the simulation:

        result <- fastei:::.random_samples(
            c(5, 50), # Ballots range
            c(2, 5), # Candidates range
            c(2, 5), # Group range
            c(50, 150), # Total votes range
            seed = i
        )

        # Check that for each ballot box, candidate votes (rows of X) match demographic votes (rows of W)
        expect_equal(rowSums(result$X), rowSums(result$W))
        # Check dimensions for candidates and demographic groups
        expect_equal(ncol(result$X), result$candidates)
        expect_equal(ncol(result$W), result$groups)
        # Check that the number of ballot boxes is consistent in both matrices
        expect_equal(nrow(result$X), result$ballots)
        expect_equal(nrow(result$W), result$ballots)
    }
})
