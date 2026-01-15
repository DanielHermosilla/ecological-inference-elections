expect_prob_matrix <- function(mat, tol = 1e-6) {
    testthat::expect_true(is.matrix(mat))
    testthat::expect_true(all(is.finite(mat)))
    testthat::expect_true(all(mat >= -tol))
    testthat::expect_true(all(mat <= 1 + tol))
    rs <- rowSums(mat)
    testthat::expect_true(all(abs(rs - 1) < tol))
}

expect_prob_array <- function(arr, tol = 1e-6) {
    testthat::expect_true(is.array(arr))
    testthat::expect_true(length(dim(arr)) == 3)
    testthat::expect_true(all(is.finite(arr)))
    testthat::expect_true(all(arr >= -tol))
    testthat::expect_true(all(arr <= 1 + tol))
    sums <- apply(arr, c(1, 3), sum)
    testthat::expect_true(all(abs(sums - 1) < tol))
}

expect_dimnames_match <- function(arr, group_names, cand_names, ballot_names) {
    testthat::expect_equal(dimnames(arr)[[1]], group_names)
    testthat::expect_equal(dimnames(arr)[[2]], cand_names)
    testthat::expect_equal(dimnames(arr)[[3]], ballot_names)
}

mae <- function(x, y) {
    mean(abs(x - y))
}

expected_votes_from_q <- function(W, q) {
    B <- nrow(W)
    C <- dim(q)[2]
    out <- matrix(0, nrow = B, ncol = C)
    for (b in seq_len(B)) {
        out[b, ] <- as.numeric(W[b, ] %*% q[, , b])
    }
    out
}
