mvn_pdf_reference_loglik <- function(X, W, probabilities) {
    B <- nrow(X)
    G <- ncol(W)
    C <- ncol(X)
    n <- C - 1L

    sum(vapply(seq_len(B), function(b) {
        p_star <- probabilities[, seq_len(n), drop = FALSE]
        scale <- if (sum(W[b, ]) > 0) sum(X[b, ]) / sum(W[b, ]) else 1
        effective_w <- W[b, ] * scale
        mu <- colSums(sweep(p_star, 1L, effective_w, `*`))
        sigma <- matrix(0, nrow = n, ncol = n)
        for (g in seq_len(G)) {
            sigma <- sigma + effective_w[g] * (
                diag(p_star[g, ], nrow = n) - tcrossprod(p_star[g, ])
            )
        }

        delta <- X[b, seq_len(n)] - mu
        chol_sigma <- chol(sigma)
        standardized <- backsolve(chol_sigma, delta, transpose = TRUE)
        -0.5 * (
            n * log(2 * pi) +
                2 * sum(log(diag(chol_sigma))) +
                sum(standardized^2)
        )
    }, numeric(1)))
}

evaluate_mvn_pdf_loglik <- function(X, W, probabilities) {
    EMLogLikFromProb(
        t(X), W, probabilities, "mvn_pdf",
        3000L, 1000L, "genz", 1e-3, 1000L, 0L,
        "project_lp", FALSE
    )
}

test_that("mvn_pdf uses the unconditional Gaussian likelihood", {
    X <- matrix(c(
        43, 27, 30,
        39, 31, 30
    ), nrow = 2, byrow = TRUE)
    W <- matrix(c(
        40, 60,
        55, 45
    ), nrow = 2, byrow = TRUE)
    probabilities <- matrix(c(
        0.2, 0.5, 0.3,
        0.6, 0.1, 0.3
    ), nrow = 2, byrow = TRUE)

    expect_equal(
        evaluate_mvn_pdf_loglik(X, W, probabilities),
        mvn_pdf_reference_loglik(X, W, probabilities),
        tolerance = 1e-10
    )
})

test_that("mvn_pdf likelihood scales W for mismatched ballot totals", {
    X <- matrix(c(24, 22, 14), nrow = 1)
    W <- matrix(c(30, 20), nrow = 1)
    probabilities <- matrix(c(
        0.55, 0.30, 0.15,
        0.20, 0.50, 0.30
    ), nrow = 2, byrow = TRUE)

    expect_equal(
        evaluate_mvn_pdf_loglik(X, W, probabilities),
        mvn_pdf_reference_loglik(X, W, probabilities),
        tolerance = 1e-10
    )
})

test_that("mvn_pdf likelihood handles binary outcomes", {
    X <- matrix(c(24, 26, 20, 30), nrow = 2, byrow = TRUE)
    W <- matrix(c(30, 20, 25, 25), nrow = 2, byrow = TRUE)
    probabilities <- matrix(c(
        0.65, 0.35,
        0.25, 0.75
    ), nrow = 2, byrow = TRUE)

    expect_equal(
        evaluate_mvn_pdf_loglik(X, W, probabilities),
        mvn_pdf_reference_loglik(X, W, probabilities),
        tolerance = 1e-10
    )
})

test_that("mvn_pdf likelihood respects degenerate Gaussian support", {
    W <- matrix(c(30, 20), nrow = 1)
    probabilities <- matrix(c(
        1, 0, 0,
        0, 1, 0
    ), nrow = 2, byrow = TRUE)

    expect_equal(
        evaluate_mvn_pdf_loglik(matrix(c(30, 20, 0), nrow = 1), W, probabilities),
        0
    )
    expect_equal(
        evaluate_mvn_pdf_loglik(matrix(c(29, 21, 0), nrow = 1), W, probabilities),
        -Inf
    )
})
