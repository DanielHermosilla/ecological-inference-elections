#' eimodel: A Fast Alternative for the R x C Ecological Inference Problem
#'
#' Addressing the \eqn{R \times C} ecological inference problem is challenging due to
#' the variety of possible approaches, including parametric probability models,
#' entropy maximization, or mathematical programming. This package introduces a
#' nonparametric model, based on the Expectation-Maximization (EM) algorithm, to maximize
#' the likelihood given observed data. Although the *M*-step can be solved
#' in closed form, the *E*-step can involve an exponential number of computations.
#'
#' The `eimodel` package provides multiple methods to estimate the *E*-step
#' in an EM algorithm for the standard \eqn{R \times C} ecological inference problem.
#' It also offers a function to simulate elections under diverse parameters. These
#' methods, which can be accessed via the `eim` object, include:
#'
#' - `Hit and Run (H&R)`: Simulates a subset of samples to approximate the exact solution.
#' - `Multivariate Normal PDF (MVN PDF)`: Approximates the *E*-step using a multivariate normal PDF.
#' - `Multivariate Normal CDF (MVN CDF)`: Uses the CDF of a multivariate normal distribution to approximate probabilities.
#' - `Multinomial (Mult)`: Approximates the *E*-step with a multinomial distribution.
#' - `Exact`: Solves the determinant approach directly.
#'
#' On average, the `Multinomial` method tends to be the most efficient and precise,
#' often yielding the lowest mean absolute error (MAE). Its runtime
#' is roughly on the order of \eqn{10^{-5}} seconds.
#'
#' To learn more about `eimodel`, please consult the available vignettes:
#'
#' \preformatted{
#' browseVignettes("eimodel")
#' }
#'
#' @references
#' [Thraves, C. and Ubilla, P.: *"Fast Ecological Inference Algorithm for the R×C Case"*.](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834)
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
#' @useDynLib eimodel, .registration = TRUE
# @importFrom Rcpp sourceCpp
NULL
