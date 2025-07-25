#' fastei: Methods for "A Fast Ecological Inference Algorithm for the R\eqn{\times}C case"
#'
#' Package that implements the methods of [Thraves, C.,Ubilla, P. and Hermosilla, D. (2024): *"A Fast Ecological Inference Algorithm for the R×C Case"*.](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834)
#'
#' Includes a method ([run_em]) to solve the R\eqn{\times}C Ecological Inference problem for the non-parametric case by using the EM algorithm with different approximation methods for the E-Step. The standard deviation of the estimated probabilities can be computed using bootstrapping ([bootstrap]).
#'
#' It also provides a function that generates synthetic election data ([simulate_election]) and a function that imports real election data ([chile_election_2021]) from the Chilean first-round presidential election of 2021.
#'
#' The setting in which the documentation presents the Ecological Inference problem is an election context where for a set of ballot-boxes we observe (i) the votes obtained by each candidate and (ii) the number of voters of each demographic group (for example, these can be defined by age ranges or sex). See [Thraves, C.,Ubilla, P. and Hermosilla, D. (2024): *"A Fast Ecological Inference Algorithm for the R×C Case"*.](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834)
#'
#' The methods to compute the conditional probabilities of the E-Step included in this package are the following:
#' - **Markov Chain Monte Carlo** `(mcmc)`: Performs MCMC to sample vote outcomes for each ballot-box consistent with the observed data. This sample is used to estimate the conditional probability of the *E*-Step.
#' - **Multivariate Normal PDF** `(mvn_pdf)`: Uses the PDF of a Multivariate Normal to approximate the conditional probability.
#' - **Multivariate Normal CDF** `(mvn_cdf)`: Uses the CDF of a Multivariate Normal to approximate the conditional probability.
#' - **Multinomial** `(mult)`: A single Multinomial is used to approximate the sum of Multinomial distributions.
#' - **Exact** `(exact)`: Solves the E-Step exactly using the Total Probability Law, which requires enumerating an exponential number of terms.
#'
#' On average, the **Multinomial** method is the most efficient and precise. Its precision matches the **Exact** method.
#'
#' The documentation uses the following notation:
#' - `b`: number of ballot-boxes.
#' - `g`: number of demographic groups.
#' - `c`: number of candidates.
#' - `a`: number of aggregated macro-groups.
#'
#' To learn more about `fastei`, please consult the available vignettes:
#'
#' \preformatted{
#' browseVignettes("fastei")
#' }
#'
#' @references
#' [Thraves, C., Ubilla, P and Hermosilla D. (2024): *"A Fast Ecological Inference Algorithm for the R×C Case"*.](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834)
#' @aliases NULL
#' @concept package
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
#' @useDynLib fastei, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats pt rgamma rmultinom pnorm
#' @importFrom utils modifyList tail write.csv
NULL
