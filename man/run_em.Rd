% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eim-class.R
\name{run_em}
\alias{run_em}
\alias{run_em()}
\title{Compute the Expected-Maximization Algorithm}
\usage{
run_em(
  object = NULL,
  X = NULL,
  W = NULL,
  json_path = NULL,
  method = "mult",
  initial_prob = "group_proportional",
  allow_mismatch = TRUE,
  maxiter = 1000,
  miniter = 0,
  maxtime = 3600,
  param_threshold = 0.001,
  ll_threshold = as.double(-Inf),
  seed = NULL,
  verbose = FALSE,
  group_agg = NULL,
  mcmc_samples = 1000,
  mcmc_stepsize = 3000,
  mvncdf_method = "genz",
  mvncdf_error = 0.00001,
  mvncdf_samples = 5000,
  metropolis_iter = 5,
  burn_in = 10000,
  ...
)
}
\arguments{
\item{object}{An object of class \code{eim}, which can be created using the \link{eim} function. This parameter should not be used if either (i) \code{X} and \code{W} matrices or (ii) \code{json_path} is supplied. See \strong{Note}.}

\item{X}{A \verb{(b x c)} matrix representing candidate votes per ballot box.}

\item{W}{A \verb{(b x g)} matrix representing group votes per ballot box.}

\item{json_path}{A path to a JSON file containing \code{X} and \code{W} fields, stored as nested arrays. It may contain additional fields with other attributes, which will be added to the returned object.}

\item{method}{An optional string specifying the method used for estimating the E-step. Valid
options are:
\itemize{
\item \code{mult}: The default method, using a single sum of Multinomial distributions.
\item \code{mvn_cdf}: Uses a Multivariate Normal CDF distribution to approximate the conditional probability.
\item \code{mvn_pdf}: Uses a Multivariate Normal PDF distribution to approximate the conditional probability.
\item \code{mcmc}: Uses MCMC to sample vote outcomes. This is used to estimate the conditional probability of the E-step.
\item \code{metropolis}: Uses the Metropolis-Hastings algorithm to sample vote outcomes. This is used to estimate the conditional probability of the E-step.
\item \code{exact}: Solves the E-step using the Total Probability Law.
}

For a detailed description of each method, see \link{fastei-package} and \strong{References}.}

\item{initial_prob}{An optional string specifying the method used to obtain the initial
probability. Accepted values are:
\itemize{
\item \code{uniform}: Assigns equal probability to every candidate within each group.
\item \code{proportional}: Assigns probabilities to each group based on the proportion of candidates votes.
\item \code{group_proportional}: Computes the probability matrix by taking into account both group and candidate proportions. This is the default method.
\item \code{random}: Use randomized values to fill the probability matrix.
}}

\item{allow_mismatch}{Boolean, if \code{TRUE}, allows a mismatch between the voters and votes for each ballot-box, only works if \code{method} is \code{"mvn_cdf"}, \code{"mvn_pdf"}, \code{"mult"} and \code{"mcmc"}. If \code{FALSE}, throws an error if there is a mismatch. By default it is \code{TRUE}.}

\item{maxiter}{An optional integer indicating the maximum number of EM iterations.
The default value is \code{1000}.}

\item{miniter}{An optional integer indicating the minimum number of EM iterations. The default value is \code{0}.}

\item{maxtime}{An optional numeric specifying the maximum running time (in seconds) for the
algorithm. This is checked at every iteration of the EM algorithm. The default value is \code{3600}, which corresponds to an hour.}

\item{param_threshold}{An optional numeric value indicating the minimum difference between
consecutive probability values required to stop iterating. The default value is \code{0.001}. Note that the algorithm will stop if either \code{ll_threshold} \strong{or} \code{param_threshold} is accomplished.}

\item{ll_threshold}{An optional numeric value indicating the minimum difference between consecutive log-likelihood values to stop iterating. The default value is \code{inf}, essentially deactivating
the threshold. Note that the algorithm will stop if either \code{ll_threshold} \strong{or} \code{param_threshold} is accomplished.}

\item{seed}{An optional integer indicating the random seed for the randomized algorithms. This argument is only applicable if \code{initial_prob = "random"} or \code{method} is either \code{"mcmc"} or \code{"mvn_cdf"}.}

\item{verbose}{An optional boolean indicating whether to print informational messages during the EM
iterations. The default value is \code{FALSE}.}

\item{group_agg}{An optional vector of increasing integers from 1 to the number of columns in \code{W}, specifying how to aggregate groups in \code{W} before running the EM algorithm. Each value represents the highest column index included in each aggregated group. For example, if \code{W} has four columns, \code{group_agg = c(2, 4)} indicates that columns 1 and 2 should be combined into one group, and columns 3 and 4 into another. Defaults to \code{NULL}, in which case no group aggregation is performed.}

\item{mcmc_samples}{An optional integer indicating the number of samples to generate for the
\strong{MCMC} method. This parameter is only relevant when \code{method = "mcmc"}.
The default value is \code{1000}.}

\item{mcmc_stepsize}{An optional integer specifying the step size for the \code{mcmc}
algorithm. This parameter is only applicable when \code{method = "mcmc"} and will
be ignored otherwise. The default value is \code{3000}.}

\item{mvncdf_method}{An optional string specifying the method used to estimate the \code{mvn_cdf} method
via a Monte Carlo simulation. Accepted values are \code{genz} and \code{genz2}, with \code{genz}
set as the default. This parameter is only applicable when \code{method = "mvn_cdf"}. See \strong{References} for more details.}

\item{mvncdf_error}{An optional numeric value defining the error threshold for the Monte Carlo
simulation when estimating the \code{mvn_cdf} method. The default value is \code{1e-6}. This parameter is only relevant
when \code{method = "mvn_cdf"}.}

\item{mvncdf_samples}{An optional integer specifying the number of Monte Carlo
samples for the \code{mvn_cdf} method. The default value is \code{5000}. This argument is only applicable when \code{method = "mvn_cdf"}.}

\item{metropolis_iter}{An optional integer specifying the amount of iterations to run the Metropolis-Hastings algorithm for the \code{metropolis} method. The default value is \code{5}. This argument is only applicable when \code{method = "metropolis"}.}

\item{burn_in}{Amount of iterations to discard before starting the MCMC sampling. This is only applicable when \code{method = "mcmc"} or \code{method = "metropolis"}. The default value is \code{10000}.}

\item{...}{Added for compability}
}
\value{
The function returns an \code{eim} object with the function arguments and the following attributes:
\describe{
\item{prob}{The estimated probability matrix \verb{(g x c)}.}
\item{cond_prob}{A \verb{(g x c x b)} 3d-array with the probability that a at each ballot-box a voter of each group voted for each candidate, given the observed outcome at the particular ballot-box.}
\item{expected_outcome}{A \verb{(g x c x b)} 3d-array with the expected votes cast for each ballot box.}
\item{logLik}{The log-likelihood value from the last iteration.}
\item{iterations}{The total number of iterations performed by the EM algorithm.}
\item{time}{The total execution time of the algorithm in seconds.}
\item{status}{
The final status ID of the algorithm upon completion:
\itemize{
\item \code{0}: Converged
\item \code{1}: Maximum time reached.
\item \code{2}: Maximum iterations reached.
}
}
\item{message}{The finishing status displayed as a message, matching the status ID value.}
\item{method}{The method for estimating the conditional probability in the E-step.}
}
Aditionally, it will create \code{mcmc_samples} and \code{mcmc_stepsize} parameters if the specified \code{method = "mcmc"}, or \code{mvncdf_method}, \code{mvncdf_error} and \code{mvncdf_samples} if \code{method = "mvn_cdf"}.

Also, if the eim object supplied is created with the function \link{simulate_election}, it also returns the real probability with the name \code{real_prob}. See \link{simulate_election}.

If \code{group_agg} is different than \code{NULL}, two values are returned: \code{W_agg} a \verb{(b x a)} matrix with the number of voters of each aggregated group o each ballot-box, and \code{group_agg} the same input vector.
}
\description{
Executes the Expectation-Maximization (EM) algorithm indicating the approximation method to use in the E-step.
Certain methods may require additional arguments, which can be passed through \code{...} (see \link{fastei-package} for more details).
}
\note{
This function can be executed using one of three mutually exclusive approaches:
\enumerate{
\item By providing an existing \code{eim} object.
\item By supplying both input matrices (\code{X} and \code{W}) directly.
\item By specifying a JSON file (\code{json_path}) containing the matrices.
}

These input methods are \strong{mutually exclusive}, meaning that you must provide \strong{exactly one} of
these options. Attempting to provide more than one or none of these inputs will result in an error.

When called with an \code{eim} object, the function updates the object with the computed results.
If an \code{eim} object is not provided, the function will create one internally using either the
supplied matrices or the data from the JSON file before executing the algorithm.
}
\examples{
\donttest{
# Example 1: Compute the Expected-Maximization with default settings
simulations <- simulate_election(
    num_ballots = 300,
    num_candidates = 5,
    num_groups = 3,
)
model <- eim(simulations$X, simulations$W)
model <- run_em(model) # Returns the object with updated attributes

# Example 2: Compute the Expected-Maximization using the mvn_pdf method
model <- run_em(
    object = model,
    method = "mvn_pdf",
)

# Example 3: Run the mvn_cdf method with default settings
model <- run_em(object = model, method = "mvn_cdf")
}
\dontrun{
# Example 4: Perform an Exact estimation using user-defined parameters

run_em(
    json_path = "a/json/file.json",
    method = "exact",
    initial_prob = "uniform",
    maxiter = 10,
    maxtime = 600,
    param_threshold = 1e-3,
    ll_threshold = 1e-5,
    verbose = TRUE
)
}

}
\references{
\href{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4832834}{Thraves, C., Ubilla, P. and Hermosilla, D.: \emph{"Fast Ecological Inference Algorithm for the RxC Case"}}. Aditionally, the MVN CDF is computed by the methods introduced in \href{https://www.researchgate.net/publication/2463953_Numerical_Computation_Of_Multivariate_Normal_Probabilities}{Genz, A. (2000). Numerical computation of multivariate normal probabilities. \emph{Journal of Computational and Graphical Statistics}}
}
\seealso{
The \link{eim} object implementation.
}
