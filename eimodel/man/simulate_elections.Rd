\name{simulate_elections}
\alias{simulate_elections}
\title{Simulates an election by creating the candidate and group matrices and their results.}
\usage{
simulate_elections(
  num_ballots,
  num_candidates,
  num_groups,
  ballot_voters,
  lambda = 0.5,
  seed = NULL
)
}
\arguments{
\item{num_ballots}{\emph{(integer(1))} The number of ballot boxes \emph{("b")}.}
\item{num_candidates}{\emph{(integer(1))} The number of candidates \emph{("c")}.}
\item{num_groups}{\emph{(integer(1))} The number of demographic groups \emph{("g")}.}
\item{ballot_voters}{\emph{(integer(num_ballots))} A vector of length \code{num_ballots} with the number of votes per ballot box.}
\item{lambda}{\emph{(numeric(1))} A value between 0 and 1 representing the heterogeneity of the groups. Values near 0 yield more heterogeneous results, but may be less realistic; the opposite is true for values closer to 1. (default: 0.5)}
\item{seed}{\emph{(numeric(1))} Optional. If provided, it overrides the current global seed. (default: \code{NULL})}
}
\value{
A list with components:
\item{X}{A matrix (b x c) with candidate votes per ballot box.}
\item{W}{A matrix (b x g) with demographic votes per ballot box.}
\item{real_p}{A matrix (g x c) with the estimated (unobserved) probabilities that a demographic group votes for a given candidate.}
}
\description{
Given the number of ballots, groups, candidates, and votes per ballot, this function simulates an election. Additionally, it generates a resulting matrix (of dimension g x c) representing the unobserved probabilities that a demographic group votes for a given candidate. These probabilities are drawn from a Dirichlet distribution and adjusted by a \code{lambda} value, which represents the heterogeneity of the groups.
}
\note{
Note: For the Dirichlet distribution, an alpha value of 1 yields a distribution that is uniform around the mean. Conversely, alpha values less than 1 tend to produce a sparser probability vector, with many values near zero and one or a few larger values. When alpha is greater than 1, the distribution becomes more concentrated around the mean.
}
\examples{
result <- simulate_elections(num_ballots = 10, 
							 num_candidates = 5, 
							 num_groups = 3, 
							 ballot_voters = c(100, 10))
result$X # Candidate matrix (b x c)
result$W # Group matrix (b x g)
result$real_p # Probability matrix (g x c)
}
\seealso{
[random_samples()]
}
