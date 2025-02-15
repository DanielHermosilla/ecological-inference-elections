\name{sum.eim}
\alias{sum.eim}
\title{Returns the total number of voters in the system.}
\usage{
\method{sum}{eim}(object)
}
\arguments{
\item{object}{An eim object.}
}
\value{
The total number of voters.
}
\description{
Given an initialized object, it returns the total number of voters. It is equivalent to running \code{sum(object$W)} or \code{sum(object$X)}.
}
\examples{
simulations <- simulate_elections(
  num_ballots = 20,
  num_candidates = 5,
  num_groups = 3,
  ballot_voters = rep(100, 20)
)

model <- eim$new(simulations$X, simulations$W)
sum(model) # 2000

}
