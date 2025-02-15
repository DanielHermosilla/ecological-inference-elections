\name{update.eim}
\alias{update.eim}
\title{Update an existing eim model with a new EM algorithm computation}
\usage{
\method{update}{eim}(object, ...)
}
\arguments{
\item{...}{New parameters to be passed to \code{compute()}.}
\item{object}{An eim object.}
}
\value{
The updated object with the changes applied.
}
\description{
This function updates an object by performing a new Expectation-Maximization computation with different parameters.
}
\examples{
simulations <- simulate_elections(
  num_ballots = 20,
  num_candidates = 5,
  num_groups = 3,
  ballot_voters = rep(100, 20)
)

model <- eim$new(simulations$X, simulations$W)

predict(model)

update(model, "MVN PDF")

model$method # Returns 'MVN PDF'
}
