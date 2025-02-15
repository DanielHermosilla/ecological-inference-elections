\name{random_samples}
\alias{random_samples}
\title{Randomly create a voting instance by defining an interval}
\usage{
random_samples(
  ballots_range,
  candidates_range,
  group_range,
  voting_range,
  seed = NULL
)
}
\arguments{
\item{ballots_range}{(integer) A vector of length 2 specifying the lower and upper bounds for the number of ballot boxes.}

\item{candidates_range}{(integer) A vector of length 2 specifying the lower and upper bounds for the number of candidates to draw.}

\item{group_range}{(integer) A vector of length 2 specifying the lower and upper bounds for the number of demographic groups to draw.}

\item{voting_range}{(integer) A vector of length 2 specifying the lower and upper bounds for the number of votes per ballot box.}

\item{seed}{\emph{(numeric(1))} Optional. If provided, it overrides the current global seed.}
}
\value{
A list with components:
\item{X}{A matrix (b x c) containing candidate votes per ballot box.}
\item{W}{A matrix (b x g) containing demographic votes per ballot box.}
\item{real_p}{A matrix (g x c) containing the estimated \strong{(unobserved)} probabilities that a demographic group votes for a given candidate.}
\item{ballots}{The number of ballot boxes that were drawn.}
\item{candidates}{The number of candidates that were drawn.}
\item{groups}{The number of demographic groups that were drawn.}
\item{total_votes}{A vector containing the total number of votes per ballot box.}
}
\description{
Given a range of possible \strong{observed} outcomes (such as the number of ballot boxes, candidates, etc.), it creates a completely random voting instance, simulating the unobserved results as well.
}
\examples{

bal_range <- c(30, 50)
can_range <- c(2, 4)
group_range <- c(2, 6)
voting_range <- c(50, 100)
results <- random_samples(bal_range, can_range, group_range, voting_range)

# X matrix
results$X # A randomly generated matrix of dimension (b x c)
ncol(results$X <= can_range[2]) # Always TRUE
ncol(results$X >= can_range[1]) # Always TRUE
nrow(results$X <= bal_range[2]) # Always TRUE
nrow(results$X >= bal_range[1]) # Always TRUE

# W matrix
results$W # A randomly generated matrix of dimension (b x g)
ncol(results$W <= group_range[2]) # Always TRUE
ncol(results$W >= group_range[1]) # Always TRUE
nrow(results$W <= bal_range[2]) # Always TRUE
nrow(results$W >= bal_range[1]) # Always TRUE

# Probability matrix
results$real_p # A matrix (g x c) that summarizes the unobserved outcomes
ncol(results$real_p) == ncol(results$X) # Always TRUE
nrow(results$real_p) == ncol(results$W) # Always TRUE

}
\seealso{
[simulate_elections()]
}
