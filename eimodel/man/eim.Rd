\name{eim}
\alias{eim}
\alias{summary.eim}
\alias{predict.eim}
\alias{as.matrix.eim}
\alias{sd.eim}
\title{An R6 Class for estimating an probability matrix from an R x C ecological inference problem.}
\usage{
\method{summary}{eim}(object)

\method{predict}{eim}(object, ...)

\method{as.matrix}{eim}(object)

sd.eim(object, ...)
}
\description{
This class implements an EM algorithm using different methods to approximate the E-step, such as \code{"Multinomial"}, \code{"Hit and Run"}, \code{"MVN CDF"}, \code{"MVN PDF"}, and \code{"Exact"}.
}
\note{
Precomputing can eventually accelerate the computation time of the Expectation-Maximization algorithm.

The \code{$compute, $bootstrap, $print} and \code{$summary} methods can also be called as an S3 method with \code{predict(), sd(), print()} and \code{summary()}.
}
\examples{

## ------------------------------------------------
## Method `eim$new`
## ------------------------------------------------

# Example 1: Create a eim object from a matrix
model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))

# Example 2: Create a eim object from a JSON file
\dontrun{
model2 <- eim$new(json_path = "a/file/path/to/a/file.json")
}

## ------------------------------------------------
## Method `eim$precompute`
## ------------------------------------------------

# Example 1: Precompute the Hit and Run method
simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)

model$precompute("Hit and Run",
					step_size = 1000,
					samples = 5000)
# Changes are made to the C internals API

# Example 2: Precompute the Exact method
model$precompute("Exact")

## ------------------------------------------------
## Method `eim$compute`
## ------------------------------------------------

# Example 1: Compute the Expectation-Maximization 
# with default values

simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute() # Retrieves the object with updated 
# attributes

# Example 2: Compute the Expectation-Maximization 
# for the Hit and Run method
\donttest{
	model$compute(method = "Hit and Run",
					step_size = 3000,
					samples = 1000)
	}

# Example 3: Omit arguments to the Hit and Run 
# method
\dontrun{
	model$compute(method = "Hit and Run",
					step_size = 3000)
	# Error; must supply the samples parameter as well
	}

# Example 4: Run the MVN CDF with default values
model$compute(method = "MVN CDF")

# Example 5: Run an Exact estimation with 
# user-defined parameters
model$compute(method = "Exact",
				probability_method = "Uniform",
 			iterations = 5,
				stopping_threshold = 1e-3)
# Verbose was omitted

## ------------------------------------------------
## Method `eim$bootstrap`
## ------------------------------------------------

simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(X = simulations$X, W = simulations$W)

model$bootstrap(30, 10)

model$sd # An estimate of the probabilities' standard 
# deviation.

## ------------------------------------------------
## Method `eim$print`
## ------------------------------------------------

simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
print(model) # Prints the X and W matrices.

model$compute()
print(model) # Prints the X and W matrices along with 
# the EM results.

## ------------------------------------------------
## Method `eim$summary`
## ------------------------------------------------

simulations <- simulate_elections(num_ballots = 5,
				num_candidates = 3,
				num_groups = 2,
				ballot_voters = rep(100, 5))

model <- eim$new(simulations$X, simulations$W)
a_list <- model$summary()
a_list$method # Not computed yet
a_list$groups # 2
a_list$ballots # 5
names(a_list)
# "candidates" "groups" "ballots" "method" "probabilities" 
# "logLikelihood" "finish_state" "sd"

## ------------------------------------------------
## Method `eim$save_results`
## ------------------------------------------------

simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute()
model$save_results("results.rds")  # Save as RDS
model$save_results("results.json") # Save as JSON
model$save_results("results.csv")  # Save the final 
# probability matrix as CSV
}
\references{
Thraves, C. and Ubilla, P.: 'Fast Ecological Inference Algorithm for the R×C Case'
}
\author{
Charles Thraves, Pablo Ubilla
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\tabular{ll}{
\code{X} \emph{(matrix)} \tab An (b x c) matrix with the observed results of the candidate votes (c) \cr \tab for each ballot box (b). Provided manually or loaded from JSON. \cr \cr

\code{W} \emph{(matrix)} \tab An (b x g) matrix with the observed results of the demographic group \cr \tab votes (g) for each ballot box (b). Provided manually or loaded from JSON. \cr \cr

\code{method} \emph{(character)} \tab A string indicating the EM method. One of: \code{"Multinomial"}, \code{"Hit and Run"}, \cr \tab \code{"MVN CDF"},
\code{"MVN PDF"}, or \code{"Exact"}. \cr \cr

\code{probability} \emph{(matrix)} \tab An (g x c) matrix that will store the final estimated probabilities of a given \cr \tab group (g) voting for a candidate (c). \cr \cr

\code{logLikelihood} \emph{(numeric)} \tab A numeric vector that will store the log-likelihood values over all iterations \cr \tab of the Expectation-Maximization algorithm. \cr \cr

\code{total_iterations} \emph{(integer(1))} \tab An integer indicating the total number of iterations performed by the \cr \tab Expectation-Maximization algorithm. \cr \cr

\code{total_time} \emph{(numeric(1))} \tab The time taken by the EM algorithm. \cr \cr

\code{finish_state} \emph{(character)} \tab The reason for the algorithm's termination. It can be either \cr \tab "Maximum iterations", "Log-likelihood decrease", "Convergence", \cr \tab or "Early exit".
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\tabular{ll}{
\code{samples} \emph{(integer(1))} \tab Active variable to show the \code{Hit and Run} samples if and \cr \tab only if the method is \code{"Hit and Run"}. \cr \cr

\code{step_size} \emph{(integer(1))} \tab Active variable to show the \code{Hit and Run} step size if and \cr \tab only if the method is \code{"Hit and Run"}. \cr \cr

\code{multivariate_method} \emph{(character)} \tab Active variable to show the method used to estimate \cr \tab the Multivariate
Normal CDF, if and only if \code{self$method} is \cr \tab \code{"MVN CDF"}. \cr \cr

\code{multivariate_error} \emph{(numeric(1))} \tab Active variable to show the error threshold for the Monte Carlo \cr \tab simulation of the Multivariate Normal CDF. \cr \cr

\code{multivariate_iterations} \emph{(numeric(1))} \tab Active variable to show the number of iterations for the Monte \cr \tab Carlo simulation of the Multivariate Normal CDF. \cr \cr

\code{sd} \emph{(matrix)} \tab Active variable to show an estimate of the probability standard \cr \tab deviation if the bootstrapping method has been called.
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-eim-new}{\code{eim$new()}}
\item \href{#method-eim-precompute}{\code{eim$precompute()}}
\item \href{#method-eim-compute}{\code{eim$compute()}}
\item \href{#method-eim-bootstrap}{\code{eim$bootstrap()}}
\item \href{#method-eim-print}{\code{eim$print()}}
\item \href{#method-eim-summary}{\code{eim$summary()}}
\item \href{#method-eim-save_results}{\code{eim$save_results()}}
\item \href{#method-eim-clone}{\code{eim$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-new"></a>}}
\if{latex}{\out{\hypertarget{method-eim-new}{}}}
\subsection{Method \code{new()}}{
Creates the object by defining the X and W matrix attributes.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$new(X = NULL, W = NULL, json_path = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{X} \emph{(matrix)} \tab A matrix (c x b) of observed candidate votes per ballot box (b) \cr \tab (optional; required if \code{json_path} is \code{NULL}). \cr \cr

\code{W} \emph{(matrix)} \tab A matrix (b x g) of observed demographic group votes per \cr \tab ballot box (b) (optional; required if \code{json_path} is \code{NULL}). \cr \cr

\code{json_path} \emph{(character)} \tab A string containing a path to a JSON file with "X" and "W" \cr \tab matrices. (optional; required if \code{X} or \code{W} are \code{NULL})
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An initialized eim object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Example 1: Create a eim object from a matrix
model <- eim$new(X = matrix(1:9, 3, 3), W = matrix(1:9, 3, 3))

# Example 2: Create a eim object from a JSON file
\dontrun{
model2 <- eim$new(json_path = "a/file/path/to/a/file.json")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-precompute"></a>}}
\if{latex}{\out{\hypertarget{method-eim-precompute}{}}}
\subsection{Method \code{precompute()}}{
Calculates the EM-independent variables for the \code{Hit and Run} and \code{Exact} methods.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$precompute(method, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{method} \emph{(character)} \tab The method for precomputing. Options: \cr \tab \code{"Hit and Run"} or \code{"Exact"}. \cr \cr

\code{...} \tab Additional arguments required by specific methods: \cr \cr
\strong{"Hit and Run" Method:} \tab  \cr \cr
      \code{step_size} \emph{(integer(1))} \tab The step size (M) for the Hit and Run algorithm. \cr \tab Must be a positive integer. \cr \cr
      \code{samples} \emph{(integer(1))} \tab The number of samples (S) to generate. Must be an \cr \tab integer.
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining). Updates are made on the C internal memory.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Example 1: Precompute the Hit and Run method
simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)

model$precompute("Hit and Run",
					step_size = 1000,
					samples = 5000)
# Changes are made to the C internals API

# Example 2: Precompute the Exact method
model$precompute("Exact")
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-compute"></a>}}
\if{latex}{\out{\hypertarget{method-eim-compute}{}}}
\subsection{Method \code{compute()}}{
Executes the Expectation-Maximization (EM) algorithm based on the selected method. Additional parameters may be required depending on the method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$compute(
  method = "Multinomial",
  probability_method = "Group proportional",
  iterations = 1000,
  maximum_time = 1440,
  stopping_threshold = 0.001,
  verbose = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{method} \emph{(character)} \tab The method for estimating the \cr \tab Expectation-Maximization (EM) algorithm. \cr \tab Options:
\code{"Multinomial"}, \code{"MVN CDF"}, \cr \tab \code{"MVN PDF"}, \code{"Hit and Run"}, \cr \tab and \code{"Exact"} \cr \tab (default: \code{"Multinomial"}). \cr \cr

\code{probability_method} \emph{(character)} \tab The method for obtaining the initial \cr \tab probability. Options: \code{"Group proportional"},
\cr \tab \code{"Proportional"}, or \code{"Uniform"}. \cr \tab (default: \code{"Group proportional"}). \cr \cr

\code{iterations} \emph{(integer(1))} \tab The maximum number of iterations to \cr \tab perform on the EM algorithm. \cr \tab (default: \code{1000}). \cr \cr

\code{maximum_time} \emph{(integer(1))} \tab The maximum time (in minutes) that the \cr \tab algorithm will  iterate. (default: \code{1440}). \cr \cr

\code{stopping_threshold} \emph{(numeric(1))} \tab The minimum difference between successive \cr \tab probabilities for stopping the iterations. \cr \tab 
(default: \code{0.001}). \cr \cr

\code{verbose} \emph{(boolean(1))} \tab A boolean indicating whether to print \cr \tab informative messages while iterating. \cr \tab (default: \code{FALSE}). \cr \cr

\code{...} \tab Additional arguments required by specific \cr \tab methods: \cr \cr

\strong{"Hit and Run" Method:} \cr \cr
  \code{step_size} \emph{(integer(1))}: \tab The step size (M) for the Hit and Run \cr \tab algorithm. \cr \cr
  \code{samples} \emph{(integer(1))}: \tab The number of samples (S) to generate. \cr \cr
\strong{"MVN CDF" Method:}  \cr \cr
  \code{multivariate_method} \emph{(character)}: \tab The integration method. Can be chosen \cr \tab between \code{"Genz"} and \code{"Genz2"}. Default is \cr \tab \code{"Genz2"}. \cr \cr
  \code{multivariate_error} \emph{(numeric(1))}: \tab The integration error threshold. Default is \cr \tab \code{1e-6}. \cr \cr
  \code{multivariate_iterations} \emph{(integer(1))}: \tab The number of Monte Carlo iterations. \cr \tab Default is \code{5000}.
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Example 1: Compute the Expectation-Maximization with default values

simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute() # Retrieves the object with updated attributes

# Example 2: Compute the Expectation-Maximization for the Hit and Run method
\donttest{
	model$compute(method = "Hit and Run",
					step_size = 3000,
					samples = 1000)
	}

# Example 3: Omit arguments to the Hit and Run method
\dontrun{
	model$compute(method = "Hit and Run",
					step_size = 3000)
	# Error; must supply the samples parameter as well
	}

# Example 4: Run the MVN CDF with default values
model$compute(method = "MVN CDF")

# Example 5: Run an Exact estimation with user-defined parameters
model$compute(method = "Exact",
				probability_method = "Uniform",
 			iterations = 5,
				stopping_threshold = 1e-3)
# Verbose was omitted
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-bootstrap"></a>}}
\if{latex}{\out{\hypertarget{method-eim-bootstrap}{}}}
\subsection{Method \code{bootstrap()}}{
Runs a bootstrap algorithm to estimate the \strong{standard deviation} of probabilities.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$bootstrap(
  bootstrap_iterations,
  ballot_boxes,
  method = "Multinomial",
  probability_method = "Group proportional",
  iterations = 1000,
  stopping_threshold = 0.001,
  verbose = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{bootstrap_iterations} \emph{(integer(1))} \tab The number of EM computations. \cr \cr

\code{ballot_boxes} \emph{(integer(1))} \tab The number of ballot boxes to use \cr \tab as a sample. It must be \cr \tab strictly less than the total \cr \tab number of ballot boxes. \cr \cr

\emph{Arguments for calling the} \code{$compute()} \cr \emph{method}: \cr \cr

\code{method} \emph{(character)} \tab The method for estimating the \cr \tab Expectation-Maximization (EM) algorithm. \cr \tab Options:
\code{"Multinomial"}, \code{"MVN CDF"}, \cr \tab \code{"MVN PDF"}, \code{"Hit and Run"}, \cr \tab and \code{"Exact"} \cr \tab (default: \code{"Multinomial"}). \cr \cr

\code{probability_method} \emph{(character)} \tab The method for obtaining the initial \cr \tab probability. Options: \code{"Group proportional"},
\cr \tab \code{"Proportional"}, or \code{"Uniform"}. \cr \tab (default: \code{"Group proportional"}). \cr \cr

\code{iterations} \emph{(integer(1))} \tab The maximum number of iterations to \cr \tab perform on the EM algorithm. \cr \tab (default: \code{1000}). \cr \cr

\code{maximum_time} \emph{(integer(1))} \tab The maximum time (in minutes) that the \cr \tab algorithm will  iterate. (default: \code{1440}). \cr \cr

\code{stopping_threshold} \emph{(numeric(1))} \tab The minimum difference between successive \cr \tab probabilities for stopping the iterations. \cr \tab 
(default: \code{0.001}). \cr \cr

\code{verbose} \emph{(boolean(1))} \tab A boolean indicating whether to print \cr \tab informative messages while iterating. \cr \tab (default: \code{FALSE}). \cr \cr

\code{...} \tab Additional arguments required by specific \cr \tab methods: \cr \cr

\strong{"Hit and Run" Method:} \cr \cr
  \code{step_size} \emph{(integer(1))}: \tab The step size (M) for the Hit and Run \cr \tab algorithm. \cr \cr
  \code{samples} \emph{(integer(1))}: \tab The number of samples (S) to generate. \cr \cr
\strong{"MVN CDF" Method:}  \cr \cr
  \code{multivariate_method} \emph{(character)}: \tab The integration method. Can be chosen \cr \tab between \code{"Genz"} and \code{"Genz2"}. Default is \cr \tab \code{"Genz2"}. \cr \cr
  \code{multivariate_error} \emph{(numeric(1))}: \tab The integration error threshold. Default is \cr \tab \code{1e-6}. \cr \cr
  \code{multivariate_iterations} \emph{(integer(1))}: \tab The number of Monte Carlo iterations. \cr \tab Default is \code{5000}.

}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object, under the field \code{$sd}.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)

model$bootstrap(30, 10)

model$sd # An estimate of the standard deviation of the probabilities.
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-print"></a>}}
\if{latex}{\out{\hypertarget{method-eim-print}{}}}
\subsection{Method \code{print()}}{
Depending on the state of the algorithm (computed or not), it prints a message with its most relevant parameters.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$print()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The object itself (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
print(model) # Prints the X and W matrices.

model$compute()
print(model) # Prints the X and W matrices along with the EM results.
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-summary"></a>}}
\if{latex}{\out{\hypertarget{method-eim-summary}{}}}
\subsection{Method \code{summary()}}{
Shows, in the form of a list, a selection of the most important attributes. It retrieves the method, number of candidates, ballots, and groups, as well as the principal results of the EM algorithm.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$summary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\emph{(list)} A list with the \code{method, candidates, ballots, groups, probabilities, log-likelihood, finish_state} and \code{standard deviation.}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_elections(num_ballots = 5,
				num_candidates = 3,
				num_groups = 2,
				ballot_voters = rep(100, 5))

model <- eim$new(simulations$X, simulations$W)
a_list <- model$summary()
a_list$method # Not computed yet
a_list$groups # 2
a_list$ballots # 5
names(a_list)
# "candidates" "groups" "ballots" "method" "probabilities"
# "logLikelihood" "finish_state" "sd"
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-save_results"></a>}}
\if{latex}{\out{\hypertarget{method-eim-save_results}{}}}
\subsection{Method \code{save_results()}}{
Saves the current eim object to a specified file. The results can be saved in:
  \itemize{
    \item \strong{RDS (Binary format):} Preserves object structure for future use in R. \cr \cr
    \item \strong{JSON:} Saves model data in a human-readable format. \cr \cr
    \item \strong{CSV:} Saves the probability matrix in a tabular format.
  }
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$save_results(filename)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{filename} \emph{(character)} \tab The file name where the results should be saved, \cr \tab including its extension.
The file \cr \tab extension determines the format: \cr \cr
\var{*.rds} \tab Saves as \file{RDS} (default, binary format). \cr \cr
\var{*.json} \tab Saves as \file{JSON} (readable and shareable). \cr \cr
\var{*.csv} \tab Saves the \emph{final probability matrix} as \file{CSV}.
}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_elections(num_ballots = 20,
				num_candidates = 5,
				num_groups = 3,
				ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute()
model$save_results("results.rds")  # Save as RDS
model$save_results("results.json") # Save as JSON
model$save_results("results.csv")  # Save the final probability matrix as CSV
}
\if{html}{\out{</div>}}

}

\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-clone"></a>}}
\if{latex}{\out{\hypertarget{method-eim-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class can be cloned using this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
\code{deep} \emph{(bool(1))} \tab Whether to make a deep clone.
}
\if{html}{\out{</div>}}
}
}
}
