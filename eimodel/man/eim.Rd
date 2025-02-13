\name{eim}
\alias{eim}
\alias{summary.eim}
\alias{predict.eim}
\alias{as.matrix.eim}
\alias{std.eim}
\title{An R6 Class for estimating an probability matrix from an R x C ecological inference problem.}
\usage{
\method{summary}{eim}(object)

\method{predict}{eim}(object, ...)

\method{as.matrix}{eim}(object)

\method{std.eim}{eim}(object, ...)
}
\description{
This class implements an EM algorithm using different methods for approximating the E-step such as \code{"Multinomial",
"Hit and Run", "MVN CDF", "MVN PDF"}, and \code{"Exact"}.
}
\note{
Precomputing would eventually accelerate the Expected Maximization algorithm computation time.

The methods \code{compute, bootstrap} and \code{summary} can also be called as S3 methods with \code{predict()}, \code{std()} and \code{summary()}.
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
simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)

model$precompute("Hit and Run",
          step_size = 1000,
          samples = 5000)
# Changes are made on C internals API

# Example 2: Precompute the Exact method
model$precompute("Exact")

## ------------------------------------------------
## Method `eim$compute`
## ------------------------------------------------

# Example 1: Compute the Expected Maximization with default values

simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute() # Retrieves the object with updated attributes

# Example 2: Compute the Expected Maximization for the Hit and Run method
\donttest{
  model$compute(main_method = "Hit and Run",
          step_size = 3000,
          samples = 1000)
  }

# Example 3: Omit arguments to the Hit and Run method
\dontrun{
  model$compute(main_method = "Hit and Run",
          step_size = 3000)
  # Error; must hand in samples parameter too
  }

# Example 4: Run the MVN CDF with default values
model$compute(main_method = "MVN CDF")

# Example 5: Run a Exact estimation with user defined parameters
model$compute(main_method = "Exact",
        probability_method = "Uniform",
        iterations = 5,
        stopping_threshold = 1e-3)
# Verbose was omitted

## ------------------------------------------------
## Method `eim$bootstrap`
## ------------------------------------------------


simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(X = simulations$X, W = simulations$W)

model$bootstrap(30, 10)

model$standard_deviation # An estimate of the probabilities sd.

## ------------------------------------------------
## Method `eim$print`
## ------------------------------------------------


simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
print(model) # Will print the X and W matrix.

model$compute()
print(model) # Will print the X and W matrix among the EM results.

## ------------------------------------------------
## Method `eim$summary`
## ------------------------------------------------


simulations <- simulate_votation(num_ballots = 5,
        num_candidates = 3,
        num_demographics = 2,
        ballot_voters = rep(100, 5))

model <- eim$new(simulations$X, simulations$W)
a_list <- model$summary()
a_list$method # Not computed yet
a_list$groups # 2
a_list$ballots # 5
names(a_list)
# "candidates" "groups" "ballots" "method" "probabilities" "logLikelihood"

## ------------------------------------------------
## Method `eim$save_results`
## ------------------------------------------------


simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute()
model$save_results("results.rds")  # Save as RDS
model$save_results("results.json") # Save as JSON
model$save_results("results.csv")  # Save final probability as CSV
}
\references{
Check \emph{Thraves, C. and Ubilla, P.: 'Fast Ecological Inference Algorithm for the R×C Case'} for additional information aswell as benchmarking results.
}

\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
    \tabular{ll}{
      \code{X} \emph{(matrix)} \tab An (b x c) matrix with the observed results of the candidate votes (c) on a \cr \tab given
      ballot box (b). Provided manually or loaded from JSON. \cr \cr

      \code{W} \emph{(matrix)} \tab An (b x g) matrix with the observed results of the demographical group \cr \tab votes (g) on
      a given ballot box (b). Provided manually or loaded from JSON. \cr \cr

      \code{method}\emph{(character)} \tab A string indicating the EM method. One of:
      \code{"Multinomial"}, \code{"Hit and Run"}, \cr \tab \code{"MVN CDF"},
      \code{"MVN PDF"} or \code{"Exact"}. \cr \cr

      \code{probability} \emph{(matrix)} \tab An (g x c) matrix that stores the final estimated probabilities of having
      a \cr \tab given group (g) voting for a candidate (c). \cr \cr

      \code{logLikelihood} \emph{(numeric)} \tab A numeric vector that stores the log-likelihood values across
      the total EM \cr \tab iterations. \cr \cr

      \code{total_iterations}\emph{(integer(1))} \tab An integer indicating the total EM iterations performed. \cr \cr

      \code{total_time}\emph{(numeric(1))} \tab The time (in seconds) that the EM algorithm took.
    }
\if{html}{\out{</div>}}
}

\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
    \tabular{ll}{
      \code{samples} \emph{(integer(1))} \tab Active variable to show the \code{Hit and Run} samples if and only if \cr \tab the method is \code{"Hit and Run"}. \cr \cr

      \code{step_size} \emph{(integer(1))} \tab Active variable to show the \code{Hit and Run} step size if and only if \cr \tab the method is \code{"Hit and Run"}. \cr \cr

      \code{multivariate_method} \emph{(character)} \tab Active variable to show the method used to estimate the \cr \tab Multivariate
      Normal CDF, if and only if \code{self$method} is \code{"MVN CDF"}. \cr \cr

      \code{multivariate_error} \emph{(numeric(1))} \tab Active variable to show the error threshold for the Monte Carlo \cr
      \tab simulation of the Multivariate Normal CDF. \cr \cr

      \code{multivariate_iterations} \emph{(numeric(1))} \tab Active variable to show the number of iterations for the Monte \cr \tab Carlo simulation
      of the Multivariate Normal CDF. \cr \cr

      \code{standard_deviation} \emph{(matrix)} \tab Active variable to show an estimate of the probability standard \cr \tab deviation if the
      bootstrapping method has been called.
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
  \code{X} \emph{(matrix)} \tab A matrix (c x b) of observed candidates (c) votes per ballot \cr \tab boxes (b) (optional, required if \code{json_path} is \code{NULL}). \cr \cr
  \code{W} \emph{(matrix)} \tab A matrix (b x g) of observed demographic (g) group votes \cr \tab per ballot boxes (b) (optional, required if \code{json_path} is \code{NULL}). \cr \cr
  \code{json_path} \emph{(character)} \tab A string containing a path to a JSON file with "X" and "W" \cr \tab matrices. (optional, required if \code{X} or \code{W} are \code{NULL}). 
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
Calculates the EM-independent variables for the \code{Hit and Run} and \code{Exact} method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$precompute(method, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
  \code{method} \emph{(character)} \tab The method for precomputing. Options: \code{"Hit and Run"} or \code{"Exact"}. \cr \cr
  \code{...} \emph{(optional)} \tab Additional arguments required by specific methods: \cr \cr
               \tab \strong{"Hit and Run" Method:} \cr \cr
               \tab \code{step_size} \emph{(integer(1))} The step size (M) for the Hit and Run algorithm. \cr \tab Must be positive. \cr \cr
               \tab \code{samples} \emph{(integer(1))} The number of samples (S) to generate.
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining). Updates are made on the C internal memory.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Example 1: Precompute the Hit and Run method
simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)

model$precompute("Hit and Run",
          step_size = 1000,
          samples = 5000)
# Changes are made on C internals API

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
  main_method = "Multinomial",
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
  \code{main_method} \emph{(character)} \tab The method for EM. Options: \code{"Multinomial"}, \code{"MVN CDF"}, \cr \tab \code{"MVN PDF"}, \code{"Hit and Run"}, \code{"Exact"}. \cr \tab (default: \code{"Multinomial"}). \cr \cr
  \code{probability_method} \emph{(character)} \tab Method to obtain the initial probability. Options: \code{"Group"} \cr \tab \code{"proportional"}, \code{"Proportional"}, \code{"Uniform"}. \cr \tab (default: \code{"Group proportional"}). \cr \cr
  \code{iterations} \emph{(integer(1))} \tab The maximum number of EM iterations. \cr \tab (default: 1000). \cr \cr
  \code{stopping_threshold} \emph{(numeric(1))} \tab The minimum difference between successive \cr \tab probabilities for stopping.\cr \tab (default: 0.001). \cr \cr
  \code{verbose} \emph{(logical(1))} \tab Whether to print progress. (default: FALSE). \cr \cr
  \code{...} \emph{(optional)} \tab Additional arguments for specific methods: \cr \cr
               \tab \strong{"Hit and Run"}: \cr \cr
               \tab \code{step_size} \emph{(integer(1))} The step size (M). \cr \cr
               \tab \code{samples} \emph{(integer(1))} The number of samples (S). \cr \cr
               \tab \strong{"MVN CDF"}: \cr \cr
               \tab \code{multivariate_method} \emph{(character)} The integration method \cr \tab (default: \code{"Genz2"}). \cr \cr
               \tab \code{multivariate_error} \emph{(numeric(1))} The integration error \cr \tab threshold (default: \code{1e-6}). \cr \cr
               \tab \code{multivariate_iterations} \emph{(integer(1))} The Monte Carlo \cr \tab iterations (default: \code{5000}).
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{# Example 1: Compute the Expected Maximization with default values

simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute() # Retrieves the object with updated attributes

# Example 2: Compute the Expected Maximization for the Hit and Run method
\donttest{
  model$compute(main_method = "Hit and Run",
          step_size = 3000,
          samples = 1000)
  }

# Example 3: Omit arguments to the Hit and Run method
\dontrun{
  model$compute(main_method = "Hit and Run",
          step_size = 3000)
  # Error; must hand in samples parameter too
  }

# Example 4: Run the MVN CDF with default values
model$compute(main_method = "MVN CDF")

# Example 5: Run a Exact estimation with user defined parameters
model$compute(main_method = "Exact",
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
Runs a bootstrap algorithm to estimate the \strong{standard deviation} of probabilities
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$bootstrap(
  bootstrap_iterations,
  ballot_boxes,
  main_method = "Multinomial",
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
  \code{bootstrap_iterations} \emph{(integer(1))} \tab Number of EM computations to run in the bootstrap. \cr \cr
  \code{ballot_boxes} \emph{(integer(1))} \tab Number of ballot boxes to resample. Must be \cr \tab smaller than the total ballot boxes. \cr \cr
  \code{main_method} \emph{(character)} \tab The EM method. (default \code{"Multinomial"}). \cr \cr
  \code{probability_method} \emph{(character)} \tab The initial probability approach. Defaults to \code{"Group"} \cr \tab \code{"proportional"}. \cr \cr
  \code{iterations} \emph{(integer(1))} \tab Maximum EM iterations. (default: 1000). \cr \cr
  \code{stopping_threshold} \emph{(numeric(1))} \tab Convergence threshold. (default: 0.001). \cr \cr
  \code{verbose} \emph{(logical(1))} \tab Print useful progress (default: FALSE). \cr \cr
  \code{...} \emph{(optional)} \tab Extra arguments for methods computation (e.g., \cr \tab \code{step_size}, \code{samples} for \code{"Hit and Run"}). \cr \tab Refer to them in \code{$compute()}.
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object, storing \code{standard_deviation} in its fields.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(X = simulations$X, W = simulations$W)

model$bootstrap(30, 10)

model$standard_deviation # An estimate of the probabilities' sd.
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-print"></a>}}
\if{latex}{\out{\hypertarget{method-eim-print}{}}}
\subsection{Method \code{print()}}{
According to the state of the algorithm (either computed or not), it prints a message with its most relevant parameters
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$print()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
The object itself (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
print(model) # Will print the X and W matrix.

model$compute()
print(model) # Will print the X, W, and some EM results.
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-summary"></a>}}
\if{latex}{\out{\hypertarget{method-eim-summary}{}}}
\subsection{Method \code{summary()}}{
Returns a list of the most important attributes: \code{method}, #\code{candidates}, #\code{ballots}, #\code{groups}, \code{final probabilities}, and \code{log-likelihood.}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$summary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
\emph{(list)} with keys: \code{"candidates", "groups", "ballots", "method", "probabilities", "logLikelihood".}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_votation(num_ballots = 5,
        num_candidates = 3,
        num_demographics = 2,
        ballot_voters = rep(100, 5))

model <- eim$new(simulations$X, simulations$W)
a_list <- model$summary()
a_list$method
a_list$groups
a_list$ballots
names(a_list)
# "candidates" "groups" "ballots" "method" "probabilities" "logLikelihood"
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-save_results"></a>}}
\if{latex}{\out{\hypertarget{method-eim-save_results}{}}}
\subsection{Method \code{save_results()}}{
Saves the current eim object to a specified file. Accepted formats are:
\itemize{
  \item \strong{RDS} (\code{*.rds}): Binary format preserving object structure. \cr
  \item \strong{JSON} (\code{*.json}): Human-readable JSON. \cr
  \item \strong{CSV} (\code{*.csv}): The final probability matrix only.
}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$save_results(filename)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
  \code{filename} \tab A string with the file name (with extension). The extension decides the format.
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
The modified eim object (for method chaining).
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{
simulations <- simulate_votation(num_ballots = 20,
        num_candidates = 5,
        num_demographics = 3,
        ballot_voters = rep(100, 20))

model <- eim$new(simulations$X, simulations$W)
model$compute()
model$save_results("results.rds")  # Save as RDS
model$save_results("results.json") # Save as JSON
model$save_results("results.csv")  # Save final probability as CSV
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-eim-clone"></a>}}
\if{latex}{\out{\hypertarget{method-eim-clone}{}}}
\subsection{Method \code{clone()}}{
Objects of this class are cloneable via this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{eim$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\tabular{ll}{
  \code{deep} \tab Whether to make a deep clone.
}
\if{html}{\out{</div>}}
}
}
