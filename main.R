dyn.load("util.so") # TODO: Corregir NAMESPACE, no importa bien util.so
library("Rcpp")
library("infPackage")
# Llamar a la función definida en util.so
print_message <- function() {
  .Call("hello_gsl")
}

# Function for creating the EM-algorithm object
eco <- function(X, W, method, ..., convergence = 0.001, max_iterations = 1000, verbose = TRUE) {
  # Check for valid methods
  valid_methods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
  if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
    stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
  }

  if (!is.matrix(X)) stop("The candidate argument isn't a matrix.")
  if (!is.matrix(W)) stop("The demographical group argument isn't a matrix.")
  # Parse the optional parameters
  args <- list(...)

  # Parse the iterations parameters
  if (!is.null(args$iterations)) {
    max_iterations <- args$iterations
    if (!is.integer(max_iterations) || maxIterations < 0) stop("The maximum iterations provided isn't an valid integer")
  }

  # Parse the convergence parameters
  if (!is.null(args$convergence)) {
    convergence <- args$convergence
    if (!is.numeric(convergence) || convergence < 0) stop("The convergence value isn't a valid number")
  }

  RsetParameters()

  # The custom object parameters
  obj <- list(
    method = ...,
    time = ...,
    threshold = ...,
    iterations = ...,
    X = ...,
    W = ...,
    probabilities = ...,
  )
  class(eco) <- "EcologicalInference"
  obj
}

# Print method
print.eco <- function(object, ...) {

}

# Summary method
summary.eco <- function(object, ...) {
  list(
    method = object$method,
    probability = object$probabilities,
    time = object$time
  )
}
