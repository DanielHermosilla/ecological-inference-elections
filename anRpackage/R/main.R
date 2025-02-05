dyn.load("src/infPackage.so")
library(jsonlite)

# Function for creating the EM-algorithm object

lmCustom <- function(X = NULL, W = NULL, jsonPath = NULL, ...) {
  # Check input: if a JSON file is provided, ignore X and W
  if (!is.null(jsonPath) && nzchar(jsonPath)) {
    # Read the JSON file. This assumes the JSON file has keys "X" and "W".
    jsonData <- fromJSON(jsonPath)
    if (!("X" %in% names(jsonData)) || !("W" %in% names(jsonData))) {
      stop("JSON file must contain elements 'X' and 'W'.")
    }
    X <- as.matrix(jsonData$X)
    W <- as.matrix(jsonData$W)
  } else {
    # Ensure that X and W are provided as matrices
    if (is.null(X) || is.null(W)) {
      stop("Either matrices X and W or a JSON file path must be provided.")
    }
    X <- as.matrix(X)
    W <- as.matrix(W)
	RsetParameters(&X, &W)
	params <- list($X = X, $W = W)
	class(params) <- "EcologicalInference"
	return param
  }
}

lmCustom.compute <- function(object, )

eco <- function(X, W, method, ..., verbose = TRUE) {
    # Check for valid methods
    validMethods <- c("Hit and Run", "Exact", "MVN CDF", "MVN PDF", "Multinomial")
    if (!is.character(method) || length(method) != 1 || !(method %in% valid_methods)) {
        stop("Invalid method. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    if (!is.matrix(X)) stop("The candidate argument isn't a matrix.")
    if (!is.matrix(W)) stop("The demographical group argument isn't a matrix.")
    # Parse the optional parameters
    args <- list(...)

    # Parse the iterations parameters
    if (!is.null(args$iterations)) {
        maxIterations <- args$iterations
        if (!is.integer(maxIterations) || maxIterations < 0) stop("The maximum iterations provided isn't an valid integer")
    } else {
        maxIterations <- 1000
    }

    # Parse the convergence parameters
    if (!is.null(args$convergence)) {
        convergence <- args$convergence
        if (!is.numeric(convergence) || convergence < 0) stop("The convergence value isn't a valid number")
    } else {
        convergence <- 0.001
    }

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
