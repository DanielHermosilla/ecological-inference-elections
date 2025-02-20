#' Internal function!
#'
#' Validate the 'eim' object inputs
#'
#' @param X A matrix representing candidate votes per ballot box.
#' @param W A matrix representing group votes per ballot box.
#' @return Stops execution if validation fails.
#' @noRd
.validate_eim <- function(X, W, mismatch = FALSE) {
    # Ensure X and W are provided
    if (is.null(X) || is.null(W)) {
        stop("Either provide X and W matrices, or a valid JSON path containing them.")
    }

    # Ensure they are matrices
    X <- as.matrix(X)
    W <- as.matrix(W)

    # Check matching dimensions
    if (!mismatch && nrow(X) != nrow(W)) {
        stop(
            "Mismatch in the number of ballot boxes: 'X' has ", nrow(X),
            " rows, but 'W' has ", nrow(W), " rows."
        )
    }

    # Check minimum column constraints
    if (ncol(X) < 2) {
        stop("Candidate matrix 'X' must have at least 2 columns.")
    }
    if (ncol(W) < 2) {
        stop("Group matrix 'W' must have at least 2 columns.")
    }

    # Check for missing values
    if (any(is.na(X)) || any(is.na(W))) {
        stop("Matrices 'X' and 'W' cannot contain missing values (NA).")
    }

    TRUE
}

#' Validates an input for a given method
#' @noRd
.validate_arg_compute <- function(params, param_name, default_val, int = TRUE) {
    if (!param_name %in% names(params)) {
        return(default_val)
    }

    value <- params[[param_name]]

    # Validate non-negative value
    if (value < 0) {
        stop("Compute:\tThe argument '", param_name, "' isn't positive.")
    }

    # Validate integer
    if (int && value != as.integer(value)) {
        stop("Compute:\tThe argument '", param_name, "' must be a whole number.")
    }

    return(value)
}

#' Internal function!
#'
#' Validate the 'eim' object JSON path
#'
#' @param json_path A path to a JSON file containing `"X"` and `"W"`.
#' @return A list with the `"X"` and `"W"` matrix. Stops execution if validation fails.
#' @noRd
.validate_json_eim <- function(json_path) {
    if (!file.exists(json_path)) {
        stop("The specified JSON file does not exist: ", json_path)
    }

    data <- tryCatch(
        jsonlite::fromJSON(json_path),
        error = function(e) stop("Failed to read JSON file: ", e$message)
    )

    # Validate JSON contents
    if (!all(c("X", "W") %in% names(data))) {
        stop("JSON file must contain the keys 'X' (candidate matrix) and 'W' (group matrix)")
    }

    if (is.null(data$X) || is.null(data$W)) {
        stop("'X' and 'W' cannot be NULL in the JSON file")
    }

    list(
        X = as.matrix(data$X),
        W = as.matrix(data$W)
    )
}

#' Internal function!
#'
#' Randomly create a voting instance by defining an interval
#'
#' @description
#' Given a range of possible \strong{observed} outcomes (such as ballot boxes, number of candidates, etc.),
#' it creates a completely random voting instance, simulating the unobserved results as well.
#'
#' @param ballots_range (integer) A vector of size 2 with the lower and upper bound of ballot boxes.
#'
#' @param candidates_range (integer) A vector of size 2 with the lower and upper bound of candidates to draw.
#'
#' @param demographic_range (integer) A vector of size 2 with the lower and upper bound of demographic groups
#' to draw.
#'
#' @param voting_range (integer) A vector of size 2 with the lower and upper bound of votes per ballot box.
#'
#' @param seed \emph{(numeric(1)} Optional. If provided, it overrides the current global seed. (default: \code{NULL})
#'
#' @return A list with components:
#' \item{X}{A matrix (b x c) with candidate votes per ballot box.}
#' \item{W}{A matrix (b x g) with demographic votes per ballot box.}
#' \item{real_p}{A matrix (g x c) with the estimated \strong{(unobserved)} probabilities that a demographic group votes for a given candidate.}
#' \item{ballots}{The number of ballot boxes that were drawn.}
#' \item{candidates}{The number of candidates that were drawn.}
#' \item{groups}{The number of demographic groups that were drawn.}
#' \item{total_votes}{A vector with the number of total votes per ballot box.}
#'
#' @seealso [simulate_elections()]
#' @examples
#'
#' bal_range <- c(30, 50)
#' can_range <- c(2, 4)
#' group_range <- c(2, 6)
#' voting_range <- c(50, 100)
#' results <- random_samples(bal_range, can_range, group_range, voting_range)
#'
#' # X matrix
#' results$X # A randomly generated matrix of dimension (b x c)
#' ncol(results$X <= can_range[2]) # Always TRUE
#' ncol(results$X >= can_range[1]) # Always TRUE
#' nrow(results$X <= bal_range[2]) # Always TRUE
#' nrow(results$X >= bal_range[1]) # Always TRUE
#'
#' # W matrix
#' results$W # A randomly generated matrix of dimension (b x g)
#' ncol(results$W <= group_range[2]) # Always TRUE
#' ncol(results$W >= group_range[1]) # Always TRUE
#' nrow(results$W <= bal_range[2]) # Always TRUE
#' nrow(results$W >= bal_range[1]) # Always TRUE
#'
#' # Probability matrix
#' results$real_p # A matrix (g x c) that summarizes the unobserved outcomes
#' ncol(results$real_p) == ncol(results$X) # Always TRUE
#' nrow(results$real_p) == ncol(results$W) # Always TRUE
#'
#' @noRd
.random_samples <- function(ballots_range, # Arguments must be vectors of size 2
                            candidates_range,
                            demographic_range,
                            voting_range,
                            seed = NULL) {
    param_list <- list(ballots_range, candidates_range, demographic_range, voting_range)
    if (!(all(sapply(param_list, length) == 2))) {
        stop("The vectors must be of size 2.")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    # Randomly choose a ballot box
    num_ballots <- sample(ballots_range[1]:ballots_range[2], 1)
    # Randomly choose demographic groups
    num_groups <- sample(demographic_range[1]:demographic_range[2], 1)
    # Randomly choose candidates
    num_candidates <- sample(candidates_range[1]:candidates_range[2], 1)
    # Randomly choose the total amount of votes per ballot box
    total_votes <- sample(
        seq.int(voting_range[1], voting_range[2]),
        size = num_ballots,
        replace = TRUE
    )

    choosen_values <- list(
        ballots = num_ballots,
        candidates = num_candidates,
        groups = num_groups,
        total_votes = total_votes
    )

    result <- simulate_elections(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = total_votes,
        seed = seed
    )

    appended_list <- c(result, choosen_values)
    appended_list
}
