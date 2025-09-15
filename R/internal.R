#' Internal function!
#' Validates all of the 'run_em' arguments
#' @noRd
.validate_compute <- function(args) {
    # General checks: Vectors aren't accepted
    # if (any(sapply(args, function(x) length(x) > 1))) {
    #    stop("run_em:\tInvalid input: no vector inputs allowed")
    # }

    object_provided <- "object" %in% names(args) || "object1" %in% names(args)
    x_provided <- "X" %in% names(args) || "X1" %in% names(args)
    w_provided <- "W" %in% names(args) || "W1" %in% names(args)
    xw_provided <- x_provided || w_provided
    json_provided <- "json_path" %in% names(args)

    if (x_provided + w_provided == 1) {
        stop("If providing a matrix, 'X' and 'W' must be provided.")
    }

    if (sum(object_provided, xw_provided, json_provided) != 1) {
        stop(
            "You must provide exactly one of the following:\n",
            "(1)\tan `eim` object (initialized with `eim`)\n",
            "(2)\t`X` and `W`\n",
            "(3)\ta `json_path`"
        )
    }

    # Mismatch argument
    if ("allow_mismatch" %in% names(args)) {
        if (!is.logical(args$allow_mismatch)) {
            stop("run_em: Invalid 'allow_mismatch'. It has to be a boolean")
        }
    }

    # Method argument
    valid_methods <- c("mcmc", "exact", "mvn_cdf", "mvn_pdf", "mult", "metropolis")
    if ("method" %in% names(args) &&
        (!is.character(args$method) || length(args$method) != 1 || !(args$method %in% valid_methods))) {
        stop("Invalid 'method'. Must be one of: ", paste(valid_methods, collapse = ", "))
    }

    # Initial prob argument
    valid_p_methods <- c("group_proportional", "proportional", "uniform", "random", "mult", "mcmc", "mvn_cdf", "mvn_pdf", "exact")
    if ("initial_prob" %in% names(args) &&
        (!is.character(args$initial_prob) || length(args$initial_prob) != 1 || !(args$initial_prob %in% valid_p_methods))) {
        stop("Invalid 'initial_prob'. Must be one of: ", paste(valid_p_methods, collapse = ", "))
    }

    if ("maxiter" %in% names(args)) {
        if (!is.numeric(args$maxiter) || as.integer(args$maxiter) != args$maxiter || args$maxiter < 1) { # Infinite are valid, skip this case
            stop("Invalid 'maxiter'. Must be a positive integer.")
        }
    }

    if ("ll_threshold" %in% names(args)) {
        if (!is.infinite(args$ll_threshold) && (!is.numeric(args$ll_threshold) || args$ll_threshold < 0)) { # Infinite are valid, skip this case
            stop("Invalid 'll_threshold'. Must be a positive numeric or infinite value.")
        }
    }

    # Maxtime argument
    if ("maxtime" %in% names(args) &&
        (!is.numeric(args$maxtime) || args$maxtime < 0)) {
        stop("Invalid 'maxtime'. Must be positive.")
    }

    # Stop threshold argument
    if ("param_threshold" %in% names(args)) {
        if (!is.infinite(args$param_threshold) && (!is.numeric(args$param_threshold) || args$param_threshold < 0)) {
            stop("run_em: Invalid 'param_threshold'. Must be a positive numeric or an infinite value.")
        }
        if (args$param_threshold >= 1) {
            warning("Warning: A 'param_threshold' greater or equal than one will always be true after the first iteration.")
        }
        if ("compute_ll" %in% names(args) && !args$compute_ll && is.infinite(args$param_threshold)) {
            stop("You must provide a parameter threshold if 'compute_ll' is FALSE.")
        }
    }

    # Verbose argument
    if ("verbose" %in% names(args) && !is.logical(args$verbose)) {
        stop("run_em: Invalid 'verbose'. It has to be a boolean.")
    }

    # mcmc: mcmc_stepsize argument
    if ("mcmc_stepsize" %in% names(args)) {
        if (!is.numeric(args$mcmc_stepsize) || as.integer(args$mcmc_stepsize) != args$mcmc_stepsize || args$mcmc_stepsize < 0) {
            stop("run_em: Invalid 'mcmc_stepsize'. Must be a positive integer.")
        }
        if (args$mcmc_stepsize < 15) {
            warning("Warning: A small 'mcmc_stepsize' could lead to highly correlated samples.")
        }
    }

    # mcmc: Samples argument
    if ("samples" %in% names(args) &&
        (!is.numeric(args$mcmc_samples) || as.integer(args$mcmc_samples) != args$mcmc_samples || args$mcmc_samples < 0)) {
        stop("run_em: Invalid 'mcmc_samples'. Must be a positive integer.")
    }

    # metropolis: Samples argument
    if ("metropolis_iter" %in% names(args) &&
        (!is.numeric(args$metropolis_iter) || as.integer(args$metropolis_iter) != args$metropolis_iter || args$metropolis_iter < 0)) {
        stop("run_em: Invalid 'metropolis_iter'. Must be a positive integer.")
    }

    # CDF: Mc_method argument
    valid_cdf_methods <- c("genz", "genz2")
    if ("mvncdf_method
" %in% names(args) &&
        (!is.character(args$mvncdf_method) || !args$mvncdf_method
            %in% valid_cdf_methods)) {
        stop("run_em: Invalid 'mvncdf_method
'. Must be one of: ", paste(valid_cdf_methods, collapse = ", "))
    }

    # CDF: Mc_error argument
    if ("mvncdf_error" %in% names(args) &&
        (!is.numeric(args$mvncdf_error) || args$mvncdf_error <= 0)) {
        stop("run_em: Invalid 'mvncdf_error'. Must be a positive number.")
    }

    # CDF: Mc_error argument
    if ("mvncdf_samples" %in% names(args) &&
        (!is.numeric(args$mvncdf_samples) || as.integer(args$mvncdf_samples) != args$mvncdf_samples || args$mvncdf_samples < 0)) {
        stop("run_em: Invalid 'mvncdf_samples'. Must be a positive integer.")
    }

    # Check mismatch
    if ("mismatch" %in% names(args)) {
        if (!is.logical(args$mismatch)) {
            stop("run_em: Invalid 'mismatch'. Must be a boolean value.")
        }
        # if ("method" %in% names(args) && "method" %in% c("exact")) {
        #    stop("run_em: Mismatched results are not supported when using 'exact'.")
        # }
    }

    # Include nboot aswell if bootstrapping is provided
    if ("nboot" %in% names(args) &&
        (!is.numeric(args$nboot) || as.integer(args$nboot) != args$nboot || args$nboot < 0)) {
        stop("Bootstrap: Invalid 'nboot'. Must be a positive integer.")
    }

    valid_sd_methods <- c("maximum", "average")
    if ("sd_statistic" %in% names(args) &&
        (!is.character(args$sd_statistic) || length(args$sd_statistic) != 1 || !(args$sd_statistic %in% valid_sd_methods))) {
        stop("Invalid 'sd_statistic'. Must be one of: ", paste(valid_sd_methods, collapse = ", "))
    }

    if ("sd_threshold" %in% names(args) &&
        (!is.numeric(args$sd_threshold) || args$sd_threshold <= 0)) {
        stop("Invalid 'sd_threshold'. Must be a positive number.")
    }

    if ("alternative" %in% names(args) && !args$alternative %in% c("two.sided", "greater", "less")) {
        stop("Invalid 'alternative'. Must be one of: two.sided, greater, less")
    }

    valid_lp_methods <- c("", "lp", "project_lp")
    if ("adjust_prob_cond_method" %in% names(args) &&
        (!is.character(args$adjust_prob_cond_method) || !(args$adjust_prob_cond_method %in% valid_lp_methods))) {
        stop("Invalid 'adjust_prob_cond_method'. Must be one of: ", paste(valid_sd_methods, collapse = ", "))
    }

    if ("adjust_prob_cond_every" %in% names(args)) {
        if (!is.logical(args$adjust_prob_cond_every)) {
            stop("Invalid 'adjust_prob_cond_every'. Must be a boolean value.")
        }
        if ("adjust_prob_cond_method" %in% names(args) && args$adjust_prob_cond_method == "") {
            warning("You provided 'adjust_prob_cond_every' but not 'adjust_prob_cond_method'. The former will be ignored.")
        }
    }
}


#' Internal function!
#'
#' Validate the 'eim' object inputs
#'
#' @param X A matrix representing candidate votes per ballot box.
#' @param W A matrix representing group votes per ballot box.
#' @return Stops execution if validation fails.
#' @noRd
.validate_eim <- function(X, W) {
    # Ensure X and W are provided
    if (is.null(X) || is.null(W)) {
        stop("Either provide X and W matrices, or a valid JSON path containing them.")
    }

    if (!is.matrix(X) || !is.matrix(W)) {
        stop("'X' and 'W' must be matrices.")
    }

    # Ensure they are matrices
    X <- as.matrix(X)
    W <- as.matrix(W)

    # Check matching dimensions
    if (nrow(X) != nrow(W)) {
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
        # stop("Group matrix 'W' must have at least 2 columns.")
    }

    # Check for missing values
    if (any(is.na(X)) || any(is.na(W))) {
        stop("Matrices 'X' and 'W' cannot contain missing values (NA).")
    }

    TRUE
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
#' @seealso [simulate_election()]
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
    # Randomly choose the group proportions
    group_prop <- rgamma(num_groups, shape = 1, rate = 1)
    group_prop <- group_prop / sum(group_prop)

    choosen_values <- list(
        ballots = num_ballots,
        candidates = num_candidates,
        groups = num_groups,
        total_votes = total_votes
    )

    result <- simulate_election(
        num_ballots = num_ballots,
        num_candidates = num_candidates,
        num_groups = num_groups,
        ballot_voters = total_votes,
        seed = seed,
        group_proportions = group_prop
    )

    appended_list <- c(result, choosen_values)
    appended_list
}
