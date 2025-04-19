#' Extracts voting and demographic data matrices for a given electoral district in Chile.
#'
#' This function retrieves the voting results and demographic covariates for a given electoral district from the 2021 Chilean election dataset included in this package. The function returns an [`eim`] object that can be directly used in [`run_em`] or other estimation functions.
#'
#' The function builds the `X` matrix using the number of votes per candidate, and the `W` matrix using the number of voters in each demographic group (e.g., age ranges). Optionally, blank and null votes can be merged into a single additional column (considered as another candidate).
#'
#' Additionally, ballot boxes where the number of votes does not match the number of registered voters (i.e., those where `MISMATCH == TRUE`) can be excluded from the dataset by setting `remove_mismatch = TRUE`.
#'
#' @param elect_district A string indicating the name of the electoral district to extract (e.g., `"NIEBLA"`).
#'
#' @param merge_blank_null Logical indicating whether blank and null votes should be merged into a single column. Defaults to `TRUE`.
#'
#' @param remove_mismatch Logical indicating whether to remove ballot boxes with mismatched vote totals (where `MISMATCH == TRUE`). Defaults to `TRUE`.
#'
#' @return
#' An [`eim`] object with the following attributes:
#' - **X**: A matrix `(b x c)` with the number of votes per candidate (including a column for blank + null votes if `merge_blank_null = TRUE`).
#' - **W**: A matrix `(b x g)` with the number of voters per group (e.g., age ranges) for each ballot box.
#'
#' This object can be passed to functions like [`run_em`] or [`get_agg_proxy`] for estimation and group aggregation.
#'
#' @examples
#' # Load data and create an eim object for the electoral district of "NIEBLA"
#' eim_obj <- get_XW_Chile("NIEBLA", remove_mismatch = FALSE)
#'
#' # Use it to run the EM algorithm
#' result <- run_em(eim_obj)
#'
#' # Use it with group aggregation
#' agg_result <- get_agg_proxy(
#'     object = eim_obj,
#'     sd_threshold = 0.03,
#'     allow_mismatch = TRUE,
#'     seed = 123
#' )
#'
#' agg_result$group_agg
#'
#' @seealso [chile_election_2021]
#' @aliases get_XW_chile()
#' @export
get_XW_chile <- function(elect_district, merge_blank_null = TRUE, remove_mismatch = TRUE) {
    df <- get("chile_election_2021")
    # filter a specific electoral district
    df_ed <- df[df$ELECTORAL.DISTRICT == elect_district, ]
    rownames(df_ed) <- df_ed$BALLOT.BOX
    # remove ballot boxes with mismatch of votes and voters
    if (remove_mismatch) {
        df_ed <- df_ed[df_ed$MISMATCH == FALSE, ]
    }
    # get columns of candidates' votes
    X <- df_ed[df_ed$ELECTORAL.DISTRICT == elect_district, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "BLANK.VOTES", "NULL.VOTES")]
    # merge blank and null votes
    if (merge_blank_null) {
        X$C8 <- X$BLANK.VOTES + X$NULL.VOTES
    }
    X <- X[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")]
    # convert X to matrix
    X <- as.matrix(X[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")])

    # get voters of each demographic group, and convert to matrix
    W <- as.matrix(df_ed[df_ed$ELECTORAL.DISTRICT == elect_district, c("X18.19", "X20.29", "X30.39", "X40.49", "X50.59", "X60.69", "X70.79", "X80.")])

    obj <- eim(X = X, W = W)
    return(obj)
}
