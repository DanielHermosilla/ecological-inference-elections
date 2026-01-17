#' Extracts voting and demographic data matrices for a given electoral district in Chile.
#'
#' This function retrieves the voting results and demographic covariates for a given electoral district from the 2021 Chilean election dataset included in this package. The function returns an [`eim`] object that can be directly used in [`run_em`] or other estimation functions.
#'
#' The function builds the `X` matrix using the number of votes per candidate, and the `W` matrix using the number of voters in each demographic group (e.g., age ranges). Optionally, blank and null votes can be merged into a single additional column (considered as another candidate).
#'
#' Additionally, ballot boxes where the number of votes does not match the number of registered voters (i.e., those where `MISMATCH == TRUE`) can be excluded from the dataset by setting `remove_mismatch = TRUE`.
#'
#' @param elect_district A string indicating the name of the electoral district to extract (e.g., "NIEBLA"). See **Note**.
#'
#' @param region A string indicating the name of the region to extract (e.g, "DE TARAPACA"). See **Note**.
#'
#' @param merge_blank_null Logical indicating whether blank and null votes should be merged into a single column. Defaults to `TRUE`.
#'
#' @param remove_mismatch Logical indicating whether to remove ballot boxes with mismatched vote totals (where `MISMATCH == TRUE`). Defaults to `FALSE`.
#'
#' @param use_sex Logical indicating whether to use the sex from the voters instead of the age ranges. Defaults to `FALSE`.
#'
#' @return
#' An [`eim`] object with the following attributes:
#' - **X**: A matrix `(b x c)` with the number of votes per candidate (including a column for blank + null votes if `merge_blank_null = TRUE`).
#' - **W**: A matrix `(b x g)` with the number of voters per group (e.g., age ranges) for each ballot box.
#'
#' This object can be passed to functions like [`run_em`] or [`get_agg_proxy`] for estimation and group aggregation. See **Example**.
#'
#' @note
#' This function returns an `eim` object with no covariates, i.e., `V=NULL`.
#'
#' Only one parameter is accepted among `elect_district` and `region`. If either both parameters are given, it will return an error. If neither of these two inputs is supplied, it will return an `eim` object with an aggregation corresponding to the whole dataset. To see all electoral districts and regions names, see the function [chile_election_2021].
#'
#' @examples
#' # Load data and create an eim object for the electoral district of "NIEBLA"
#' eim_obj <- get_eim_chile(elect_district = "NIEBLA", remove_mismatch = FALSE)
#'
#' # Use it to run the EM algorithm
#' result <- run_em(eim_obj, allow_mismatch = TRUE)
#'
#' # Use it with group aggregation
#' agg_result <- get_agg_proxy(
#'     object = eim_obj,
#'     sd_threshold = 0.05,
#'     allow_mismatch = TRUE,
#'     seed = 123
#' )
#'
#' agg_result$group_agg
#'
#' @seealso [chile_election_2021]
#' @aliases get_eim_chile()
#' @export
get_eim_chile <- function(elect_district = NULL,
                          region = NULL,
                          merge_blank_null = TRUE,
                          remove_mismatch = FALSE,
                          use_sex = FALSE) {
    df <- get("chile_election_2021")

    # Apply filtering only if exactly one of region or elect_district is provided
    if (!is.null(region) && is.null(elect_district)) {
        df_ed <- df[df$REGION == toupper(region), ]
        rownames(df_ed) <- paste(df_ed$ELECTORAL.DISTRICT, df_ed$BALLOT.BOX, sep = " - ")
    } else if (!is.null(elect_district) && is.null(region)) {
        df_ed <- df[df$ELECTORAL.DISTRICT == toupper(elect_district), ]
        rownames(df_ed) <- df_ed$BALLOT.BOX
    } else if (is.null(elect_district) && is.null(region)) {
        # If both are provided or both are NULL, use full data
        # Generate unique composite key for rownames
        df_ed <- df
        df_ed$row_id <- paste(df_ed$REGION, df_ed$ELECTORAL.DISTRICT, df_ed$BALLOT.BOX, sep = " - ")
        # Remove duplicates based on this composite key
        df_ed <- df_ed[!duplicated(df_ed$row_id), ]
        # Assign rownames
        rownames(df_ed) <- df_ed$row_id
    } else {
        stop("You cannot provide an electoral district and a region simultaneously.")
    }

    # Remove mismatches if applicable
    if (remove_mismatch && "MISMATCH" %in% names(df_ed)) {
        df_ed <- df_ed[df_ed$MISMATCH == FALSE, ]
    }

    if (nrow(df_ed) == 0) {
        # Handle the empty case
        stop("No rows matched the filter.")
    }

    # Extract candidate votes
    X <- df_ed[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "BLANK.VOTES", "NULL.VOTES")]

    # Merge blank and null votes if requested
    if (merge_blank_null) {
        X$C8 <- X$BLANK.VOTES + X$NULL.VOTES
    }

    # Keep only candidate columns including merged C8
    X <- X[, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")]
    X <- as.matrix(X)

    # Extract demographic groups
    if (!use_sex) {
        W <- as.matrix(df_ed[, c("X18.19", "X20.29", "X30.39", "X40.49", "X50.59", "X60.69", "X70.79", "X80.")])
    } else {
        W <- as.matrix(df_ed[, c("M", "F")])
    }

    return(eim(X = X, W = W))
}
