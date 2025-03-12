#' Chilean 2021 First Round Presidential Election
#'
#' This dataset contains the results of the first round of the 2021 Chilean presidential elections. It provides
#' 9 possible voting options (7 candidates, blank, and null). Each ballot-box is identified by its id (`BALLOT.BOX`) and an electoral circumscription
#' (`ELECTORAL.DISTRICT`). Additionally, it provides demographic information on voters' age range for each ballot box.
#'
#' @docType data
#' @usage data("chile_election_2021")
#'
#' @format A data frame with 46,639 rows and 14 variables:
#' \describe{
#'   \item{\code{ELECTORAL.DISTRICT}}{The electoral circumscription of the ballot box.}
#'   \item{\code{BALLOT.BOX}}{An identifier for the ballot box within the \code{ELECTORAL.DISTRICT}.}
#'   \item{\code{C1}}{The number of votes cast for the candidate *Gabriel Boric*.}
#'   \item{\code{C2}}{The number of votes cast for the candidate *José Antonio Kast*.}
#'   \item{\code{C3}}{The number of votes cast for the candidate *Yasna Provoste*.}
#'   \item{\code{C4}}{The number of votes cast for the candidate *Sebastián Sichel*.}
#'   \item{\code{C5}}{The number of votes cast for the candidate *Eduardo Artés*.}
#'   \item{\code{C6}}{The number of votes cast for the candidate *Marco Enríquez-Ominami*.}
#'   \item{\code{C7}}{The number of votes cast for the candidate *Franco Parisi*.}
#'   \item{\code{BLANK.VOTES}}{The number of blank votes.}
#'   \item{\code{NULL.VOTES}}{The number of null votes.}
#'   \item{\code{X18.19}}{Number of voters aged 18--19.}
#'   \item{\code{X20.29}}{Number of voters aged 20--29.}
#'   \item{\code{X30.39}}{Number of voters aged 30--39.}
#'   \item{\code{X40.49}}{Number of voters aged 40--49.}
#'   \item{\code{X50.59}}{Number of voters aged 50--59.}
#'   \item{\code{X60.69}}{Number of voters aged 60--69.}
#'   \item{\code{X70.79}}{Number of voters aged 70--79.}
#'   \item{\code{X80.}}{Number of voters aged 80 and older.}
#'   \item{\code{MISMATCH}}{Boolean that takes value `TRUE` if the ballot-box has a mismatch between the total number of votes and the total number of voters. If this is not the case, its value is `FALSE`.}
#' }
#'
#' @keywords datasets
#' @source [Chilean Electoral Service (SERVEL)](https://www.servel.cl/)
#' @examples
#' data("chile_election_2021")
#' head(chile_election_2021)
#'
"chile_election_2021"
