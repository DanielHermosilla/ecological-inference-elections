#' Chilean 2021 Presidential Election Votes
#'
#' This dataset contains the results of the first round of the 2021 Chilean presidential elections.
#' It provides 10 possible voting options. Each ballot box is identified by a ballot box number
#' (`MESA`), voting location (`LOCAL`), electoral circumscription
#' (`CIRCUNSCRIPCION.ELECTORAL`), and region (`REGION`).
#'
#' @docType data
#' @usage data("chile_candidate_votes_2021")
#'
#' @format A data frame with 46,639 rows and 14 variables:
#' \describe{
#'   \item{\code{REGION}}{Geographic region of the ballot box.}
#'   \item{\code{CIRCUNSCRIPCION.ELECTORAL}}{The electoral circumscription of the ballot box.}
#'   \item{\code{LOCAL}}{The name of the location where the voting took place (usually a school).}
#'   \item{\code{MESA}}{An identifier for the ballot box within the \code{LOCAL}.}
#'   \item{\code{EDUARDO.ARTES}}{The number of votes cast for the candidate *Eduardo Artés*.}
#'   \item{\code{FRANCO.ALDO.PARISI}}{The number of votes cast for the candidate *Franco Parisi*.}
#'   \item{\code{GABRIEL.BORIC}}{The number of votes cast for the candidate *Gabriel Boric*.}
#'   \item{\code{JOSE.ANTONIO.KAST}}{The number of votes cast for the candidate *José Antonio Kast*.}
#'   \item{\code{MARCO.ENRIQUEZ.OMINAMI}}{The number of votes cast for the candidate *Marco Enríquez-Ominami*.}
#'   \item{\code{SEBASTIA.SICHEL}}{The number of votes cast for the candidate *Sebastián Sichel*.}
#'   \item{\code{VOTOS.EN.BLANCO}}{The number of blank votes.}
#'   \item{\code{VOTOS.NULO}}{The number of null votes.}
#'   \item{\code{YASNA.PROVOSTE}}{The number of votes cast for the candidate *Yasna Provoste*.}
#'   \item{\code{NULO.BLANCO}}{The sum of \code{VOTOS.NULO} and \code{VOTOS.EN.BLANCO}.}
#' }
#'
#' @keywords datasets
#' @source [Chilean Electoral Service (SERVEL)](https://www.servel.cl/)
#' @seealso chile_groups_2021
#' @examples
#' data("chile_candidate_votes_2021")
#' head(chile_candidate_votes_2021)
#'
"chile_candidate_votes_2021"


#' Chilean 2021 Presidential Election Age Distribution
#'
#' This dataset provides age-related voting group information for the first round of the
#' 2021 Chilean presidential election. Each ballot box is characterized by its region
#' (`REGION`), electoral circumscription (`CIRCUNSCRIPCION.ELECTORAL`),
#' voting location (`LOCAL`), and ballot box identifier (`MESA`).
#' The age distribution of voters is given across several ranges (e.g.,
#' `X18.19` for ages 18--19, `X80.` for ages 80 and older).
#'
#' @docType data
#' @usage data("chile_groups_2021")
#'
#' @format A data frame with 46,606 rows and 11 variables:
#' \describe{
#'   \item{\code{REGION}}{Geographic region of the ballot box.}
#'   \item{\code{CIRCUNSCRIPCION.ELECTORAL}}{The electoral circumscription of the ballot box.}
#'   \item{\code{LOCAL}}{The name of the location where the voting took place.}
#'   \item{\code{MESA}}{An identifier for the ballot box within the \code{LOCAL}.}
#'   \item{\code{X18.19}}{Number of voters aged 18--19.}
#'   \item{\code{X20.29}}{Number of voters aged 20--29.}
#'   \item{\code{X30.39}}{Number of voters aged 30--39.}
#'   \item{\code{X40.49}}{Number of voters aged 40--49.}
#'   \item{\code{X50.59}}{Number of voters aged 50--59.}
#'   \item{\code{X60.69}}{Number of voters aged 60--69.}
#'   \item{\code{X70.79}}{Number of voters aged 70--79.}
#'   \item{\code{X80.}}{Number of voters aged 80 and older.}
#' }
#'
#' @keywords datasets
#' @source [Chilean Electoral Service (SERVEL)](https://www.servel.cl/)
#' @seealso chile_candidate_votes_2021
#' @examples
#' data("chile_groups_2021")
#' head(chile_groups_2021)
"chile_groups_2021"
