## Load database
#' Load the IPD IMGT/HLA database
#' 
#' @param version Either a valid version of the IPD IMGT/HLA database or 
#'   "Latest" to fetch the latest version
#' @return an \code{\link{IpdDb}} object containing the database.
#' @usage loadHlaData(version = "Latest")
#' @examples 
#' ## Load the HLA database
#' hla <- loadHlaData()
#' @export
loadHlaData <- function(version = "Latest") {
  accNumber <- switch(version,
                      "Latest" = "AH63658",
                      "3.32.0"  = "AH63658",
                      NULL)
  if (is.null(accNumber))
    stop(sprintf("IPD IMGT/HLA version %s in not yet included or exist at all", 
                 version))
  .getData(accNumber)
}
#' Load the IPD KIR database
#' 
#' @param version Either a valid version of the IPD KIR database or "Latest" to
#'   fetch the latest version
#' @return an \code{\link{IpdDb}} object containing the database.
#' @usage loadKirData(version = "Latest")
#' @examples
#' ## Load the KIR database
#' kir <- loadKirData()
#' @export
loadKirData <- function(version = "Latest") {
  accNumber <- switch(version,
                      "Latest" = "AH63659",
                      "2.7.1"  = "AH63659",
                      NULL)
  if (is.null(accNumber))
    stop(sprintf("IPD KIR version %s in not yet included or exist at all", 
                 version))
  .getData(accNumber)
}


### Helper ###
.getData <- function(accNumber) {
  hub <- AnnotationHub()
  hub[[accNumber]]
}
