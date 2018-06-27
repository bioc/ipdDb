## Load database
#' Load the IPD IMGT/HLA database
#' 
#' @return an \code{\link{IpdDb}} object containing the database.
#' @usage loadHlaData()
#' @examples 
#' ## Load the HLA database
#' hla <- loadHlaData()
#' @export
loadHlaData <- function() {
  accNumber <- "AH63658"
  .getData(accNumber)
}

#' Load the IPD KIR database
#' 
#' @return an \code{\link{IpdDb}} object containing the database.
#' @usage loadKirData()
#' @examples
#' ## Load the KIR database
#' kir <- loadKirData()
#' @export
loadKirData <- function() {
  accNumber <- "AH63659"
  .getData(accNumber)
}


### Helper ###
.getData <- function(accNumber) {
  hub <- AnnotationHub()
  hub[[accNumber]]
}
