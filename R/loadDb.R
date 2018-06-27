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
  accNumber <- "NotYetKnown"
  dbfile <- .getData(accNumber)
  .loadDatabase(dbfile)
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
  accNumber <- "notYetKnown"
  dbfile <- .getData(accNumber)
  .loadDatabase(dbfile)
}


### Helper ###
.loadDatabase <- function(dbfile) {
  assert_that(file.exists(dbfile))
  loadDb(dbfile)
}
.getData <- function(accNumber) {
  hub <- AnnotationHub()
  hub[["accNumber"]]
}
