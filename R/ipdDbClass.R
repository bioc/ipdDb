## define IpdDb class
##

#' The database class for storing allele data from IPD.
#' 
#' This class extends the 
#' \code{\link[AnnotationDbi:AnnDbObj-class]{AnnotationDbi::AnnDbObj-class}}
#' object by higher level methods for sequence and annotation retrieval.
#' blubb
#' 
#' @aliases IpdDb 
#' @rdname IpdDb
#' @seealso 
#'   \code{\link[AnnotationDbi:AnnDbObj-class]{AnnotationDbi::AnnDbObj-class}}
#' @slot getDbVersion() Get the version of the original ipd database
#' @slot getLoci() get all loci from a database, see 
#'   \code{\link{getLoci}}.
#' @slot getReference(alleles) Get the reference sequence for alleles, see 
#'   \code{\link{getReference}}.
#' @slot getStructure(alleles) Get the structures of alleles, see 
#'   \code{\link{getStructure}}.
#' @slot getClosestComplete(allele) Get the closest full-length reference 
#'   sequence of one allele, see \code{\link{getClosestComplete}}.
#' @slot getAlleles(locus) Get all alleles of a locus, see 
#'   \code{\link{getAlleles}}.
#' 
#' @examples 
#' ## load the data 
#' hla <- loadHlaData()
#' ## get all valid keytypes
#' kts <- keytypes(hla)
#' ## get all valid columns
#' cols <- columns(hla)
#' ## get the keys of one keytype
#' kt <- kts[1]
#' keys <- keys(hla, kt)
#' ## Get data of the two first columns for the first 10 keys 
#' cols <- cols[1:10]
#' res <- select(hla, keys, cols, kt)
#'   
.IpdDb <- setRefClass("IpdDb", contains="AnnotationDb", 
                      methods = list(
                        getDbVersion = function(.self) {
                          info <- dbInfo(.self$conn)
                          info[info$name == "DB_VERSION", "value"]
                        },
                        getLoci = function(.self) {
                          getLoci(.self)
                        },
                        getReference = function(.self, allele) {
                          getReference(.self, allele)
                        },
                        getStructure = function(.self, allele) {
                          getStructure(.self, allele)
                        },
                        getClosestComplete = function(.self, allele, 
                                                      locus = NULL) {
                          getClosestComplete(.self, allele, locus)
                        },
                        getAlleles = function(.self, locus) {
                          getAlleles(.self, locus)
                        }
                      ))


## Set generics and methods
## Select interface
#' Columns method queries all columns present in the database
#' @rdname IpdDb 
#' @aliases columns
#' @usage columns(x)
#' @export 
setMethod("columns", "IpdDb", function(x){.cols(x)})
#' The keytypes method returns all columns that can be used as keytypes to 
#'   retrieve keys.
#' @rdname IpdDb
#' @aliases keytypes
#' @usage keytypes(x)
#' @export 
setMethod("keytypes", "IpdDb", function(x){.cols(x)})
#' The keys method gets all keys for a give keytype, i.e. all unique entries of
#'   the keytype.
#' @rdname IpdDb
#' @aliases keys
#' @usage keys(x, keytype, ...)
#' @export 
setMethod("keys", "IpdDb", function(x, keytype, ...){.keys(x, keytype)})
#' The select method queries the database for the values of all give columns for
#'   a set of keys that match a specified keytype.
#' @aliases select
#' @rdname IpdDb 
#' @param x the IpdDb object
#' @param keys The keys for which columns should be selected by select()
#' @param columns The columns to retrieve by select
#' @param keytype The keytype for which the keys are retrieved
#' @param ... Additional arguments. Not used now.
#' @return character vector (keys, columns, keytypes) or a data.frame (select).
#' @usage select(x, keys, columns, keytype, ...)
#' @export 
setMethod("select", "IpdDb", function(x, keys, columns, keytype, ...){
  .select(x, keys, columns, keytype)
})

## Sequence retrieval
setGeneric("getLoci", signature = "x", function(x) {
  standardGeneric("getLoci")
})
#' Get loci
#' 
#' Get all available loci of the KIR or HLA database
#' @aliases getLoci
#' @usage getLoci(x)
#' @param x The database connection; an \code{\link{IpdDb}} object.
#' @return A vector of available loci in the database.
#' @examples
#' ## Load the database 
#' hla <- loadHlaData()
#' ## Get the loci
#' loci <- getLoci(hla)
#' 
#' @export 
setMethod("getLoci", signature = "IpdDb", function(x){
  .getLoci(x)
})

setGeneric("getReference", signature = "x", function(x, allele) {
  standardGeneric("getReference")
})

#' Get reference sequences
#' 
#' Get the reference sequences for alleles.
#' 
#' @aliases getReference
#' @usage getReference(x, allele)
#' @param x The database connection; an \code{\link{IpdDb}} object.
#' @param allele The alleles of interest as a character vector.
#' @return A \code{\link[Biostrings:XStringSet-class]{Biostrings:DNAStringSet}}
#' object with all references.
#' @examples
#' ## Load the database 
#' hla <- loadHlaData()
#' ## Get the loci
#' loci <- getLoci(hla)
#' ## Get alleles of a locus
#' alleles <- getAlleles(hla, loci[1])
#' allelesOfInterest <- alleles[1:10]
#' ## Get the sequences
#' seqs <- getReference(hla, allelesOfInterest)
#' 
#' @export 
setMethod("getReference", signature = "IpdDb", function(x, allele){
  .getReference(x, allele)
})

setGeneric("getStructure", signature = "x", function(x, allele) {
  standardGeneric("getStructure")
})
#' Get gene structures
#' 
#' Get the gene structures for alleles.
#' 
#' @aliases getStructure
#' @usage getStructure(x, allele)
#' @param x The database connection; an \code{\link{IpdDb}} object.
#' @param allele The alleles of interest as a character vector.
#' @return A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges:GRanges}} 
#' object with all gene
#'   structures.
#' @examples
#' ## Load the database 
#' hla <- loadHlaData()
#' ## Get the loci
#' loci <- getLoci(hla)
#' ## Get alleles of a locus
#' alleles <- getAlleles(hla, loci[1])
#' allelesOfInterest <- alleles[1:10]
#' ## Get the structures
#' seqs <- getStructure(hla, allelesOfInterest)
#' 
#' @export 
setMethod("getStructure", signature = "IpdDb", function(x, allele){
  .getStructure(x, allele)
})

setGeneric("getClosestComplete", signature = "x", function(x, allele, 
                                                           locus = NULL) {
  standardGeneric("getClosestComplete")
})
#' Get closest full-length sequence
#' 
#' Get the sequence of the closest allele which for which a full-length 
#'   sequence is available.
#' 
#' @aliases getClosestComplete
#' @usage getClosestComplete(x, allele, locus = NULL)
#' @param x The database connection; an \code{\link{IpdDb}} object.
#' @param allele A single allele as a string.
#' @param locus optional parameter used if the allele identifier is not found.
#' @return A \code{\link[Biostrings:XStringSet-class]{Biostrings:DNAStringSet}}
#' object with the sequence 
#' of the closest full-length allele.
#' @examples
#' ## Load the database 
#' hla <- loadHlaData()
#' ## Get the loci
#' loci <- getLoci(hla)
#' ## Get alleles of a locus
#' alleles <- getAlleles(hla, loci[1])
#' alleleOfInterest <- alleles[1]
#' ## Get the closest complete sequence
#' seqs <- getClosestComplete(hla, alleleOfInterest, loci[1])
#' 
#' @export 
setMethod("getClosestComplete", signature = "IpdDb", function(x, allele, 
                                                              locus = NULL){
  .getClosestComplete(x, allele, locus)
})

setGeneric("getAlleles", signature = "x", function(x, locus) {
  standardGeneric("getAlleles")
})
#' Get alleles
#' 
#' Get all alleles of a given locus.
#' 
#' @aliases getAlleles
#' @usage getAlleles(x, locus)
#' @param x The database connection; an \code{\link{IpdDb}} object.
#' @param locus A single locus as a string.
#' @return A character vector with all alleles of the give locus.
#' @examples
#' ## Load the database 
#' hla <- loadHlaData()
#' ## Get the loci
#' loci <- getLoci(hla)
#' ## Get alleles of a locus
#' alleles <- getAlleles(hla, loci[1])
#' 
#' @export 
setMethod("getAlleles", signature = "IpdDb", function(x, locus){
  .getAlleles(x, locus)
})
