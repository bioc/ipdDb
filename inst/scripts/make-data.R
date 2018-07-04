## This script holds condensed information on how the hla and kir databases
## for the ipdDb package are created. This requires four steps:
## 1. download the data from github:
##   Hla: https://github.com/ANHIG/IMGTHLA We use the hla/hla.xml.zip in that
##     repo
##   Kir: https://github.com/ANHIG/IPDKIR We use the KIR.dat, as there is not
##    yet a reliable xml file
## 2. Read the file and store it in a data.frame. Extract the sequences as
##   DNAStringSet and the structures as GRanges objects.
## 3. Find the closest full-length allele. All alleles of one locus
##   are put into groups of "feature compositions". We have a group for each
##   kind of partial allele, e.g. all alleles were only exon 2 and 3 are known
##   are put together. All full-length alleles are also grouped and known part
##   of the "feature composition" is extracted, e.g. only exon 2 and 3.
##   A distance matrix is build from all sequences in the group. The full-length
##   allele with the minimal distance is considered the closest full-length
##   allele.
## 4. The dataframes are put into an SQLite database based on the AnnotationDbi
##   vignette.
##




library(BSgenome)
library(dplyr)
library(foreach)
library(Biostrings)
library(BiocParallel)
### Get IPD Database ####
###
###
#' Clone the IMGTHLA repo on github
#'
#' @note Due to limitations in the \pkg{git2r} this will not work
#' behind a proxy.
#' @param family one of "kir" or "hla". If none is provided it will fetch HLA.
#' @param url url to fetch from
#' @return The path to the local repo
fetchIPD <- function(family="hla", url=NULL) {
  hla_url <- "https://github.com/ANHIG/IMGTHLA"
  kir_url <- "https://github.com/ANHIG/IPDKIR"
  if (is.null(url))
    url <- ifelse(tolower(family) == "kir", kir_url, hla_url)
  assertthat::assert_that(RCurl::url.exists(url))
  assertthat::assert_that(all(
    requireNamespace("git2r", quietly = TRUE)
  ))
  tmp <- tempfile()
  # Currently not working behind proxy
  git2r::clone(url, tmp)
  tmp
}

#' Fetch the IMGT/HLA hla.xml file
#'
#' @param dbpath the path to the xml file
#' @return An object of class (S3) \code{XMLInternalDocument}.
readHlaXml <- function(dbpath) {
  dbfile <- normalizePath(file.path(dbpath, "xml", "hla.xml.zip"),
                          mustWork = TRUE)
  tfile <- unzip(zipfile = dbfile, exdir = dbpath)
  tfile <- tfile[grepl(pattern = "^hla.xml", x = basename(tfile))]
  doc <- xml2::read_xml(tfile)
  unlink(dbpath, force = TRUE)
  doc
}

#' Fetch the KIR.da file in EMBL format ##
#'
#' @param dbpath the path to the KIR.dat file
#' @return A text representation of the KIR.dat file
readKirDat <- function(dbpath) {
  dbfile <- normalizePath(file.path(dbpath, "KIR.dat"), mustWork = TRUE)
  kirDat <- readLines(file(dbfile))
  unlink(dbpath, force = TRUE)
  kirDat
}

#' Get clostest complete neighbour allele
#'
#' @importFrom Biostrings getSeq
#' @importFrom Biostrings DNAStringSet
#' @import BSgenome

getClosest <- function(x){
  message("Infer the closest complete neighbours")
  # assertthat::assert_that(is("ipdDF", x))
  cmpl <- x$annot[x$annot$complete,]$GID
  ncmpl <- x$annot[!x$annot$complete,]$GID
  ann_ncmpl <- x$features[x$features$GID %in% ncmpl,]

  seqs <- DNAStringSet(x$seq$seq)
  names(seqs) <- x$seq$GID

  ## We need to know which partial genes are there
  featureGroups <- .getFeatureGroups(ann_ncmpl)
  message(sprintf("Found %s distinct feature compositions",
                  length(unique(featureGroups$feat))))

  ## Iterate over all possible partials and calc the corresponding dists
  allNeighbours <- lapply(unique(featureGroups$feat), function(featurePattern) {
    # featurePattern <- unique(featureGroups$feat)[1]
    # loi
    message(sprintf("Process feature composition %s", featurePattern))
    partialSeqs <- .getPartialSeqs(x$features, cmpl, featurePattern, seqs)
    neighbours <- lapply(unique(x$annot$locus), function(loi) {
      message(sprintf("Get closest neighbours for Locus %s", loi))
      .getDistMatrices(x$annot,
                       loi,
                       featureGroups,
                       featurePattern,
                       partialSeqs,
                       seqs)
    })
    do.call(rbind, neighbours)
  })
  do.call(rbind, allNeighbours)
}


## Get the distance matrix for each locus and report the nearest neighbour
.getDistMatrices <- function(annot, loi, featureGroups, featurePattern,
                             partialSeqs, seqs){
  allele_locus <- filter(annot, locus == loi)$GID
  aoi <- allele_locus[
    allele_locus %in%
      filter(featureGroups, feat == featurePattern)$GID]
  caoi <- names(partialSeqs)[names(partialSeqs) %in% allele_locus]
  allPartialSeqs <- c(partialSeqs[caoi],
                      seqs[names(seqs) %in% aoi])
  allPartialAln <- DECIPHER::AlignSeqs(allPartialSeqs,
                                       processors = getOption("ipd.ncores"))
  allPartialDm  <- DECIPHER::DistanceMatrix(
    allPartialAln, processors = getOption("ipd.ncores"))
  do.call(rbind, lapply(aoi, function(x) {
    data_frame(GID = x,
               closest_complete = names(which.min(allPartialDm[x,caoi])))
  }))
}

## Get the partial seqs composed of the desired features
.getPartialSeqs <- function(features, cmpl, featPattern, seqs){
  ranges <- GenomicRanges::GRangesList(BiocParallel::bplapply(cmpl, function(s)
    .extractRanges(s, featPattern, features)))
  partialSeqs <- getSeq(seqs, ranges)
  partialSeqs <- DNAStringSet(lapply(partialSeqs, unlist))
  names(partialSeqs) <- cmpl
  partialSeqs
}

.extractRanges <- function(allele, featPattern, features){
  featPattern <- strsplit(featPattern, ":")[[1]]
  alleleFeatures <- filter(features,
                           GID == allele,
                           feat_name %in% featPattern)
  GenomicRanges::GRanges(seqnames = alleleFeatures$GID,
                         ranges = IRanges::IRanges(start = alleleFeatures$start,
                                                   end = alleleFeatures$end,
                                                   names = alleleFeatures$GID))
}
## extract all possible and valid intron/exon combinations
.getFeatureGroups <- function(x){
  total_feats <- x %>%
    group_by(GID) %>%
    summarize(feat = paste(feat_name, collapse = ":"))
}

### Create HLA data without hlatools
###

#' Parse all HLA alleles for a locus from hla.xml
#'
#' @param doc \file{hla.xml} as an \code{xml_document} object.
#' @param locusname One of [\code{HLA-A, HLA-B, HLA-C, HLA-DPA1, HLA-DPB1,
#' HLA-DQA1, HLA-DQB1, HLA-DRB}]
#' @param ncores The number of compute cores to use.
#'
#' @return A \code{\linkS4class{HLAAllele}} object.
#'
parseHlaAlleles <- function(doc, locusnames = NULL,
                            ncores = getOption("ipd.ncores")) {
  if (is.null(locusnames))
    locusnames <- .VALID_HLA_LOCI()
  locusnames <- sapply(locusnames, .matchHlaLocus)
  if (is.null(ncores))
    ncores <- parallel::detectCores()/2
  assertthat::is.count(ncores)

  ns <- xml2::xml_ns(doc)
  xpaths1 <- sapply(locusnames, function(locusname)
    paste0("/d1:alleles/d1:allele/d1:locus[@locusname='",
           locusname,
           "']/parent::node()"))
  nodess1 <- lapply(xpaths1, function(xpath1) xml2::xml_find_all(doc,
                                                                 xpath1, ns))
  xpath2 <- paste0(".//d1:releaseversions[not(starts-with(@releasestatus,
                   'Allele Deleted'))]/parent::node()")
  nodes2 <- lapply(nodess1, function(nodes1) xml2::xml_find_all(nodes1,
                                                                xpath2, ns))

  releaseVersion <- xml2::xml_find_chr(
    nodes2[[1]][1], "string(./d1:releaseversions/@currentrelease)", ns)

  seq <- do.call(rbind, lapply(nodes2, .parseSequence))
  annot <- do.call(rbind,
                  lapply(1:length(nodes2),
                         function(lc)
                           .parseAnnotations(nodes2[[lc]],
                                             names(nodes2[lc]))))
  annot <- bind_cols(GID = annot$allele_name, annot)
  features <- do.call(rbind,
                      lapply(nodes2, function(x)
                        .parseFeatures(x, ncores = ncores)))

  list(
    seq      = seq,
    annot    = annot,
    features = features,
    version  = releaseVersion)
}

############# Helper to extract the data ###########
##### Extract sequence #####
.parseSequence = function(nodes) {
  ns     <- xml2::xml_ns(nodes)
  xpath  <- './d1:sequence/d1:nucsequence'
  data_frame(
    GID = xml2::xml_attr(nodes, "name"),
    seq     = xml2::xml_text(xml2::xml_find_all(nodes, xpath, ns)),
  )
}

## Extract annotations
.parseAnnotations = function(nodes, locus) {
  ns      <- xml2::xml_ns(nodes)
  res <- data_frame(
    ##
    ## Allele designation
    ##
    locus         = locus,
    allele_name   = xml2::xml_attr(nodes, "name"),
    allele_id     = xml2::xml_attr(nodes, "id"),
    g_group       = xml2::xml_find_chr(
      nodes, "string(./d1:hla_g_group/@status)", ns),
    p_group       = xml2::xml_find_chr(
      nodes, "string(./d1:hla_p_group/@status)", ns),
    date_assigned = xml2::xml_attr(nodes, "dateassigned"),
    ##
    ## Release
    ##
    first_released = xml2::xml_find_chr(
      nodes, "string(./d1:releaseversions/@firstreleased)", ns),
    last_updated   = xml2::xml_find_chr(
      nodes, "string(./d1:releaseversions/@lastupdated)", ns),
    release_status = xml2::xml_find_chr(
      nodes, "string(./d1:releaseversions/@releasestatus)", ns),
    confirmed      = xml2::xml_find_lgl(
      nodes, "string(./d1:releaseversions/@confirmed)=\"Confirmed\"", ns),
    ##
    ## CWD status and Completeness (we consider as complete alleles for which
    ## both UTRs are present)
    ##
    cwd_status    = xml2::xml_find_chr(
      nodes, "string(./d1:cwd_catalogue/@cwd_status)", ns),
    complete      = xml2::xml_find_lgl(
      nodes, "count(./d1:sequence/d1:feature[@featuretype=\"UTR\"])=2", ns),
    ##
  )
  res
}

##### Extract gene features ######
.parseFeatures = function(nodes, ncores) {
  doParallel::registerDoParallel(cores = ncores)
  ns       <- xml2::xml_ns(nodes)
  nodeset  <- xml2::xml_find_all(nodes, "./d1:sequence", ns)
  xpath    <- "./d1:feature[not(@featuretype=\"Protein\")]"
  seqnames <- xml2::xml_attr(nodes, "name")
  res <- foreach(node = nodeset, seqname = seqnames, .combine = rbind) %dopar% {
    data_frame(
      GID = seqname,
      start    = as.integer(xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/d1:SequenceCoordinates/@start"), ns))),
      end      = as.integer(xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/d1:SequenceCoordinates/@end"), ns))),
      feat_name = xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/@name"), ns)),
      id       = xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/@id"), ns)),
      feat_order    = as.integer(xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/@order"), ns))),
      type     = xml2::xml_text(xml2::xml_find_all(
        node, paste0(xpath, "/@featuretype"), ns)),
      frame    = vapply(xml2::xml_find_all(node, xpath, ns), function(node) {
        frame <- xml2::xml_find_chr(
          node, "string(./d1:cDNACoordinates/@readingframe)", ns)
        as.integer(ifelse(frame == "", 0, frame))
      }, FUN.VALUE = 0L)
      )
  }
  res
}





#' Parse all KIR alleles for a locus from KIR.dat
#'
#' @param doc \file{KIR.dat} as an line vector object.
#' @return A list of R6 class EmblEntry objects.
#'
readEmblEntries <- function(dat){
  iLine <- itertools::ihasNext(dat)
  emblEntryList <- c()
  idx <- 0
  while (TRUE){
    if (!itertools::hasNext(iLine))
      return(emblEntryList)
    line <- iterators::nextElem(iLine)
    idx <- idx+1
    line
    if (startsWith(line, "ID")){
      line = strsplit(strsplit(line, ";")[[1]][1], "\\s+")[[1]][2]
      eEntry <- EmblEntry(line)
      seqHeader <- ""
      deleted <- FALSE
    } else if (startsWith(line, "DT")){
      line <- strsplit(line, "\\s+")[[1]]
      if (line[5] == "Created,"){
        eEntry$setDateAssigned(line[2])
      } else if (line[6] == "Updated,"){
        eEntry$setLastUpdated(line[2])
      } else if (line[5] == "Current"){
        eEntry$setReleaseStatus(line[2])
      }
    } else if (startsWith(line, "FT")){
      line <- strsplit(line, "\\s+")[[1]]
      if (startsWith(line[2], "/allele")){
        features <- dplyr::data_frame(feat=character(),
                                      start=integer(),
                                      end=integer(),
                                      count=integer())
        allele <- strsplit(line[2], '\"')[[1]][2]
        eEntry$setAlleleName(allele)
      } else if (line[2] %in% c("UTR", "intron", "exon")){
        positions <- strsplit(line[3], "\\.\\.")[[1]]
        feature <- list(feat = line[2],
                        start = as.integer(positions[1]),
                        end = as.integer(positions[2]))
        if (line[2] == "UTR"){
          feature["count"] <- ifelse(NROW(features) == 0, 1L,2L)
          features <- dplyr::bind_rows(features, feature)
        }
      } else if (startsWith(line[2], "/number")){
        count = strsplit(line[2], '\"')[[1]][2]
        feature["count"] <- as.integer(gsub("/","", count))
        features <- dplyr::bind_rows(features, feature)
      }
    } else if (startsWith(line, "SQ")){
      if (is.null(eEntry$getAlleleName()))
        deleted <- TRUE
      eEntry$setFeat(features)

      seq
      seqHeader = eEntry$getAlleleName()
      seq <- Biostrings::DNAString()
    } else if (line == "//"){
      seq <- Biostrings::DNAStringSet(seq)
      names(seq) <- seqHeader
      eEntry$setSeq(seq)
      if (!is.null(eEntry$getAlleleName()))
        emblEntryList <- c(emblEntryList, eEntry)
    } else if (!seqHeader == "" && !deleted){
      line <- strsplit(line, "\\s+")[[1]]
      seqP <- Biostrings::DNAString(paste0(line[1:(length(line)-1)],
                                           collapse = ""))
      seq <- c(seq, seqP)
    }
    TRUE
  }
  emblEntryList
}

.emblEntryListToDataFrame <- function(embl){
  allele <- embl$getAlleleName()
  seq <- data_frame(GID = allele %||% NA,
                    seq = as.character(embl$getSeq()))
  features <- embl$getFeat()
  features$feat <- simpleCap(features$feat)
  features <- data_frame(GID        = allele %||% NA,
                              start      = features$start,
                              end        = features$end,
                              feat_name  = features$feat_name,
                              id         = "NA",
                              feat_order = features$count,
                              type       = features$feat,
                              frame      = 0
                              )
  annot <- data_frame(GID = allele %||% NA,
                      locus = embl$getLocus(),
                      allele_name = allele %||% NA ,
                      allele_id = embl$getId(),
                      g_group = "NA",
                      p_group = "NA",
                      date_assigned = embl$getDateAssigned() %||% "NA",
                      first_released = embl$getFirstReleased() %||% "NA",
                      last_updated = embl$getLastUpdated(),
                      release_status = embl$getReleaseStatus(),
                      confirmed = embl$getConfirmed() %||% NA,
                      cwd_status = "NA",
                      complete = embl$getComplete())

 list(seq = seq, features = features, annot = annot)
}


# Class: EmblEntry ----------------------------------------------------------

#' Constructor for EmblEntry objects.
#'
#' @param allele_id the id of a kir allele
#' @return A EmblEntry object

EmblEntry <- function(allele_id) {
  EmblEntry_$new(allele_id)
}

#' Class \code{"EmblEntry"}
#'
#' @docType class
#' @keywords data internal
#' @return Object of \code{\link{R6Class}} representing an Embl Entry.
#' @section Methods:
#' \describe{
#'   \item{\code{x$new(locusname, db_version, ncores = parallel::detectCores(),
#'     with_dist = FALSE)}}{Create an object of this class.}
#'   \item{\code{x$get_db_version()}}{Get the IMGT/HLA database version.}
#'   \item{\code{x$get_locusname()}}{Get the name of the locus.}
#'   \item{\code{x$get_alleles(allele)}}{Get alleles.}
#'   \item{\code{x$get_closest_complete_neighbor(allele)}}{Get the complete
#'     allele that is closest at exon 2 to the query allele.}
#'   \item{\code{x$get_reference_sequence(allele)}}{Get the (imputed) reference
#'     sequence for allele.}
#' }

EmblEntry_ <- R6::R6Class(
  classname = "EmblEntry",
  public = list(
    initialize = function(allele_id) {
      private$id <- allele_id
    },
    print = function() {
      fmt0 <- "EmblEntry <%s>;\nLocus: <%s>\nAllele: <%s>\nComplete: <%s>"
      cat(sprintf(fmt0,
                  self$getId(),
                  self$getLocus(),
                  self$getAlleleName(),
                  self$getComplete()))
      invisible(self)
    },
    ## getters and setters
    ## get Allele name
    getAlleleName = function() {
      private$name
    },
    ## set Allele name
    setAlleleName = function(allele_name){
      assertthat::is.string(allele_name)
      private$name = allele_name
    },
    ## set sequence
    setSeq = function(seq){
      assertthat::assert_that(is(object = seq, class2 = "DNAStringSet"))
      private$seq = seq
      return(self)
    },
    ## get Sequence
    getSeq = function(){
      private$seq
    },
    ## set features
    setFeat = function(features){
      private$feat = features
      return(self)
    },
    ## get features
    getFeat = function(){
      .getFeatures(private$feat)
    },
    ## get Locus
    getLocus = function(){
      strsplit(private$name, "\\*")[[1]][1]
    },
    ## set Id
    setId = function(id){
      assertthat::is.string(id)
      private$id = id
      return(self)
    },
    ## get allele id
    getId = function(){
      private$id
    },
    ## set date assigned
    setDateAssigned = function(date){
      ## Todo add test for object
      private$date_assigned = date
      return(self)
    },
    ## get date assigned
    getDateAssigned = function(){
      private$date_assigned
    },
    ## set first_released
    setFirstReleased = function(date){
      ## Todo add test for object
      private$first_released = date
      return(self)
    },
    ## get first released
    getFirstReleased = function(){
      private$first_released
    },
    ## set release status
    setReleaseStatus = function(date){
      ## Todo add test for object
      private$release_status = date
      return(self)
    },
    ## get release status
    getReleaseStatus = function(){
      private$release_status
    },
    ## set last updated
    setLastUpdated = function(date){
      ## Todo add test for object
      private$last_updated = date
      return(self)
    },
    ## get last updated
    getLastUpdated= function(){
      private$last_updated
    },
    ## get confirmed
    setConfirmed = function(confirmed){
      assertthat::assert_that(is.logical(confirmed))
      private$confirmed = confirmed
      return(self)
    },
    ## get confirmed
    getConfirmed = function(){
      private$confirmed
    },
    ## get complete
    getComplete = function(){
      sum(self$getFeat()$feat == "UTR") == 2
    }
  ),
  private = list(
    name = NULL, # [character]; allele name
    seq  = NULL, # [DNAStringSet]: allele sequence
    feat = NULL, # [GRanges]: feature ranges
    lcn  = NULL, # [character]; locus name
    id   = NULL, # [character]: allele id
    date_assigned = NULL,
    first_released = NULL,
    last_updated = NULL,
    release_status = NULL,
    confirmed = NULL,
    cwd_status = NULL
  )
)
##### Helper #####
#####

.getFeatures <- function(feats){
  feats %>%
    dplyr::mutate(feat_name = ifelse(
      feat == "UTR",
      ifelse(count == 1,
             "3' UTR",
             "5' UTR"),
      paste(paste0(toupper(substring(feat, 1, 1)),
                   substring(feat, 2, nchar(feat))),
            count, sep = " " )))
}

#### Utils ####

.VALID_HLA_LOCI <- function(){
  c("HLA-A", "HLA-B", "HLA-C", "HLA-DPB1", "HLA-DQB1", "HLA-DRB1")
}
.VALID_KIR_LOCI <- function(){
    c("KIR2DL1",  "KIR2DL2",  "KIR2DL3",  "KIR2DL4", "KIR2DL5A", "KIR2DL5B",
      "KIR2DP1",  "KIR2DS1", "KIR2DS2",  "KIR2DS3",  "KIR2DS4",  "KIR2DS5",
      "KIR3DL1",  "KIR3DL2",  "KIR3DL3",  "KIR3DP1", "KIR3DS1")
  }

.expandHlaAllele <- function(x, locus = NULL) {
  if (is.null(locus)) {
    ifelse(!grepl("^HLA-\\S+", x), paste0("HLA-", x), x)
  } else {
    locus <- sub("HLA-", "", toupper(locus))
    pattern1 <- paste0("^HLA-", locus, "[*]\\d\\d\\d?:.+$")
    pattern2 <- paste0("^", locus, "[*]\\d\\d\\d?:.+$")
    pattern3 <- "^\\d\\d\\d?:.+$"
    ifelse(grepl(pattern1, x),
           x,
           ifelse(grepl(pattern2, x),
                  paste0("HLA-", x),
                  ifelse(grepl(pattern3, x),
                         paste0("HLA-", locus, "*", x), x)))
  }
}

.matchHlaLocus <- function(locus) {
  locus <- .expandHlaAllele(locus)
  match.arg(locus, .VALID_HLA_LOCI())
}

.expandKirAllele <- function(x, locus = NULL) {
  if (is.null(locus)) {
    ifelse(!grepl("^HLA-\\S+", x), paste0("HLA-", x), x)
  } else {
    locus <- sub("HLA-", "", toupper(locus))
    pattern1 <- paste0("^HLA-", locus, "[*]\\d\\d\\d?:.+$")
    pattern2 <- paste0("^", locus, "[*]\\d\\d\\d?:.+$")
    pattern3 <- "^\\d\\d\\d?:.+$"
    ifelse(grepl(pattern1, x),
           x,
           ifelse(grepl(pattern2, x),
                  paste0("HLA-", x),
                  ifelse(grepl(pattern3, x),
                         paste0("HLA-", locus, "*", x), x)))
  }
}

.matchKirLocus <- function(locus) {
  locus <- .expandKirAllele(locus)
  match.arg(locus, .VALID_KIR_LOCI())
}
simpleCap <- function(x) {
  sapply(x, function(e){
    s <- strsplit(e, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
  })
}

`%||%` <- function(a, b) {
  if (length(a) == 0) b else a
}

#' gather the Hla data

#'  fetch every information and create the dataframes
#' Get the HLA info
#' @return A data.frame with all information to be stored in the database
getHlaDf <- function(){
  hlaDb <- fetchIPD(family = "HLA")
  hladoc <- readHlaXml(hlaDb)
  hla <- parseHlaAlleles(hladoc)
  closestComplete <- getClosest(hla)
  hla$annot <- dplyr::left_join(hla$annot, closestComplete, by="GID")
  hla
}

#' gather the KIR data
#'
#' fetch every information and create the dataframes
#' Get the HLA info
#' @return A data.frame with all information to be stored in the database
getKirDf <- function() {
  kirDb <- fetchIPD(family = "KIR")
  kirdoc <- readKirDat(kirDb)
  pattern <- "^.*(\\d+\\.\\d+\\.\\d+).*"
  version <- sub(pattern, "\\1",
                 kirdoc[grep("CC   IPD-KIR Release Version", kirdoc)[1]])
  kirObj <- readEmblEntries(kirdoc)
  kirFullDFList <- lapply(kirObj, .emblEntryListToDataFrame)
  kirDf <- list()
  kirDf$seq <- do.call(rbind, lapply(kirFullDFList, function(x) x$seq))
  kirDf$features <- do.call(rbind,
                            lapply(kirFullDFList, function(x) x$features))
  kirDf$annot <- do.call(rbind, lapply(kirFullDFList, function(x) x$annot))
  closestComplete <- getClosest(kirDf)
  kirDf$annot <- dplyr::left_join(kirDf$annot, closestComplete, by="GID")
  kirDf$version <- version
  kirDf
}

#' Prepare the final database
#'
#' Some restructuring of the data to fit into the SQLite database.
#'
#' @param df the dataframe from either getKirDf() or getHlaDf().
#' @return a list containing the final data.frame and version information.
#'
prepareDb <- function(df){
  ipd <- list()
  ipd$annot <- df$annot %>%
    mutate(closest_complete = ifelse(is.na(closest_complete), allele_name,
                                     closest_complete)) %>%
    mutate(g_group = ifelse(is.na(g_group), "none", g_group)) %>%
    mutate(p_group = ifelse(is.na(p_group), "none", p_group)) %>%
    mutate(confirmed = ifelse(is.na(confirmed), "NA", confirmed)) %>%
    mutate(first_released = ifelse(is.na(first_released), "NA",
                                   first_released)) %>%
    mutate(date_assigned = ifelse(is.na(date_assigned), "NA",
                                  date_assigned)) %>%
    mutate(cwd_status = ifelse(is.na(cwd_status), "NA", cwd_status))
  ipd$features <- df$features
  ipd$seq <- df$seq
  list(ipd = ipd, version = df$version)
}

#' Get the Kir data ready
#'
#' @return a list containing the final data.frame and version information.
getKir <- function() {
  df <- getKirDf()
  prepareDb(df)
}

#' Get the Hla Data ready
#'
#' @return a list containing the final data.frame and version information.
getHla <- function() {
  df <- getHlaDf()
  prepareDb(df)
}

## Run the functions and get the actual data
hla <- getHla()
kir <- getKir()
