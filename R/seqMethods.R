## Get the reference sequences as Biostrings::DNAStringSets
#
.getReference <- function(x, alleles) {
  cols <- c("seq")
  keytype <- "allele_name"
  res <- suppressMessages(
    select(x, keys = alleles, columns = cols, keytype = keytype))
  seqs <- DNAStringSet(res$seq)
  names(seqs) <- res$allele_name
  seqs
}

## Get gene structure as GenomicRanges::GRanges objects
##
.getStructure <- function(x, alleles) {
  cols <- c("start", "end", "type", "feat_name",
            "complete", "frame")
  keytype <- "allele_name"
  res <- suppressMessages(
    select(x, keys = alleles, columns = cols, keytype =keytype ))
  GRanges(seqnames = res$allele_name,
          ranges = IRanges(start = as.integer(res$start),
                           end = as.integer(res$end),
                           names = res$allele_name),
          strand = "+",
          type = res$type,
          name = res$feat_name,
          complete = as.logical(as.integer(res$complete)),
          frame  = as.integer(res$frame))
}

## Get the clostest complete reference
#
.getClosestComplete <- function(x, allele, locus = NULL) {
  assert_that(is.string(allele))
  if (!is.null(locus))
    allele <- .matchAllele(x, allele, locus)
  cols <- c("complete", "closest_complete")
  keytype <- "allele_name"
  res <- suppressMessages(
    select(x, keys = allele, columns = cols, keytype = keytype))
  if (NROW(res) == 0) {
    stop("Allele not found. Try again with adding the locus to the call")
  }
  closestComplete <- ifelse(length(allele) == 1,
                            res$closest_complete,
                            res[res$complete == 1, ][1, ]$closest_complete)

  getReference(x, closestComplete)
}

## Get all alleles of a locus
##
.getAlleles <- function(x, locus) {
  keytypes(x)
  keys(x, keytype = "TYPE")
  cols <- c("allele_name")
  res <- suppressMessages(
    select(x, keys = locus, keytype = "LOCUS", columns = cols))
  res$allele_name
}

## Match partial alleles
##
.matchAllele <- function(x, allele, locus) {
  assert_that(locus %in% getLoci(x))
  assert_that(is.character(allele))
  validAlleles <- getAlleles(x, locus)
  if (allele %in% validAlleles) {
    message(sprintf("Found allele %s", allele ))
    return(allele)
  } else if (any(matchingAlleles <- startsWith(validAlleles, allele))) {
    allele <- validAlleles[matchingAlleles]
    message(sprintf("Found allele %s by extension", allele))
    return(allele)
  } else if (any(matchingAlleles <- grepl(allele, validAlleles))) {
    allele <- validAlleles[matchingAlleles]
    message(sprintf("Found allele %s by regex matching", allele))
    return(allele)
  } else {
    message(
      sprintf("Could not find allele %s. Match to all of the locus", allele))
    return(validAlleles)
  }
}

## Get all loci
##
.getLoci <- function(x) {
  assert_that(is(x, "IpdDb"))
  unique(keys(x, keytype = "locus"))
}
