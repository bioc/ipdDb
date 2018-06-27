#' Get allele information from IPD
#' 
#' This package holds the IPD IMGT/HLA and IPD KIR database.
#' All alleles are accessible using the select, columns, keys and keytypes
#' methods of the AnnoatationDbi package of bioconductor.
#' 
#' Included data are:
#'   
#' Allele names
#' 
#' p-groups
#' 
#' g-groups
#' 
#' cwd_status
#' 
#' completeness status
#' 
#' gene structure
#' 
#' reference sequences
#' 
#' closest full-length allele
#' 
#' @docType package
#' @name ipdDbPackage
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats setNames
#' @import methods
#' @import AnnotationDbi
#' @import RSQLite
#' @import DBI
#' @import assertthat
#' @import AnnotationHub

NULL
