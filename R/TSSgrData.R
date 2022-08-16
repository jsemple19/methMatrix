#' Genomic Ranges for methMats data TSS
#'
#' Genomic ranges object for the transcription start site (TSS) that is within
#' the regions for which there is methylation matrix data in methMats data object
#' The metadata column "ID" contains the gene name for each TSS. This genomic ranges
#' object is of width 1 (just the TSS), and it should be extended to +- 250bp to define
#' the full region for which the methMats have been defined.
#'
#' @docType data
#'
#' @usage data(TSS)
#'
#' @format An object of class genomicRanges
#'
#' @keywords datasets
#'
#' @references Semple.J.I. & Meister P.
#'
#' @examples
#' data(TSS)
"TSS"
