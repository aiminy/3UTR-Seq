#' Loadlibrary
#'
#' @return
#' @export
#'
#' @examples
#'
#' Loadlibrary()
#'
Loadlibrary <- function() {
  suppressPackageStartupMessages(library(Homo.sapiens))
  suppressMessages(library(hgu95av2.db))
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(rtracklayer))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(annotate))
  suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
  suppressPackageStartupMessages(library("sleuth"))
  suppressPackageStartupMessages(library(corrplot))
  suppressPackageStartupMessages(library(sva))
  suppressPackageStartupMessages(library(DESeq2))
}
