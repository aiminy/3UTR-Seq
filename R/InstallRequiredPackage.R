#' InstallRequiredPackage
#'
#' @return
#' @export
#'
#' @examples
#' InstallRequiredPackage()
#'
InstallRequiredPackage<- function() {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Homo.sapiens")
  biocLite("hgu95av2.db")
  biocLite("rhdf5")
  library("Homo.sapiens")
  suppressMessages(library(hgu95av2.db))
  install.packages("plyr")
  library(plyr)
  biocLite("rtracklayer")
  library(rtracklayer)
  library(org.Hs.eg.db)
  library(annotate)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  install.packages("devtools")
  devtools::install_github("pachterlab/sleuth")
  library("sleuth")
  biocLite("biomaRt")
  install.packages("corrplot")
  library(corrplot)
  biocLite("sva")
  library(sva)
  #library(plyr)
}
