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
  install.packages("plyr")
  biocLite("rtracklayer")
  install.packages("devtools")
  devtools::install_github("pachterlab/sleuth")
  biocLite("biomaRt")
  install.packages("corrplot")
  biocLite("sva")
  source("https://bioconductor.org/biocLite.R")
  biocLite("ShortRead")
  install.packages("doMC")

  library(Homo.sapiens)
  library(plyr)
  library(rtracklayer)
  library(annotate)

}
