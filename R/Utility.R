
installercran <- function(.cran_packages) {
  .inst <- .cran_packages %in% installed.packages()
  if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
  }
  # Load packages into session, and print package version
  sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
}

installerbioc <- function(.bioc_packages){
  .inst <- .bioc_packages %in% installed.packages()
  if(any(!.inst)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(.bioc_packages[!.inst], ask = F)
  }
  # Load packages into session, and print package version
  sapply(.bioc_packages, require, character.only = TRUE)
}
