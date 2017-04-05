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

#' parserreadfiles
#'
#' @param input.file.dir
#' @param input.file.type
#' @param output.file.dir
#'
#' @return
#' @export
#'
#' @examples
#' input.file.dir <- "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq"
#' output.file.dir <- "/Volumes/Bioinformatics$/Aimin_project"
#'
#' res <- parserreadfiles(input.file.dir=input.file.dir,"bam",output.file.dir=output.file.dir)
#'
parserreadfiles <- function (input.file.dir, input.file.type)
{
  dir.name = input.file.dir
  dir.name = reformatPath(dir.name)
  file.name = file.path(dir.name, dir(dir.name, recursive = TRUE))
  file.name.2 <- as.list(file.name)
  file.name.3 <- lapply(1:length(file.name.2), function(u,
                                                        input.file.type, file.name.2) {
    tmp = file.name.2
    x = tmp[[u]]
    path_name = dirname(x)
    file_name = basename(x)
    n <- length(grep(input.file.type, file_name))
    if (n == 1) {
      if(file_ext(file_name)==input.file.type){
      re <- file.path(path_name, file_name)
      } else{
        re <- NULL
    }

    }
    else {
      re <- NULL
    }
    re
  }, input.file.type, file.name.2)

  file.name.4 <- file.name.3[lapply(file.name.3, length) >
                               0]
  names(file.name.4) = unlist(lapply(1:length(file.name.4),
                                     function(u, file.name.4) {
                                       tmp = file.name.4
                                       x = tmp[[u]]
                                       path_name = dirname(x)
                                       file_name = basename(x)
                                       file_name
                                     }, file.name.4))

  #output.dir.name = reformatPath(output.file.dir)
  #temp3 = output.dir.name
  re2 <- list(input = file.name.4, input.file.type = input.file.type)
  return(re2)
}

#' convertbam2bed
#'
#' @param input.bamfile.dir
#' @param output.bedfile.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' input.bamfile.dir <- "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq"
#' output.bedfile.dir <- "/Volumes/Bioinformatics$/Aimin_project"
#'
#' res <- convertbam2bed(input.bamfile.dir,output.bedfile.dir)
#'
convertbam2bed <- function(input.bamfile.dir,output.bedfile.dir){

  res <- parserreadfiles(input.file.dir=input.bamfile.dir,
                         "bam")

  res <- res$input

  cmd0 <- "bedtools bamtobed -i"
  cmd1 <- ">"

  output.bedfile.dir <-file.path(output.bedfile.dir,"BedFileFromBam")

  if (!dir.exists(output.bedfile.dir)) {
    dir.create(output.bedfile.dir)
  }

  cmd.l <- lapply(res,function(u,output.bedfile.dir){

    # cat(u,"\n")
    # cmd9 <- "grep"
    # cmd10 <- "~/PathwaySplice/inst/extdata/"
    # cmd11 <- "/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
    # cmd12 <- ">"
    # cmd13 <- paste0("/Counts.",n,".genes.txt")
    # xxx <- gsub(";","",xx)
    #
    file_name = file_path_sans_ext(basename(u))

    cmd2 <- paste(cmd0,u,cmd1,file.path(output.bedfile.dir,paste0(file_name,".bed")),sep = " ")

    system(cmd2)

    cmd2
  },output.bedfile.dir)

  re <- list(cmdl = cmd.l,output.bedfile.dir = output.bedfile.dir)

  re

}
