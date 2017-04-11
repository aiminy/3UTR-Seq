installercran <- function(.cran_packages)
{
    .inst <- .cran_packages %in% installed.packages()
    if (any(!.inst))
    {
        install.packages(.cran_packages[!.inst])
    }
    # Load packages into session, and print package version
    sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
}

installerbioc <- function(.bioc_packages)
{
    .inst <- .bioc_packages %in% installed.packages()
    if (any(!.inst))
    {
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
#' input.file.dir <- '/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq'
#' output.file.dir <- '/Volumes/Bioinformatics$/Aimin_project'
#'
#'
#' res <- parserreadfiles(input.file.dir,'bam',filter.sample="emp1")
#'
parserreadfiles <- function(input.file.dir,input.file.type,sample.group=NULL,filter.sample=NULL)
{
    dir.name = input.file.dir
    dir.name = reformatPath(dir.name)
    file.name = file.path(dir.name, dir(dir.name, recursive = TRUE))
    file.name.2 <- as.list(file.name)
    file.name.3 <- lapply(1:length(file.name.2), function(u, input.file.type, file.name.2)
    {
        tmp = file.name.2
        x = tmp[[u]]
        path_name = dirname(x)
        file_name = basename(x)
        n <- length(grep(input.file.type, file_name))
        if (n == 1)
        {
            if (file_ext(file_name) == input.file.type)
            {
                re <- file.path(path_name, file_name)
            } else
            {
                re <- NULL
            }

        } else
        {
            re <- NULL
        }
        re
    }, input.file.type, file.name.2)

    file.name.4 <- file.name.3[lapply(file.name.3, length) > 0]
    names(file.name.4) = unlist(lapply(1:length(file.name.4), function(u, file.name.4)
    {
        tmp = file.name.4
        x = tmp[[u]]
        path_name = dirname(x)
        path_name2 <- basename(path_name)

        file_name = basename(x)

        file_name <- paste0(path_name2,"-",file_name)

        file_name

    }, file.name.4))

    if(!is.null(filter.sample)){

     file.name.5 <- file.name.4[-grep(filter.sample,names(file.name.4))]

    }else{

     file.name.5 <- file.name.4
     }


    if(!is.null(sample.group)){
    g1 <- grep(sample.group[1],toupper(names(res$input)))

    g2 <- grep(sample.group[2],toupper(names(res$input)))

    # output.dir.name = reformatPath(output.file.dir) temp3 = output.dir.name
    re2 <- list(input = file.name.5,g1=g1,g2=g2,input.file.type = input.file.type)}

    else
    {
      re2 <- list(input = file.name.5,input.file.type = input.file.type)
    }

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
#' input.bamfile.dir <- '/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq'
#' output.bedfile.dir <- '/Volumes/Bioinformatics$/Aimin_project'
#'
#' res <- convertbam2bed(input.bamfile.dir,output.bedfile.dir)
#'
convertbam2bed <- function(input.bamfile.dir, output.bedfile.dir)
{
    res <- parserreadfiles(input.bamfile.dir,"bam")

    res <- res$input

    cmd0 <- "bedtools bamtobed -i"
    cmd1 <- ">"

    output.bedfile.dir <- file.path(output.bedfile.dir, "BedFileFromBam")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir)
    }

    cmd.l <- lapply(res, function(u, output.bedfile.dir)
    {

        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11 <-
        # '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13 <-
        # paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)

         path_name = dirname(u)
         path_name2 <- basename(path_name)

         file_name = file_path_sans_ext(basename(u))

         file_name <- paste0(path_name2,"-",file_name)

        cmd2 <- paste(cmd0, u, cmd1, file.path(output.bedfile.dir, paste0(file_name,
            ".bed")), sep = " ")

        system(cmd2)

        cmd2
    }, output.bedfile.dir)

    re <- list(cmdl = cmd.l, output.bedfile.dir = output.bedfile.dir)

    re

}

#' matchbed2annotation
#'
#' @param input.bedfile.dir
#' @param annotation.bed.file
#' @param ld upstream base
#' @param rd downstream base
#' @param output.matched.bed.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' res <- convertbam2bed(input.bamfile.dir,output.bedfile.dir)
#' input.bedfile.dir <- res$output.bedfile.dir
#' annotation.bed.file <- annotation.bed.file
#' ld <- ld
#' rd <- rd
#' output.matched.bed.file.dir <- output.matched.bed.file.dir
#'
#' res <- matchbed2annotation(input.bedfile.dir,annotation.bed.file,
#' ld,rd,output.matched.bed.file.dir)
#'
matchbed2annotation <- function(input.bedfile.dir, annotation.bed.file, ld, rd, output.matched.bed.file.dir)
{

    res <- parserreadfiles(input.bedfile.dir, "bed")

    res <- res$input

    cmd0 <- paste("bedtools window -a", annotation.bed.file, "-b", sep = " ")

    cmd1 <- paste("-l", ld, "-r", rd, "-sw", ">", sep = " ")

    output.bedfile.dir <- file.path(output.matched.bed.file.dir, "MatchedBedFile")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir)
    }

    cmd.l <- lapply(res, function(u, output.bedfile.dir)
    {

        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11 <-
        # '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13 <-
        # paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)
        file_name = file_path_sans_ext(basename(u))

        cmd2 <- paste(cmd0, u, cmd1, file.path(output.bedfile.dir, paste0(file_name,
            "_matched.bed")), sep = " ")

        system(cmd2)

        cmd2
    }, output.bedfile.dir)

    re <- list(cmdl = cmd.l, output.bedfile.dir = output.bedfile.dir)

    re

}

#' getcountsfromMatchedbed
#'
#' @param input.bedfile.dir
#' @param annotation.bed.file
#' @param ld upstream base
#' @param rd downstream base
#' @param output.matched.bed.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' res <- matchbed2annotation(input.bedfile.dir,annotation.bed.file,
#' ld,rd,output.matched.bed.file.dir)
#'
#' input.bedfile.dir <- res$output.bedfile.dir
#'
#' res <- getcountsfromMatchedbed (input.bedfile.dir,output.count.file.dir)
#'
getcountsfromMatchedbed <- function(input.bedfile.dir,output.count.file.dir,filter.sample)
{

    res <- parserreadfiles(input.bedfile.dir,"bed", filter.sample=filter.sample)

    res <- res$input

    # system("awk -F '\\t' '$6==\"+\" && $12==\"-\"' ~/MatchedBedFile/R1_Dox_matched.bed | awk '$8<$2&&$9>=$2' | awk '{print $4}' | sort | uniq -c | sort -nr | head")

    cmd0 <- "awk -F '\\t'"

    cmd1 <- "'$6==\"+\" && $12==\"-\"'"
    cmd11 <- "'$6==\"+\" && $12==\"+\"'"
    cmd12 <- "'$6==\"-\" && $12==\"+\"'"
    cmd13 <- "'$6==\"-\" && $12==\"-\"'"

    cmd2 <- "| awk '$8 < $2 && $9 >= $2' | awk '{print $4}' | sort | uniq -c | sort -nr"
    cmd21 <- "| awk '$8 >= $2 && $9 <= $3' | awk '{print $4}' | sort | uniq -c | sort -nr"
    cmd22 <- "| awk '$8 <= $3 && $9 > $3' | awk '{print $4}' | sort | uniq -c | sort -nr"

    cmd3 <- ">"

    cmdtemp<-rbind(
    cbind(rep(cmd1,3),c(cmd2,cmd21,cmd22),rep("plus",3),rep("minus",3),c("below.DoGs","DoGs","over.DoGs")),
    cbind(rep(cmd11,3),c(cmd2,cmd21,cmd22),rep("plus",3),rep("plus",3),c("below.DoGs","DoGs","over.DoGs")),
    cbind(rep(cmd12,3),c(cmd2,cmd21,cmd22),rep("minus",3),rep("plus",3),c("over.DoGs","DoGs","below.DoGs")),
    cbind(rep(cmd13,3),c(cmd2,cmd21,cmd22),rep("minus",3),rep("minus",3),c("over.DoGs","DoGs","below.DoGs")))

    output.count.file.dir <- file.path(output.count.file.dir, "Counts")

    if (!dir.exists(output.count.file.dir))
    {
        dir.create(output.count.file.dir)
    }

    counteachcase <- function(res, cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand, location, output.count.file.dir) {

      cmd.l <- lapply(res, function(u, output.count.file.dir)
      {

          file_name = file_path_sans_ext(basename(u))

          cmd4 <- paste(cmd0, cmd1, u, cmd2, cmd3, file.path(output.count.file.dir,
              paste0(file_name,".",gene.strand,".gene.",read.strand,".read.",
                     location,".count.txt")), sep = " ")
          cat(cmd4,"\n")

          system(cmd4)

          cmd4

      }, output.count.file.dir)

     return(cmd.l)

    }


    cmdtempres2 <- apply(cmdtemp,1,function(u,cmd0,cmd3,output.count.file.dir){
     x <- as.data.frame(t(u))

     cmd1 <- x[,1]
     cmd2 <- x[,2]
     gene.strand <- x[,3]

     read.strand <- x[,4]

     location <- x[,5]

     cmdtempres <- counteachcase(res, cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand, location, output.count.file.dir)

     cmdtempres

    },cmd0,cmd3,output.count.file.dir)

    re <- list(cmdl = cmdtempres2, output.count.file.dir = output.count.file.dir)

    re

}

#' getcounts
#'
#' @param input.bamfile.dir
#' @param annotation.bed.file
#' @param ld
#' @param rd
#' @param output.count.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' getcounts()
#'
#'
#'
getcounts <- function(input.bamfile.dir,annotation.bed.file,ld,rd,output.count.file.dir,filter.sample)
{

    res <- convertbam2bed(input.bamfile.dir, output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    annotation.bed.file <- annotation.bed.file
    ld <- ld
    rd <- rd

    res <- matchbed2annotation(input.bedfile.dir, annotation.bed.file, ld, rd, output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    res <- getcountsfromMatchedbed(input.bedfile.dir,output.count.file.dir,filter.sample)

    res
}

#' parserAnnotationFile
#'
#' @return
#' @export
#'
#' @examples
#'
#' input.annotation.file <- "/Volumes/Bioinformatics$/Aimin_project/UTR/3UTR_GFF2.txt"
#'
#' y <- parserAnnotationFile(input.annotation.file)
#'
#'
parserAnnotationFile <- function(input.annotation.file){

   dir.name <- dirname(input.annotation.file)

   file.name <- file_path_sans_ext(basename(input.annotation.file))

   x <- read.table(input.annotation.file)

   d <- x[,5]-x[,4]

   nn <- substr(x[,9],9,23)

   xx <-cbind.data.frame(x[,1],x[,4],x[,5],nn,d,x[,7])

   write.table(xx, file = file.path(dir.name, paste0(file.name,
                                                        ".bed")),
             row.names = FALSE,
             col.names =FALSE,
             quote = FALSE,
             sep ="\t"
             )

   xx
}

parsersample <- function(input.bamfile.dir) {

  xx <- parserreadfiles(input.bamfile.dir,"bam")

  cell<-factor(rep(c('emp','hela'),Nbatch))
  cell=rep(cell,2)

  a <- length(wt.index)
  b <- length(dox.index)

  sample <- c("R1_Dox.bam",
              "R2_Dox.bam",
              "R3_Dox.bam",
              "R4_Dox.bam",
              "R5_Dox.bam",
              "R1_WT.bam",
              "R2_WT.bam",
              "R3_WT.bam",
              "R4_WT.bam",
              "R5_WT.bam")

  colData <- data.frame(sample=sample,condition=factor(rep(c('Dox', 'WT'),c(5, 5))),cell=cell)

  colData
}

