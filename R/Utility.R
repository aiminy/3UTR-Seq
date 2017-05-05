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

parsersample <-function()
  {

  cell<-factor(rep(c('batch1','batch2'),c(3,2)))
  cell=rep(cell,2)

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

    pkg.env <- new.env(parent = emptyenv())
    pkg.env$sample <- ThreeUTR:::parsersample()

    return(re2)


}

useWget2Download <- function(sra.accession.number,output.dir){

  cmd0 <- "wget -c -r -nd -np -L"
  cmd1 <- "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/"

  temp <- file.path(substr(sra.accession.number,1,6),sra.accession.number)

  temp2 <- paste0(cmd1,temp)

  cmd2 <- paste(cmd0,temp2,"-P",output.dir,sep=" ")

  system(cmd2)

}

useFastqDump <- function(sra.accession.number,output.dir){

  cmd0 <- "fastq-dump -I --split-files"

  cmd1 <- paste(cmd0,sra.accession.number,"-O",output.dir,sep=" ")

  system(cmd1)

}

useFastqDumpConvertSra2Fastq <- function(sra.file.dir,output.dir){

  cmd0 <- "fastq-dump --split-3"

  re <- parserreadfiles(sra.file.dir,'sra')

  res <- re$input

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir)
  }

  cmd.l <- lapply(res, function(u, output.dir)
  {

    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    cmd1 <- paste(cmd0,u,"-O",output.dir,sep=" ")

    system(cmd1)

    cmd1
  }, output.dir)

  re <- list(cmdl = cmd.l, output.dir = output.dir)

  re

}

#infer_experiment.py
useInferExperiment<-function(input.file.dir,ref.gene.bed.file,output.dir){

  re <- parserreadfiles(input.file.dir,'bam')

  res <- re$input

  cmd0 <- "infer_experiment.py -i"
  cmd1 <- "-r"
  cmd2 <- ">"

  #output.dir <- file.path(output.dir, "BamInfo")

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir)
  }

  cmd.l <- lapply(res, function(u, output.dir)
  {

    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    cmd2 <- paste(cmd0,u,cmd1,ref.gene.bed.file,
                  cmd2,file.path(output.dir, paste0(file_name,"_infor.txt")), sep = " ")

    system(cmd2)

    cmd2
  }, output.dir)

  re <- list(cmdl = cmd.l, output.dir = output.dir)

  re
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

    cmd2 <- "| awk '$8 < $2 && $9 >= $2' | awk '{print $4}' | sort | uniq -c | sort -nr"   #below
    cmd21 <- "| awk '$8 >= $2 && $9 <= $3' | awk '{print $4}' | sort | uniq -c | sort -nr" #DoGs
    cmd22 <- "| awk '$8 <= $3 && $9 > $3' | awk '{print $4}' | sort | uniq -c | sort -nr"  #over

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

    xx <-cbind.data.frame(x,nn,d)

    counts <- as.data.frame(table(xx$nn))

    colnames(counts)=c("nn","counts")

    xxx <- merge(xx,counts,by="nn",sort=FALSE)

    xxxx <- xxx[,c(2,5,6,10,11,8)]

    write.table(xxxx, file = file.path(dir.name, paste0(file.name,
                                                        ".bed")),
              row.names = FALSE,
              col.names =FALSE,
              quote = FALSE,
              sep ="\t"
              )

    xxxx
}

convertBam2Bw <- function(input.bam.file.dir,input.chromosome.size.file,output.bw.file.dir){

  re <- parserreadfiles(input.bam.file.dir,'bam')

  res <- re$input

  cmd0="samtools sort"

  if (!dir.exists(output.bw.file.dir))
  {
    dir.create(output.bw.file.dir)
  }

  cmd.l <- lapply(res, function(u, output.bw.file.dir)
  {
    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    cmd <- paste(cmd0,u,file.path(output.bw.file.dir, paste0(file_name,"_sorted")), sep = " ")

    system(cmd)

    cmd
  }, output.bw.file.dir)

  cmd1="samtools index"

  cmd.l <- lapply(res, function(u, output.bw.file.dir)
  {

    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    cmd <- paste(cmd1,file.path(output.bw.file.dir, paste0(file_name,"_sorted.bam")), sep = " ")

    system(cmd)

    cmd
  }, output.bw.file.dir)

  cmd3 <- "tail -n +2"

  path_name = dirname(input.chromosome.size.file)
  path_name2 <- basename(path_name)

  file_name = file_path_sans_ext(basename(input.chromosome.size.file))

  file_name <- paste0(path_name2,"-",file_name,".txt")

  cmd <- paste(cmd3,input.chromosome.size.file,">",file.path(path_name,file_name), sep = " ")

  system(cmd)

  cmd4 <- "genomeCoverageBed -ibam"
  cmd5 <- "-bg -g"
  cmd6 <- ">"

  cmd.l <- lapply(res, function(u, output.bw.file.dir)
  {
    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    x = file.path(output.bw.file.dir, paste0(file_name,"_sorted.bam"))

    cmd <- paste(cmd4,x,cmd5,input.chromosome.size.file,
                  cmd6,file.path(output.bw.file.dir, paste0(file_name,".bdg")), sep = " ")

    system(cmd)

    cmd
  }, output.bw.file.dir)

  #re <- list(cmdl = cmd.l, output.dir = output.dir)

  #input.bdg.file.dir <- re$output.dir

  cmd7 <- "LC_COLLATE=C sort -k1,1 -k2,2n"

  #re <- parserreadfiles(input.bdg.file.dir,'bdg')

  #res <- re$input

  cmd.l <- lapply(res, function(u, output.bw.file.dir)
  {
    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    x = file.path(output.bw.file.dir, paste0(file_name,".bdg"))

    cmd <- paste(cmd7,x,cmd6,file.path(output.bw.file.dir, paste0(file_name,".sorted_bdg")), sep = " ")

    system(cmd)

    cmd
  }, output.bw.file.dir)

  cmd8 = "bedGraphToBigWig"

  input.chromosome.size.file.m <- file.path(path_name,file_name)

  #re <- parserreadfiles(input.bdg.file.dir,'.sorted_bdg')

  #res <- re$input

  cmd.l <- lapply(res, function(u,input.chromosome.size.file.m,output.bw.file.dir)
  {
    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    x = file.path(output.bw.file.dir,paste0(file_name,".sorted_bdg"))

    cmd <- paste(cmd8,x,input.chromosome.size.file.m,file.path(output.bw.file.dir,paste0(file_name,".bw")), sep = " ")

    system(cmd)

    cmd
  },input.chromosome.size.file.m,output.bw.file.dir)

  cmd9 = "bigWigToWig"

  cmd.l <- lapply(res, function(u,output.bw.file.dir)
  {
    path_name = dirname(u)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    file_name <- paste0(path_name2,"-",file_name)

    x = file.path(output.bw.file.dir,paste0(file_name,".bw"))

    cmd <- paste(cmd9,x,file.path(output.bw.file.dir,paste0(file_name,".wig")), sep = " ")

    system(cmd)

    cmd
  },output.bw.file.dir)

}

prepareDaPars<-function(input.wig.file.dir,sample.group=c("Dox","WT"),output.dir){

  re <- parserreadfiles(input.wig.file.dir,'wig')

  res <- re$input

  cmd = "tail -n +2"

  cmd.l <- lapply(res, function(u,output.dir)
  {
    path_name = dirname(u)
    #path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(u))

    #file_name <- paste0(path_name2,"-",file_name)

    x = u

    cmd <- paste(cmd,x,">",file.path(output.dir,paste0(file_name,"2.wig")), sep = " ")

    system(cmd)

    cmd
  },output.dir)

  g1 <- "1,2,3,4,5"
  g2 <- "a,b,c,d,e"

cat(c("Annotated_3UTR=hg19_refseq_extracted_3UTR.bed\n",
paste0("Group1_Tophat_aligned_Wig=",g1,"\n"),
paste0("Group2_Tophat_aligned_Wig=",g2,"\n"),
"Output_directory=DaPars_Dox_WT_data/\n",
"Output_result_file=DaPars_Dox_WT_data\n",
"#Parameters\n
Num_least_in_group1=1\n
Num_least_in_group2=1\n
Coverage_cutoff=30\n
FDR_cutoff=0.05\n
PDUI_cutoff=0.5\n
Fold_change_cutoff=0.59\n"),file=file.path(output.dir,"config_test.txt"),sep="")
}

select3UTR <- function(genome,tablename) {
  library(GenomicFeatures)
  library(dplyr)
  refSeq             <- makeTxDbFromUCSC(genome="hg19",tablename="knownGene")
  threeUTRs          <- threeUTRsByTranscript(refSeq, use.names=TRUE)
  length_threeUTRs   <- width(ranges(threeUTRs))
  the_lengths        <- as.data.frame(length_threeUTRs)
  the_lengths        <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
  the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
  colnames(the_lengths) <- c("RefSeq Transcript", "3' UTR Length")

}

subsetFastq<-function(input.fastq.files.dir,output.dir,n,gene.model.file,genome.index){

  re <- parserreadfiles(input.fastq.files.dir,'fastq')

  res <- re$input

  cmd0 = paste0("head -",n)
  cmd1 = ">"

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir)
  }

  # cmd1="bsub -P bbc -J \"STAR-alignment\" -o %J.STAR-alignment.log -e %J.STAR-alignment.err -W"
  # cmd2="72:00 -n 8 -q bigmem -R 'rusage[mem=36864] span[hosts=1]' -u aimin.yan@med.miami.edu"
  #
  # cmd3=paste("STAR",strand,"--genomeLoad NoSharedMemory --runThreadN",Ncores,"--sjdbGTFfile",collapse = " ")
  # cmd4=paste(input.gtf.file,"--outFileNamePrefix",collapse = " ")
  #
  # cmd44=paste("--genomeDir",STAR.index.file,collapse = " ")
  #
  #
  # cmd11="bsub -w \"done(\"STAR-alignment\")\" -P bbc -J \"samtools-sort\" -o %J.samtools-sort.log -e %J.samtools-sort.err -W"

  xx <-lapply(res,function(u,output.dir){

    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    sample.name.out = file.path(output.dir,paste0(file_name,"-test-",n,".fastq"))

    cmd3= paste(cmd0,u,cmd1,sample.name.out)

    cmd3

    system(cmd3)

  },output.dir)


  re <- parserreadfiles(output.dir,'fastq')

  res <- re$input

  xx <-lapply(res,function(u){

    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    if(regexpr(pattern ='_',file_name)!=-1){

      cat("match:",file_name,"\n")

      p <- regexpr(pattern ='_',file_name)
      pp <- p-1
      x <- substr(file_name,1,pp)

    }else{
      cat("no match:",file_name,"\n")
      x <- file_name
    }

    x
  })


  print(xx)

  xxx <- unique(unlist(xx))
  res2 <- unlist(res)

  print(xxx)
  #print(res2)
  #  tophat -G genes.gtf -p 4 -o "201348193-01"_tophat_out mm10_index_bt2/genome ~/RNAseqData/"nBishopric_Project1_201348193
  #-01_S_1_1.txt" ~/RNAseqData/"nBishopric_Project1_201348193-01_S_1_2.txt"

  cmd0 = "tophat -G"
  cmd1 = "-p 4 -o"

  cmd2 = "mv"

  for(i in 1:length(xxx)){

    sample.name <- xxx[i]
    sample.name.out.dir <-file.path(output.dir,sample.name)

    if (!dir.exists(sample.name.out.dir))
    {
      dir.create(sample.name.out.dir)
    }

    y <- res2[grep(xxx[i],res2)]

    print(y)
    print(length(y))

    if(length(y)==2){
      cmd3= paste(cmd0,gene.model.file,cmd1,sample.name.out.dir,genome.index,y[1],y[2],sep=" ")
    }else
    {
      cmd3= paste(cmd0,gene.model.file,cmd1,sample.name.out.dir,genome.index,y[1],sep=" ")
    }


    print(cmd3)

    }


}

useTophat4Alignment<-function(input.fastq.files.dir,output.dir,gene.model.file=NULL,genome.index=NULL){

  re <- parserreadfiles(input.fastq.files.dir,'fastq')

  res <- re$input

  xx <-lapply(res,function(u){

    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    if(regexpr(pattern ='_',file_name)!=-1){

      cat("match:",file_name,"\n")

      p <- regexpr(pattern ='_',file_name)
      pp <- p-1
      x <- substr(file_name,1,pp)

    }else{
      cat("no match:",file_name,"\n")
      x <- file_name
    }

    x
  })


  print(xx)

  xxx <- unique(unlist(xx))
  res2 <- unlist(res)

  print(xxx)
  #print(res2)
#  tophat -G genes.gtf -p 4 -o "201348193-01"_tophat_out mm10_index_bt2/genome ~/RNAseqData/"nBishopric_Project1_201348193
  #-01_S_1_1.txt" ~/RNAseqData/"nBishopric_Project1_201348193-01_S_1_2.txt"

  cmd0 = "tophat -G"
  cmd1 = "-p 4 -o"

  cmd2 = "mv"

  for(i in 1:length(xxx)){

    sample.name <- xxx[i]
    sample.name.out.dir <-file.path(output.dir,sample.name)

    if (!dir.exists(sample.name.out.dir))
    {
      dir.create(sample.name.out.dir)
    }

    y <- res2[grep(xxx[i],res2)]

    print(y)
    print(length(y))

    if(length(y)==2){
      cmd3= paste(cmd0,gene.model.file,cmd1,sample.name.out.dir,genome.index,y[1],y[2],sep=" ")
    }else
    {
      cmd3= paste(cmd0,gene.model.file,cmd1,sample.name.out.dir,genome.index,y[1],sep=" ")
    }


    print(cmd3)
    #system(cmd3)

    #cmd4=paste(cmd2,paste0(sample.name.out.dir,"accepted_hits.bam"),file.path(output.dir,paste0(xxx[i],".bam")),sep=" ")
    #system(cmd4)
  }

}
