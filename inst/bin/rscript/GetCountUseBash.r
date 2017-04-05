#!/usr/bin/Rscript

#Usage:

#Rscript  Rscript ~/R/lib64/R/library/ThreeUTR/bin/rscript/GetCountUseBash.r

R_lib=.libPaths()[1]

cat("Do you want to get count?\n")

input<-file('stdin', 'r')
row <- readLines(input, n=1)

print(row)

if(row=="Yes") {

    cat("please define the input Bam file directory:\n")

    input<-file('stdin', 'r')
    input.file.dir <- readLines(file('stdin', 'r'), n=1)

    cat("please define the output file directory:\n")

    input<-file('stdin', 'r')
    output.file.dir <- readLines(input, n=1)

    cat("please define the annotation file:\n")

    input<-file('stdin', 'r')
    annotation.file <- readLines(input, n=1)

    cat("please define the upstream base of each transcript :\n")

    input<-file('stdin', 'r')
    ld <- readLines(input, n=1)

    cat("please define the downstream base of each transcript :\n")

    input<-file('stdin', 'r')
    rd <- readLines(input, n=1)

    cmd1="bsub -P bbc -J \"DogFT\" -o %J.DogFT.log -e %J.DogFT.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd2=paste("Rscript",R_lib,"/ThreeUTR/bin/rscript/convertbam2bed.r",
               input.bamfile.dir,
               annotation.bed.file,
               ld,
               rd,
               output.count.file.dir,
               sep=" ")

    system(paste0(cmd1," ",cmd2))

    cat("Finished conversion...\n")

  }else{quit()}
