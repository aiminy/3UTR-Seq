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
    input.file.dir <- readLines(input, n=1)

    #input.file.dir=row

    cat("please define the output file directory:\n")

    input<-file('stdin', 'r')
    output.file.dir <- readLines(input, n=1)

    cmd1="bsub -P bbc -J \"DogFT\" -o %J.DogFT.log -e %J.DogFT.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd2=paste0("Rscript ",R_lib,"/ThreeUTR/bin/rscript/convertbam2bed.r ",input.file.dir,
                " ",output.file.dir)

    system(paste0(cmd1," ",cmd2))

    cat("Finished conversion...\n")

  }else{quit()}
