#!/usr/bin/Rscript

# Usage:

# Rscript ~/R/lib64/R/library/ThreeUTR/bin/rscript/utr_interactive.r

R_lib = .libPaths()[1]

get_os <- function()
{
  sysinf <- Sys.info()
  if (!is.null(sysinf))
  {
    os <- sysinf["sysname"]
    if (os == "Darwin")
      os <- "osx"
  } else
  {
    ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

os <- get_os()

cat("Your operating system is: ", os, "\n")

cat("Choose analysis: \n")
cat("\tAvaliable analysis: \n
    1. Download_SRA \n
    2. SRA2Fastq \n
    3. subset_Fastq_file \n
    4. testAlignment \n
    5. Alignment \n
    6. QC \n
    7. Bam2Bw \n
    8. Counts \n
    9. DE_analysis \n
    10. All\n")

input <- file("stdin", "r")
row <- readLines(input, n = 1)

if (row == 1)
{
  choose.type <- "Download_SRA"
}

if (row == 2)
{
  choose.type <- "SRA2Fastq"
}

if (row == 3)
{
  choose.type <- "subset_Fastq_file"
}

if (row == 4)
{
  choose.type <- "testAlignment"
}

if (row == 5)
{
  choose.type <- "Alignment"
}

if (row == 6)
{
  choose.type <- "QC"
}

if (row == 7)
{
  choose.type <- "Bam2Bw"
}

if (row == 8)
{
  choose.type <- "Counts"
}

if (row == 9)
{
  choose.type <- "DE_analysis"
}

if (row == 10)
{
  choose.type <- "All"
}

downloadSRA <- function(R_lib)
{
  cat("You choose to download SRA files, please define the following setting parameters: \n",
      "sra.accession.number (ex: SRP058633 Use an example from DoGs paper)\n",
      "output.dir (ex:DoGsExample)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 2)

  sra.accession.number <- count.file.dir[1]
  output.dir <- count.file.dir[2]

  library(ChipSeq)
  library(ThreeUTR)

  cmd1 = "bsub -P bbc -J \"downloadSRA\" -o %J.downloadSRA.log -e %J.downloadSRA.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  cmd2 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/downloadSRA.r"),
               sra.accession.number,output.dir,sep=" ")

  cmd3 = paste(cmd1,cmd2,sep=" ")

  print(cmd3)

  system(cmd3,intern= TRUE)

  cat("Finished downloadSRA\n")

}

sra2Fastq <- function(R_lib)
{
  cat("You choose to convert SRA files to Fastq files, please define the following setting parameters: \n",
      "sra.file.dir (ex: /nethome/axy148/DoGsExample)\n",
      "output.dir (ex:DoGsFastq)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 2)

  sra.file.dir <- count.file.dir[1]
  output.dir <- count.file.dir[2]

  library(ChipSeq)
  library(ThreeUTR)

  cmd1 = "bsub -P bbc -J \"sra2Fastq\" -o %J.sra2Fastq.log -e %J.sra2Fastq.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  cmd2 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/sra2Fastq.r"),
               sra.file.dir,output.dir,sep=" ")

  cmd3 = paste(cmd1,cmd2,sep=" ")

  print(cmd3)

  system(cmd3,intern= TRUE)

  cat("Finished sra2Fastq\n")

}

subsetFastqFile <- function(R_lib)
{
  cat("You choose to sample Fastq files, please define the following setting parameters: \n",
      "fastq.file.dir (ex: /scratch/projects/bbc/aiminy_project/DoGsFastq)\n",
      "output.dir (ex:/scratch/projects/bbc/aiminy_project/DoGsFastqTest)\n",
      "Number of short reads (ex: 20. it has to be 4n)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 3)

  fastq.file.dir <- count.file.dir[1]
  output.dir <- count.file.dir[2]
  num.read <- count.file.dir[3]

  library(ChipSeq)
  library(ThreeUTR)

  # cmd1 = "bsub -P bbc -J \"sra2Fastq\" -o %J.sra2Fastq.log -e %J.sra2Fastq.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  #cmd1 ="72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  #cmd2 = "bsub -P bbc -J \"subSetFastqAll\" -o %J.subSetFastqAll.log -e %J.subSetFastqAll.err -W"

  cmd3 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/subsetFastq.r"),
               fastq.file.dir,output.dir,num.read,sep=" ")

  #cmd4 = paste(cmd2,cmd1,cmd3,sep=" ")

  #print(cmd4)

  system(cmd3,intern= TRUE,wait = TRUE)

  cat("Finished subset fastq\n")

}

testAlignment <- function(R_lib) {

  cat("You choose to test Alignment, please define the following setting parameters: \n",
      "fastq.file.dir (ex: /scratch/projects/bbc/aiminy_project/DoGsFastq)\n",
      "gene.model.file(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf)\n",
      "genome.index(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome)\n")

  #cmd5="bsub -w \"done(\"subSetFastqAll\")\" -P bbc -J \"tophatALL\" -o %J.tophatALL.log -e %J.tophatALL.err -W"

  library(ChipSeq)
  library(ThreeUTR)

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 3)

  fastq.file.dir <- count.file.dir[1]
  gene.model.file <- count.file.dir[2]
  genome.index <- count.file.dir[3]

  cmd6 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/testAlignment.r"),
               fastq.file.dir,gene.model.file,genome.index,sep=" ")

  #cmd7 = paste(cmd5,cmd1,cmd6,sep=" ")

  system(cmd6,intern= TRUE)

  cat("Finished testAlignment fastq\n")
}

performAlignment <- function(R_lib)
{
  cat("You choose to perform alignment, please define the following setting parameters: \n",
      "fastq.file.dir (ex:/scratch/projects/bbc/aiminy_project/DoGsFastq)\n",
      "output.dir (ex:/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBam)\n",
      "gene.model.file(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf or /projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/refGene.txt\n",
      "genome.index(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 4)

  fastq.file.dir <- count.file.dir[1]
  output.dir <- count.file.dir[2]
  gene.model.file <- count.file.dir[3]
  genome.index <- count.file.dir[4]

  cat("Do you want to use bigmem to run your job:\n")

  input <- file("stdin", "r")
  general.or.bigmem <- readLines(input, n = 1)



  if(toupper(general.or.bigmem) == "YES")
  {
#    cmd1 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W 72:00 -n 16 -q bigmem -R 'rusage[mem=36864] span[ptile=8]' -u aimin.yan@med.miami.edu"
    cmd1 <- "BigMEM"
  }else
  {
#    cmd1 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
    cmd1 = "General"
  }

  #cmd4="bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W"

  #library(ChipSeq)
  #library(ThreeUTR)

  cmd2 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/alignment.r"),
               fastq.file.dir,output.dir,gene.model.file,genome.index,cmd1,sep=" ")

  #cmd3 = paste(cmd1,cmd2,sep=" ")

  #print(cmd3)

  cmd3 =cmd2

  system(cmd3,intern= TRUE)

  cat("Finished alignment\n")

}

bamQC <- function(R_lib)
{
  cat("You choose to perform QC on the bam files, please define the following setting parameters: \n",
      "input.bam.files.dir (ex: /scratch/projects/exempt/tr/azhang/for_bioinfo_core/3_end_seq)\n",
      "input.ref.gene.bed.file (ex: /projects/ctsi/bbc/aimin/annotation/3UTR_GFF2_2.bed)\n",
      "output.dir (ex:BAMQC)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 3)

  input.bam.file.dir <- count.file.dir[1]
  input.ref.gene.bed.file <- count.file.dir[2]
  output.dir <- count.file.dir[3]

  library(ChipSeq)
  library(ThreeUTR)


  cmd1 = "bsub -P bbc -J \"QCBam\" -o %J.QCBam.log -e %J.QCBam.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  cmd2 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/qcBAM.r"),
                input.bam.file.dir,input.ref.gene.bed.file,output.dir,sep=" ")

  cmd3 = paste(cmd1,cmd2,sep=" ")

  print(cmd3)

  system(cmd3,intern= TRUE)

  cat("Finished qcBAM\n")

}


bam2Bw <- function(R_lib)
{
  cat("You choose to convert bam files to bw files, please define the following setting parameters: \n",
      "input.bam.files.dir (ex: /scratch/projects/exempt/tr/azhang/for_bioinfo_core/3_end_seq)\n",
      "input.chrosome.size.file (ex: /nethome/axy148/hg19.genome)\n",
      "output.bw.file.dir (ex:BW)\n")

  input <- file("stdin", "r")
  count.file.dir <- readLines(input, n = 3)

  input.bam.file.dir <- count.file.dir[1]
  input.ref.gene.bed.file <- count.file.dir[2]
  output.dir <- count.file.dir[3]

  library(ChipSeq)
  library(ThreeUTR)


  cmd1 = "bsub -P bbc -J \"QCBam\" -o %J.QCBam.log -e %J.QCBam.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

  cmd2 = paste("Rscript",paste0(R_lib,"/ThreeUTR/bin/rscript/bam2Bw.r"),
               input.bam.file.dir,input.ref.gene.bed.file,output.dir,sep=" ")

  cmd3 = paste(cmd1,cmd2,sep=" ")

  print(cmd3)

  system(cmd3,intern= TRUE)

  cat("Finished bam2Bw\n")

}

analysisAll <- function(R_lib)
{
  cat("Do you want to perform peak calling, annotation, and coverage visualization?\n")

  input <- file("stdin", "r")
  row <- readLines(input, n = 1)

  print(row)

  if (row == "Yes")
  {



    cat("Do you have input bam file for control ?\n")

    input <- file("stdin", "r")
    row <- readLines(input, n = 1)

    if (row == "No")
    {

      cat("please define the input Bam file directory:\n")

      input <- file("stdin", "r")
      input.file.dir <- readLines(input, n = 1)

      # input.file.dir=row

      cat("please define the output file directory:\n")

      input <- file("stdin", "r")
      output.file.dir <- readLines(input, n = 1)

      # out.file.dir=row

      cat("please define genome name:\n")

      input <- file("stdin", "r")
      genome <- readLines(input, n = 1)

      cmd1 = "bsub -P bbc -J \"RunR\" -o %J.RunR.log -e %J.RunR.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
      cmd2 = paste0("Rscript ", R_lib, "/ChipSeq/bin/Run_Chip_Seq.r ",
                    input.file.dir, " ", output.file.dir, " ", genome)

      system(paste0(cmd1, " ", cmd2))

      # cmd2=paste0()

      # system(cmd2)

      cat("Finished peak calling, annotation, and coverage visualizatio...\n")

    } else if (row == "Yes")
    {

      cat("please define the sample information file:\n")

      input <- file("stdin", "r")
      sample.info.file <- readLines(input, n = 1)

      # input.file.dir=row

      cat("please define the Bam files information file:\n")

      input <- file("stdin", "r")
      bam.info.file <- readLines(input, n = 1)

      # out.file.dir=row

      cat("please define genome name:\n")

      input <- file("stdin", "r")
      genome <- readLines(input, n = 1)

      cat("Which peakcaller you want to use, please choose: macs14 or macs2 \n")

      input <- file("stdin", "r")
      peakcaller <- readLines(input, n = 1)

      cat("The peakcaller you choose is :", peakcaller,
          "please define p value threshold for peak calling: \n")

      input <- file("stdin", "r")
      peakPvalue <- readLines(input, n = 1)

      cat("Please define the name of output directory for peaks: \n")

      input <- file("stdin", "r")
      peak.out.dir <- readLines(input, n = 1)

      peak.out.dir <- paste0(peak.out.dir, "_", genome,
                             "_", peakcaller, "_", peakPvalue)

      cmd1 = "bsub -P bbc -J \"InputCall\" -o %J.InputCall.log -e %J.InputCall.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
      cmd2 = paste("Rscript", file.path(R_lib, "ChipSeq/bin/UseInputCall.r"),
                   sample.info.file, bam.info.file, genome, peak.out.dir,
                   peakcaller, peakPvalue, sep = " ")

      system(paste0(cmd1, " ", cmd2))

    }

  } else if (row == "No")
  {

    cat("Do you want to perform peak merge, overlap, produing Venn diagram?\n")

    input <- file("stdin", "r")
    row <- readLines(input, n = 1)

    if (row == "Yes")
    {

      cat("please defne the input file for sample information:\n")

      input <- file("stdin", "r")
      input.file.sample <- readLines(input, n = 1)

      cat("Finished to get sample information\n")


      cat("please defne the input file directory for bed files:\n")

      input <- file("stdin", "r")
      input.file.dir <- readLines(input, n = 1)

      cat("please defne the input file pattern:\n")
      input <- file("stdin", "r")
      input.file.pattern <- readLines(input, n = 1)

      cat("please defne the output file directory:\n")
      input <- file("stdin", "r")
      output.file.dir <- readLines(input, n = 1)

      R_lib = .libPaths()[1]

      cmd1 = "bsub -P bbc -J \"RunR\" -o %J.RunR.log -e %J.RunR.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
      cmd2 = paste0("Rscript ", R_lib, "/ChipSeq/bin/Run_Merge_peak.r ",
                    input.file.dir, " ", output.file.dir, " ", input.file.pattern,
                    " ", input.file.sample)

      system(paste0(cmd1, " ", cmd2))

      cat("Finished peak merge, overlap, produing Venn diagram?\n")
    } else if (row == "No")
    {

      cat("Do you want to perform  Bam files sorting, index and Visualization?\n")

      input <- file("stdin", "r")
      row <- readLines(input, n = 1)

      if (row == "Yes")
      {

        cat("please defne the input file directory for Bam files:\n")

        input <- file("stdin", "r")
        input.file.dir <- readLines(input, n = 1)

        cat("Finished to get bam files\n")

        cat("please defne the out file directory:\n")

        input <- file("stdin", "r")
        output.file.dir <- readLines(input, n = 1)

        cat("please defne genome:(example: Hs) \n")


        input <- file("stdin", "r")
        genomeID <- readLines(input, n = 1)

        cmd1 = "bsub -P bbc -J \"BamR\" -o %J.BamR.log -e %J.BamR.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

        cmd2 = paste0("Rscript", R_lib,"/ChipSeq/bin/Process_Bam.r ",
                      input.file.dir, " ", output.file.dir, " ",
                      genomeID)

        system(paste0(cmd1, " ", cmd2))

        cat("Finished Bam files sorting, index and Visualization ?\n")

      } else
      {
        quit()
      }
    }
  }
}

switch(choose.type,

Download_SRA ={
  downloadSRA(R_lib)
},
SRA2Fastq ={
  sra2Fastq(R_lib)
},
subset_Fastq_file ={
  subsetFastqFile(R_lib)
},
testAlignment ={
  testAlignment(R_lib)
},
Alignment ={
performAlignment(R_lib)
},
QC = {
  bamQC(R_lib)
},
Bam2Bw ={
  bam2Bw(R_lib)
},
peakcalling = {
},
annotation = {

},
visualization = {
  analysisVisualization(R_lib)
}
)
