#!/usr/bin/Rscript

# Usage:

# Rscript ~/R/lib64/R/library/ThreeUTR/bin/rscript/GetCountUseBash.r

R_lib = .libPaths()[1]

cat("Do you want to get count?\n")

input <- file("stdin", "r")
row <- readLines(input, n = 1)

print(row)

if (row == "Yes") {

  cat("please define the input Bam file directory:\n")

  input <- file("stdin", "r")
  input.bamfile.dir <- readLines(file("stdin", "r"), n = 1)

  cat("please define the output file directory:\n")

  input <- file("stdin", "r")
  output.count.file.dir <- readLines(input, n = 1)

  cat("please define the annotation file:\n")

  input <- file("stdin", "r")
  annotation.bed.file <- readLines(input, n = 1)

  cat("please define the upstream base of each transcript :\n")

  input <- file("stdin", "r")
  ld <- readLines(input, n = 1)

  cat("please define the downstream base of each transcript :\n")

  input <- file("stdin", "r")
  rd <- readLines(input, n = 1)

  cat("please define the sample that will not used for count :\n")

  input <- file("stdin", "r")
  no.use.sample <- readLines(input, n = 1)

  cmd1 = "bsub -P bbc -J \"DogFT\" -o %J.DogFT.log -e %J.DogFT.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
  cmd2 = paste0(R_lib, "/ThreeUTR/bin/rscript/convertbam2bed.r")

  cmd3 = paste("Rscript", cmd2, input.bamfile.dir, annotation.bed.file,
    ld, rd, output.count.file.dir,no.use.sample,sep = " ")

  system(paste0(cmd1, " ", cmd3))

  cat("Finished conversion...\n")

} else {

   cat("I already have counts, so I need to perform analysis.\n")

  library("optparse")

  option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
                help="output file name [default= %default]", metavar="character")
  );

  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);

  if (is.null(opt$file)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  }

  #  cat("Please define the following setting parameters: \n")
  #
  #
  #  input.count.file.dir = arg[1]
  #  input.file.pattern = arg[2]
  #  output.anlysis.dir = arg[3]
  #  out.file.pattern.interested = arg[4]
  #  out.file.pattern.positive.gene = arg[5]
  #  out.file.pattern.negative.gene = arg[6]
  #  out.file.pattern.all = arg[7]
  #  dir.name.gene.list = arg[8]
  #  pattern.4.gene.list = arg[9]
  #  adjust_by_batch = arg[10]
  #
  #
  #  input <- file("stdin", "r")
  #  count.file.dir <- readLines(input, n = 1)
  #
  # cmd1 = "bsub -P bbc -J \"DogFT\" -o %J.DogFT.log -e %J.DogFT.err -W 72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
  # cmd2 = paste0(R_lib, "/ThreeUTR/bin/rscript/analysis.r")
  #
  # cmd3 = paste("Rscript", cmd2, input.bamfile.dir, annotation.bed.file,
  #              ld, rd, output.count.file.dir,no.use.sample,sep = " ")
  #
  # system(paste0(cmd1, " ", cmd3))
  #
  # cat("Finished analysis...\n")

  #quit()
  }
