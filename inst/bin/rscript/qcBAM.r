
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {

  input.bamfile.dir=args[1]
  input.ref.gene.bed.file=args[2]
  output.dir=args[3]
}

#cat(input.bamfile.dir,"\t",output.bedfile.dir,"\n")

library(ChipSeq)
library(ThreeUTR)

ThreeUTR:::useInferExperiment(input.bamfile.dir,input.ref.gene.bed.file,output.dir)

