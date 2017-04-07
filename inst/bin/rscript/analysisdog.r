args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an
# error

if (length(args) == 0) {
    stop("At least one argument must be supplied (input file).n",
        call. = FALSE)
} else if (length(args) > 1) {
    input.count.file.dir = arg[1]
    input.file.pattern = arg[2]
    output.anlysis.dir = arg[3]
    out.file.pattern.interested = arg[4]
    out.file.pattern.positive.gene = arg[5]
    out.file.pattern.negative.gene = arg[6]
    out.file.pattern.all = arg[7]
    dir.name.gene.list = arg[8]
    pattern.4.gene.list = arg[9]
    adjust_by_batch = arg[10]
}

if (!dir.exists(output.count.file.dir)) {
    dir.create(output.count.file.dir)
}

library(ChipSeq)
library(ThreeUTR)

Re.unadjusted.adjusted <- ProcessOutputFilesFromDoGsOnly(input.count.file.dir,
    input.file.pattern, output.anlysis.dir, out.file.pattern.interested,
    out.file.pattern.positive.gene, out.file.pattern.negative.gene,
    out.file.pattern.all, dir.name.gene.list, pattern.4.gene.list,
    adjust_by_batch = "YES")
