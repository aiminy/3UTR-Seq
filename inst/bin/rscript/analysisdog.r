args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an
# error

if (length(args) == 0)
{
    stop("At least one argument must be supplied (input file).n",
        call. = FALSE)
} else if (length(args) > 1)
{
    input.count.file.dir = args[1]
    input.file.pattern = args[2]
    output.anlysis.dir = args[3]
    out.file.pattern.interested = args[4]
    out.file.pattern.positive.gene = args[5]
    out.file.pattern.negative.gene = args[6]
    out.file.pattern.all = args[7]
    dir.name.gene.list = args[8]
    pattern.4.gene.list = args[9]
    adjust_by_batch = args[10]
}

if (!dir.exists(output.anlysis.dir))
{
    dir.create(output.anlysis.dir)
}

library(ChipSeq)
library(ThreeUTR)
library(org.Hs.eg.db)

Re.unadjusted.adjusted <- ProcessOutputFilesFromDoGsOnly(input.count.file.dir,
    input.file.pattern, output.anlysis.dir, out.file.pattern.interested,
    out.file.pattern.positive.gene, out.file.pattern.negative.gene,
    out.file.pattern.all, dir.name.gene.list, pattern.4.gene.list,
    adjust_by_batch)
