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

parsersample <- function()
{
    cell <- factor(rep(c("batch1", "batch2"), c(3, 2)))
    cell = rep(cell, 2)

    sample <- c("R1_Dox.bam", "R2_Dox.bam", "R3_Dox.bam", "R4_Dox.bam", "R5_Dox.bam",
        "R1_WT.bam", "R2_WT.bam", "R3_WT.bam", "R4_WT.bam", "R5_WT.bam")

    colData <- data.frame(sample = sample, condition = factor(rep(c("Dox", "WT"),
        c(5, 5))), cell = cell)

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
#' res <- parserreadfiles(input.file.dir,'bam',filter.sample='emp1')
#'
parserreadfiles <- function(input.file.dir, input.file.type, sample.group = NULL,
    filter.sample = NULL)
    {
    dir.name = input.file.dir
    dir.name = reformatPath(dir.name)
    file.name = file.path(dir.name, dir(dir.name, recursive = TRUE))
    file.name.2 <- as.list(file.name)
    file.name.3 <- lapply(1:length(file.name.2), function(u, input.file.type,
        file.name.2)
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

        file_name <- paste0(path_name2, "-", file_name)

        file_name

    }, file.name.4))

    if (!is.null(filter.sample))
    {
        file.name.5 <- file.name.4[-grep(filter.sample, names(file.name.4))]

    } else
    {
        file.name.5 <- file.name.4
    }


    if (!is.null(sample.group))
    {
        g1 <- grep(sample.group[1], toupper(names(res$input)))

        g2 <- grep(sample.group[2], toupper(names(res$input)))

        # output.dir.name = reformatPath(output.file.dir) temp3 = output.dir.name
        re2 <- list(input = file.name.5, g1 = g1, g2 = g2, input.file.type = input.file.type)
    } else
    {
        re2 <- list(input = file.name.5, input.file.type = input.file.type)
    }

    pkg.env <- new.env(parent = emptyenv())
    pkg.env$sample <- ThreeUTR:::parsersample()

    return(re2)


}

useWget2Download <- function(sra.accession.number, output.dir)
{
    cmd0 <- "wget -c -r -nd -np -L"
    cmd1 <- "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/"

    temp <- file.path(substr(sra.accession.number, 1, 6), sra.accession.number)

    temp2 <- paste0(cmd1, temp)

    cmd2 <- paste(cmd0, temp2, "-P", output.dir, sep = " ")

    system(cmd2)

}

useFastqDump <- function(sra.accession.number, output.dir)
{
    cmd0 <- "fastq-dump -I --split-files"

    cmd1 <- paste(cmd0, sra.accession.number, "-O", output.dir, sep = " ")

    system(cmd1)

}

useFastqDumpConvertSra2Fastq <- function(sra.file.dir, output.dir)
{
    cmd0 <- "fastq-dump --split-3"

    re <- parserreadfiles(sra.file.dir, "sra")

    res <- re$input

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir)
    }

    cmd.l <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        cmd1 <- paste(cmd0, u, "-O", output.dir, sep = " ")

        system(cmd1)

        cmd1
    }, output.dir)

    re <- list(cmdl = cmd.l, output.dir = output.dir)

    re

}

# infer_experiment.py
useInferExperiment <- function(input.file.dir, ref.gene.bed.file, output.dir)
{
    re <- parserreadfiles(input.file.dir, "bam")

    res <- re$input

    cmd0 <- "infer_experiment.py -i"
    cmd1 <- "-r"
    cmd2 <- ">"

    # output.dir <- file.path(output.dir, 'BamInfo')

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    cmd.l <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd2 <- paste(cmd0, u, cmd1, ref.gene.bed.file, cmd2, file.path(output.dir,
            paste0(file_name, "_infor.txt")), sep = " ")

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
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertbam2bed("/scratch/projects/bbc/aiminy_project/DoGs/BAM","/scratch/projects/bbc/aiminy_project/DoGs")'
#'
convertbam2bed <- function(input.bamfile.dir, output.bedfile.dir)
{
    res <- parserreadfiles(input.bamfile.dir, "bam")

    res <- res$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    #if (!dir.exists(output.bedfile.dir))
    #{
    #  dir.create(output.bedfile.dir, recursive = TRUE)
    #}

    output.bedfile.dir <- file.path(output.bedfile.dir, "BedFileFromBam")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id,res,output.bedfile.dir)
    {
        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11
        # <- '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13
        # <- paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)

        path_name = dirname(res[[u]])
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        file_name <- paste0(path_name2, "-", file_name)
        if (m.id == 1)
        {
        job.name <- paste0("bam2bed.",u)
        cmd0 <- ChipSeq:::usePegasus('parallel', Wall.time = '72:00',cores = 32,Memory = 25000,span.ptile = 16,job.name)

        cmd1 <- "bedtools bamtobed -i"
        cmd2 <- "\\>"

        cmd3 <- paste(cmd0,cmd1, res[[u]], cmd2, file.path(output.bedfile.dir, paste0(file_name,
            ".bed")), sep = " ")
        }else
        {
        cmd3 <- paste(cmd1, res[[u]], cmd2, file.path(output.bedfile.dir, paste0(file_name,
                                                                                 ".bed")), sep = " ")
        }

        system(cmd3)

        cmd3
    }, m.id,res,output.bedfile.dir)

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
matchbed2annotation <- function(input.bedfile.dir, annotation.bed.file, ld,
    rd, output.matched.bed.file.dir)
    {
    res <- parserreadfiles(input.bedfile.dir, "bed")

    res <- res$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    output.bedfile.dir <- file.path(output.matched.bed.file.dir, "MatchedBedFile")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u,res,m.id,ld,rd,annotation.bed.file,output.bedfile.dir)
    {
        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11
        # <- '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13
        # <- paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)
        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
          job.name <- paste0("bedMannot.",u)
          wait.job.name <- paste0("bam2bed.",u)
          cmd.p <- ChipSeq:::usePegasus('parallel', Wall.time = '72:00',cores = 32,Memory = 25000,span.ptile = 16,job.name,wait.job.name)

          #cmd1 <- "bedtools bamtobed -i"
          #cmd2 <- "\\>"

          cmd0 <- paste("bedtools window -a", annotation.bed.file, "-b", sep = " ")
          cmd1 <- paste("-l", ld, "-r", rd, "-sw", "\\>", sep = " ")
          cmd2 <- paste(cmd.p,cmd0, res[[u]], cmd1, file.path(output.bedfile.dir, paste0(file_name,
                                                                            "_matched.bed")), sep = " ")

        }else
        {
          cmd0 <- paste("bedtools window -a", annotation.bed.file, "-b", sep = " ")
          cmd1 <- paste("-l", ld, "-r", rd, "-sw", ">", sep = " ")
          cmd2 <- paste(cmd.p,cmd0, res[[u]], cmd1, file.path(output.bedfile.dir, paste0(file_name,
                                                                                         "_matched.bed")), sep= " ")
        }

        system(cmd2)

        cmd2
    }, m.id,ld,rd,annotation.bed.file,output.bedfile.dir)

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
getcountsfromMatchedbed <- function(input.bedfile.dir, output.count.file.dir,
    filter.sample)
    {

    res <- parserreadfiles(input.bedfile.dir, "bed", filter.sample = filter.sample)

    res <- res$input

    # system('awk -F '\\t' '$6==\'+\' && $12==\'-\''
    # ~/MatchedBedFile/R1_Dox_matched.bed | awk '$8<$2&&$9>=$2' | awk '{print
    # $4}' | sort | uniq -c | sort -nr | head')

    cmd0 <- "awk -F '\\t'"

    cmd1 <- "'$6==\"+\" && $12==\"-\"'"
    cmd11 <- "'$6==\"+\" && $12==\"+\"'"
    cmd12 <- "'$6==\"-\" && $12==\"+\"'"
    cmd13 <- "'$6==\"-\" && $12==\"-\"'"

    cmd2 <- "| awk '$8 < $2 && $9 >= $2' | awk '{print $4}' | sort | uniq -c | sort -nr"  #below
    cmd21 <- "| awk '$8 >= $2 && $9 <= $3' | awk '{print $4}' | sort | uniq -c | sort -nr"  #DoGs
    cmd22 <- "| awk '$8 <= $3 && $9 > $3' | awk '{print $4}' | sort | uniq -c | sort -nr"  #over

    cmd3 <- ">"

    cmdtemp <- rbind(cbind(rep(cmd1, 3), c(cmd2, cmd21, cmd22), rep("plus",
        3), rep("minus", 3), c("below.DoGs", "DoGs", "over.DoGs")), cbind(rep(cmd11,
        3), c(cmd2, cmd21, cmd22), rep("plus", 3), rep("plus", 3), c("below.DoGs",
        "DoGs", "over.DoGs")), cbind(rep(cmd12, 3), c(cmd2, cmd21, cmd22), rep("minus",
        3), rep("plus", 3), c("over.DoGs", "DoGs", "below.DoGs")), cbind(rep(cmd13,
        3), c(cmd2, cmd21, cmd22), rep("minus", 3), rep("minus", 3), c("over.DoGs",
        "DoGs", "below.DoGs")))

    output.count.file.dir <- file.path(output.count.file.dir, "Counts")

    if (!dir.exists(output.count.file.dir))
    {
        dir.create(output.count.file.dir, recursive = TRUE)
    }

    counteachcase <- function(res, cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand,
        location, output.count.file.dir)
        {
        cmd.l <- lapply(1:length(res), function(u,res,cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand,
                                      location,output.count.file.dir)
        {
            file_name = file_path_sans_ext(basename(res[[u]]))

            if (m.id == 1)
            {
              job.name <- paste0("count.",u)
              wait.job.name <- paste0("bedMannot.",u)

              cmd.p <- ChipSeq:::usePegasus('parallel', Wall.time = '72:00',cores = 32,Memory = 25000,span.ptile = 16,job.name,wait.job.name)
              cmd3 <- "\\>"

            cmd4 <- paste(cmd.p,cmd0, cmd1, res[[u]], cmd2, cmd3, file.path(output.count.file.dir,
                paste0(file_name, ".", gene.strand, ".gene.", read.strand, ".read.",
                  location, ".count.txt")), sep = " ")

            }else
            {
              cmd4 <- paste(cmd0, cmd1, res[[u]], cmd2, cmd3, file.path(output.count.file.dir,                                                  paste0(file_name, ".", gene.strand, ".gene.", read.strand, ".read.",
                   location, ".count.txt")), sep = " ")

            }

            cat(cmd4, "\n")

            system(cmd4)

            cmd4

        },res,cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand,
        location,output.count.file.dir)

        return(cmd.l)

    }

    cmdtempres2 <- apply(cmdtemp, 1, function(u, cmd0, cmd3, output.count.file.dir)
    {
        x <- as.data.frame(t(u))

        cmd1 <- x[, 1]
        cmd2 <- x[, 2]
        gene.strand <- x[, 3]

        read.strand <- x[, 4]

        location <- x[, 5]

        cmdtempres <- counteachcase(res, cmd0, cmd1, cmd2, cmd3, gene.strand,
            read.strand, location, output.count.file.dir)

        cmdtempres

    }, cmd0, cmd3, output.count.file.dir)

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
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::getcounts("/scratch/projects/bbc/aiminy_project/DoGs/BAM","/projects/ctsi/bbc/aimin/annotation/hg19_DoGs_2.bed",0,0,"/scratch/projects/bbc/aiminy_project/DoGs")'
#'
getcounts <- function(input.bamfile.dir, annotation.bed.file, ld, rd, output.count.file.dir,
    filter.sample)
    {
    res <- convertbam2bed(input.bamfile.dir, output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    annotation.bed.file <- annotation.bed.file
    ld <- ld
    rd <- rd

    res <- matchbed2annotation(input.bedfile.dir, annotation.bed.file, ld, rd,
        output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    res <- getcountsfromMatchedbed(input.bedfile.dir, output.count.file.dir,
        filter.sample)

    res
}

#' parserAnnotationFile
#'
#' @return
#' @export
#'
#' @examples
#'
#' input.annotation.file <- '/Volumes/Bioinformatics$/Aimin_project/UTR/3UTR_GFF2.txt'
#'
#' y <- parserAnnotationFile(input.annotation.file)
#'
#'
parserAnnotationFile <- function(input.annotation.file)
{
    dir.name <- dirname(input.annotation.file)

    file.name <- file_path_sans_ext(basename(input.annotation.file))

    x <- read.table(input.annotation.file)

    d <- x[, 5] - x[, 4]

    nn <- substr(x[, 9], 9, 23)

    xx <- cbind.data.frame(x, nn, d)

    counts <- as.data.frame(table(xx$nn))

    colnames(counts) = c("nn", "counts")

    xxx <- merge(xx, counts, by = "nn", sort = FALSE)

    xxxx <- xxx[, c(2, 5, 6, 10, 11, 8)]

    write.table(xxxx, file = file.path(dir.name, paste0(file.name, ".bed")),
        row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

    xxxx
}

convertBam2Bw <- function(input.bam.file.dir, input.chromosome.size.file, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    cmd0 = "samtools sort"

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd <- paste(cmd0, u, file.path(output.bw.file.dir, paste0(file_name,
            "_sorted")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd1 = "samtools index"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd <- paste(cmd1, file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam")),
            sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd3 <- "tail -n +2"

    path_name = dirname(input.chromosome.size.file)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(input.chromosome.size.file))

    file_name <- paste0(path_name2, "-", file_name, ".txt")

    cmd <- paste(cmd3, input.chromosome.size.file, ">", file.path(path_name,
        file_name), sep = " ")

    system(cmd)

    cmd4 <- "genomeCoverageBed -ibam"
    cmd5 <- "-bg -g"
    cmd6 <- ">"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam"))

        cmd <- paste(cmd4, x, cmd5, input.chromosome.size.file, cmd6, file.path(output.bw.file.dir,
            paste0(file_name, ".bdg")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    # re <- list(cmdl = cmd.l, output.dir = output.dir)

    # input.bdg.file.dir <- re$output.dir

    cmd7 <- "LC_COLLATE=C sort -k1,1 -k2,2n"

    # re <- parserreadfiles(input.bdg.file.dir,'bdg')

    # res <- re$input

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".bdg"))

        cmd <- paste(cmd7, x, cmd6, file.path(output.bw.file.dir, paste0(file_name,
            ".sorted_bdg")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd8 = "bedGraphToBigWig"

    input.chromosome.size.file.m <- file.path(path_name, file_name)

    # re <- parserreadfiles(input.bdg.file.dir,'.sorted_bdg')

    # res <- re$input

    cmd.l <- lapply(res, function(u, input.chromosome.size.file.m, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".sorted_bdg"))

        cmd <- paste(cmd8, x, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, ".bw")), sep = " ")

        system(cmd)

        cmd
    }, input.chromosome.size.file.m, output.bw.file.dir)

    cmd9 = "bigWigToWig"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".bw"))

        cmd <- paste(cmd9, x, file.path(output.bw.file.dir, paste0(file_name,
            ".wig")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

}

#' @examples
#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2StrandBw('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/nethome/axy148/hg19.genome','/scratch/projects/bbc/aiminy_project/DoGs/BW2')'
#'
convertBam2StrandBw <- function(input.bam.file.dir, input.chromosome.size.file,
    output.bw.file.dir)
    {
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    # job.name=paste0('bamSort[',length(res),']')

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamSort.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))

            # job.name=paste0('bamSort[',length(res),']') cmd1 = paste0('bsub -w
            # \'done(\'bamSort[*]\')\'', 'bsub -P bbc -J \'',job.name,paste0('\'
            # -o %J.log '),paste0('-e %J.err -W')) job.name=paste0('Bdg[',u,']') cmd1 =
            # paste0('bsub -w \'done(\'bamIndex[*]\') && done(\'Chrosome\')\'',
            # 'bsub -P bbc -J \'',job.name,paste0('\' -o %J.',job.name,'.log
            # '),paste0('-e %J.',job.name,'.err -W'))

            cmd2 = "samtools sort"
            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            cmd3 = "samtools sort"
        }

        cmd <- paste(cmd3, res[[u]], file.path(output.bw.file.dir, paste0(file_name,
            "_sorted")), sep = " ")

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamIndex.", u)
            wait.job.name = paste0("bamSort.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            # cmd1 = 'bsub -w \'done(\'bamSort\')\' -P bbc -J \'bamIndex\' -o
            # %J.bamIndex.log -e %J.bamIndex.err -W'
            cmd4 = "samtools index"
            cmd5 = paste(cmd1, cmd0, cmd4, sep = " ")
        } else
        {
            cmd5 = "samtools index"
        }

        cmd <- paste(cmd5, file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam")),
            sep = " ")

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

    path_name <- dirname(input.chromosome.size.file)
    path_name2 <- basename(path_name)
    file_name <- file_path_sans_ext(basename(input.chromosome.size.file))
    file_name <- paste0(path_name2, "-", file_name, ".txt")

    if (m.id == 1)
    {
        cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
        cmd1 = "bsub -P bbc -J \"Chrosome\" -o %J.Chrosome.log -e %J.Chrosome.err -W"
        cmd6 <- "tail -n +2"
        cmd7 <- paste(cmd1, cmd0, cmd6, sep = " ")
        cmd8 <- paste(cmd7, input.chromosome.size.file, "\\>", file.path(path_name,
            file_name), sep = " ")
    } else
    {
        cmd7 <- "tail -n +2"
        cmd8 <- paste(cmd7, input.chromosome.size.file, ">", file.path(path_name,
            file_name), sep = " ")
    }

    system(cmd8)

    # cmd11='bsub -w \'done(\'STAR-alignment\')\' -P bbc -J
    # \'samtools-sort\' -o %J.samtools-sort.log -e %J.samtools-sort.err -W'

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

            # u=1
            job.name = paste0("Bdg.", u)
            wait.job.name = paste(paste0("bamIndex.", u, "\")"), "&& done(\"Chrosome\")",
                sep = " ")
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            cmd9 = "genomeCoverageBed -split -strand + -ibam"
            cmd10 = "genomeCoverageBed -split -strand - -ibam"
            cmd11 = paste(cmd1, cmd0, cmd9, sep = " ")
            cmd12 = paste(cmd1, cmd0, cmd10, sep = " ")
            cmd13 <- "\\>"
        } else
        {
            cmd11 = "genomeCoverageBed -split -strand + -ibam"
            cmd12 = "genomeCoverageBed -split -strand - -ibam"
            cmd13 <- ">"
        }

        cmd14 <- "-bg -g"

        cmd.x <- paste(cmd11, x, cmd14, input.chromosome.size.file, cmd13, file.path(output.bw.file.dir,
            paste0(file_name, "_plus.bdg")), sep = " ")

        cmd.y <- paste(cmd12, x, cmd14, input.chromosome.size.file, cmd13, file.path(output.bw.file.dir,
            paste0(file_name, "_minus.bdg")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd

    }, m.id, res, output.bw.file.dir)


    # re <- list(cmdl = cmd.l, output.dir = output.dir)

    # input.bdg.file.dir <- re$output.dir

    # re <- parserreadfiles(input.bdg.file.dir,'bdg')

    # res <- re$input

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.bdg"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.bdg"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("sortBdg.", u)
            wait.job.name = paste(paste0("Bdg.", u, "\")"), "&& done(\"Chrosome\")",
                sep = " ")
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            cmd15 <- "LC_COLLATE=C sort -k1,1 -k2,2n"
            cmd16 <- paste(cmd1, cmd0, cmd15, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd16 <- "LC_COLLATE=C sort -k1,1 -k2,2n"
            cmd17 <- ">"
        }

        cmd.x <- paste(cmd16, x, cmd17, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.sorted_bdg")), sep = " ")
        cmd.y <- paste(cmd16, y, cmd17, file.path(output.bw.file.dir, paste0(file_name,
            "_minus.sorted_bdg")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, output.bw.file.dir)


    input.chromosome.size.file.m <- file.path(path_name, file_name)

    # re <- parserreadfiles(input.bdg.file.dir,'.sorted_bdg')

    # res <- re$input

    cmd.l <- lapply(1:length(res), function(u, m.id, res, input.chromosome.size.file.m,
        output.bw.file.dir)
        {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.sorted_bdg"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.sorted_bdg"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

            job.name = paste0("BigWig.", u)
            wait.job.name = paste0("sortBdg.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))
            cmd18 = "bedGraphToBigWig"
            cmd19 <- paste(cmd1, cmd0, cmd18, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd19 <- "bedGraphToBigWig"
        }

        cmd.x <- paste(cmd19, x, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, "_plus.bw")), sep = " ")

        cmd.y <- paste(cmd19, y, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, "_minus.bw")), sep = " ")


        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, input.chromosome.size.file.m, output.bw.file.dir)


    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.bw"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.bw"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"


            job.name = paste0("Wig.", u)
            wait.job.name = paste0("BigWig.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))
            cmd20 = "bigWigToWig"
            cmd21 <- paste(cmd1, cmd0, cmd20, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd21 <- "bigWigToWig"
        }

        cmd.x <- paste(cmd21, x, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.wig")), sep = " ")

        cmd.y <- paste(cmd21, y, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.wig")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, output.bw.file.dir)

}

# HCI-Scripts/BamFile/split_bam_by_isize.pl
#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::splitBam('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs/Bam_split')'
#'
splitBam <- function(input.bam.file.dir, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    # job.name=paste0('bamSort[',length(res),']')

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)
        u <- 3
        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamSplit.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))

            # job.name=paste0('bamSort[',length(res),']') cmd1 = paste0('bsub -w
            # \'done(\'bamSort[*]\')\'', 'bsub -P bbc -J \'',job.name,paste0('\'
            # -o %J.log '),paste0('-e %J.err -W')) job.name=paste0('Bdg[',u,']') cmd1 =
            # paste0('bsub -w \'done(\'bamIndex[*]\') && done(\'Chrosome\')\'',
            # 'bsub -P bbc -J \'',job.name,paste0('\' -o %J.',job.name,'.log
            # '),paste0('-e %J.',job.name,'.err -W'))
            if (u <= 6)
            {
                cmd2 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd2 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            }
            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            if (u <= 6)
            {
                cmd3 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd3 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            }
        }

        cmd <- cmd3

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

}

#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2StrandBw2("/scratch/projects/bbc/aiminy_project/DoGs/Bam_split/","/scratch/projects/bbc/aiminy_project/DoGs/BW_split")'
#'

convertBam2StrandBw2 <- function(input.bam.file.dir, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
        span.ptile, res, output.bw.file.dir)
        {

        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            job.name = paste0("bam2wig.", u)
            cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 16000, span.ptile = 16, job.name)

            if (u <= 6)
            {
                cmd2 = paste("bam2wig.pl -pe --pos span --strand --bw  --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd2 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            }
            cmd3 = paste(cmd1, cmd2, sep = " ")
        } else
        {
            if (u <= 6)
            {
                cmd3 = paste("bam2wig.pl -pe --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd3 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            }
        }

        cmd <- cmd3

        system(cmd)

        cmd
    }, m.id, Wall.time, cores, Memory, span.ptile, res, output.bw.file.dir)

}

prepareDaPars <- function(input.wig.file.dir, sample.group = c("Dox", "WT"),
    output.dir)
    {
    re <- parserreadfiles(input.wig.file.dir, "wig")

    res <- re$input

    cmd = "tail -n +2"

    cmd.l <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)
        # path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = u

        cmd <- paste(cmd, x, ">", file.path(output.dir, paste0(file_name, "2.wig")),
            sep = " ")

        system(cmd)

        cmd
    }, output.dir)

    g1 <- "1,2,3,4,5"
    g2 <- "a,b,c,d,e"

    cat(c("Annotated_3UTR=hg19_refseq_extracted_3UTR.bed\n", paste0("Group1_Tophat_aligned_Wig=",
        g1, "\n"), paste0("Group2_Tophat_aligned_Wig=", g2, "\n"), "Output_directory=DaPars_Dox_WT_data/\n",
        "Output_result_file=DaPars_Dox_WT_data\n", "#Parameters\n
        Num_least_in_group1=1\n
        Num_least_in_group2=1\n
        Coverage_cutoff=30\n
        FDR_cutoff=0.05\n
        PDUI_cutoff=0.5\n
        Fold_change_cutoff=0.59\n"),
        file = file.path(output.dir, "config_test.txt"), sep = "")
}

select3UTR <- function(genome, tablename)
{
    library(GenomicFeatures)
    library(dplyr)
    refSeq <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
    threeUTRs <- threeUTRsByTranscript(refSeq, use.names = TRUE)
    length_threeUTRs <- width(ranges(threeUTRs))
    the_lengths <- as.data.frame(length_threeUTRs)
    the_lengths <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
    the_lengths <- unique(the_lengths[, c("group_name", "sum(value)")])
    colnames(the_lengths) <- c("RefSeq Transcript", "3' UTR Length")

}

subsetFastq <- function(input.fastq.files.dir, output.dir, n)
{
    re <- parserreadfiles(input.fastq.files.dir, "fastq")
    res <- re$input

    # cmd0 = paste0('head -',n)
    cmd0 = "seqtk sample -s100"
    # read1.fq 10000

    cmd1 = "\\>"
    # cmd1 = '>'

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd3 = "bsub -P bbc -J \"subSetFastq\" -o %J.subSetFastq.log -e %J.subSetFastq.err -W"

    # cmd4='bsub -w \'done(\'subSetFastq\')\' -P bbc -J \'tophat\' -o
    # %J.tophat.log -e %J.tophat.err -W'

    xx <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        sample.name.out = file.path(output.dir, paste0(file_name, "-test-",
            n, ".fastq"))

        cmd = paste(cmd3, cmd2, cmd0, u, n, cmd1, sample.name.out, sep = " ")

        print(cmd)

        cmd

        system(cmd, intern = TRUE)

    }, output.dir)

}

checkStrand <- function(input.alignment.dir)
{
    re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
        pattern = "junctions.bed"))

    cmd0 <- "wc -l"

    y <- lapply(re, function(u)
    {
        cmd <- paste(cmd0, u, sep = " ")

        # system(cmd,intern = TRUE)
        cmd
    })

    yyy <- lapply(y, function(u)
    {
        yy <- system(u, intern = TRUE)
        yy
    })

    yyyy <- unlist(yyy)

    return(yyyy)
}

#' @examples
#' R -e 'library(ThreeUTR);ThreeUTR:::processBamFiles('/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2','/scratch/projects/bbc/aiminy_project/DoGs/BAM')'
#'
processBamFiles <- function(input.alignment.dir, output.dir)
{
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
        pattern = "accepted_hits.bam"))

    cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd1 = "bsub -P bbc -J \"bamProcess\" -o %J.bamProcess.log -e %J.bamProcess.err -W"

    cmd2 <- "mv"

    y <- lapply(re, function(u)
    {
        # /SRR2038198/Fs12/accepted_hits.bam

        dir.name.0 <- dirname(u)
        dir.name.1 <- dirname(dir.name.0)

        x <- basename(dir.name.0)
        y <- basename(dir.name.1)
        file.name <- basename(u)

        sample.name <- paste(y, x, file.name, sep = "-")

        cmd <- paste(cmd1, cmd0, cmd2, u, file.path(output.dir, sample.name),
            sep = " ")

        # system(cmd,intern = TRUE
        cmd
    })

    print(y)

    yyy <- lapply(y, function(u)
    {
        yy <- system(u, intern = TRUE)
        yy
    })

    yyyy <- unlist(yyy)

    return(yyyy)
}

testAlignment <- function(output.dir, gene.model.file, genome.index)
{
    cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd4 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W"

    re <- parserreadfiles(output.dir, "fastq")
    res <- re$input

    xx <- lapply(res, function(u)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        if (regexpr(pattern = "_", file_name) != -1)
        {
            # cat('match:',file_name,'\n')

            p <- regexpr(pattern = "_", file_name)
            pp <- p - 1
            x <- substr(file_name, 1, pp)

        } else
        {
            # cat('no match:',file_name,'\n')
            x <- file_name
        }

        x
    })


    xxx <- unique(unlist(xx))
    res2 <- unlist(res)

    cmd5 = "tophat --library-type fr-unstranded -g 1 -G"
    cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
    cmd7 = "tophat --library-type fr-secondstrand -g 1 -G"
    cmd8 = "-p 4 -o"
    cmd9 = "mv"

    for (i in 1:length(xxx))
    {
        sample.name <- xxx[i]
        sample.name.out.dir <- file.path(output.dir, sample.name)

        if (!dir.exists(sample.name.out.dir))
        {
            dir.create(sample.name.out.dir, recursive = TRUE)
        }

        y <- res2[grep(xxx[i], res2)]

        # print(y) print(length(y))

        if (length(y) == 2)
        {
            yy1 <- basename(y[1])
            yy2 <- basename(y[2])

            p1 <- regexpr(pattern = "_", yy1)
            pp1 <- p1 + 1
            x1 <- substr(yy1, pp1, pp1)

            p2 <- regexpr(pattern = "_", yy2)
            pp2 <- p2 + 1
            x2 <- substr(yy2, pp2, pp2)


            sample.name.out.dir.1 = file.path(sample.name.out.dir, paste0("Us",
                x1, x2))
            sample.name.out.dir.2 = file.path(sample.name.out.dir, paste0("Us",
                x2, x1))
            sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                x1, x2))
            sample.name.out.dir.4 = file.path(sample.name.out.dir, paste0("Fs",
                x2, x1))
            sample.name.out.dir.5 = file.path(sample.name.out.dir, paste0("Ss",
                x1, x2))
            sample.name.out.dir.6 = file.path(sample.name.out.dir, paste0("Ss",
                x2, x1))

            if (!dir.exists(sample.name.out.dir.1))
            {
                dir.create(sample.name.out.dir.1, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.2))
            {
                dir.create(sample.name.out.dir.2, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.3))
            {
                dir.create(sample.name.out.dir.3, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.4))
            {
                dir.create(sample.name.out.dir.4, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.5))
            {
                dir.create(sample.name.out.dir.5, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.6))
            {
                dir.create(sample.name.out.dir.6, recursive = TRUE)
            }

            cmd10 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.1,
                genome.index, y[1], y[2], sep = " ")
            cmd11 = paste(cmd4, cmd2, cmd10)
            print(cmd10)


            system(cmd11)

            cmd12 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.2,
                genome.index, y[2], y[1], sep = " ")
            cmd13 = paste(cmd4, cmd2, cmd12)
            system(cmd13)

            cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                genome.index, y[1], y[2], sep = " ")
            cmd15 = paste(cmd4, cmd2, cmd14)
            system(cmd15)

            cmd16 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.4,
                genome.index, y[2], y[1], sep = " ")
            cmd17 = paste(cmd4, cmd2, cmd16)
            system(cmd17)

            cmd18 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.5,
                genome.index, y[1], y[2], sep = " ")
            cmd19 = paste(cmd4, cmd2, cmd18)
            system(cmd19)

            cmd20 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.6,
                genome.index, y[2], y[1], sep = " ")
            cmd21 = paste(cmd4, cmd2, cmd20)
            system(cmd21)
        } else
        {
            sample.name.out.dir.7 = file.path(sample.name.out.dir, "Us")
            sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")
            sample.name.out.dir.9 = file.path(sample.name.out.dir, "Ss")

            if (!dir.exists(sample.name.out.dir.7))
            {
                dir.create(sample.name.out.dir.7, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.8))
            {
                dir.create(sample.name.out.dir.8, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.9))
            {
                dir.create(sample.name.out.dir.9, recursive = TRUE)
            }

            cmd22 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.7,
                genome.index, y[1], sep = " ")
            cmd23 = paste(cmd4, cmd2, cmd22)
            # print(cmd23)

            system(cmd23)

            cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                genome.index, y[1], sep = " ")
            cmd25 = paste(cmd4, cmd2, cmd24)
            system(cmd25)

            cmd26 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.9,
                genome.index, y[1], sep = " ")
            cmd27 = paste(cmd4, cmd2, cmd26)
            system(cmd27)
        }

    }
}

useTophat4Alignment <- function(input.fastq.files.dir, output.dir, gene.model.file = NULL,
    genome.index, cmd.input)
    {
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    if (cmd.input == "General")
    {
        cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
        cmd4 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W"
    } else
    {

    }

    re <- parserreadfiles(input.fastq.files.dir, "fastq")
    res <- re$input

    xx <- lapply(res, function(u)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        if (regexpr(pattern = "_", file_name) != -1)
        {
            # cat('match:',file_name,'\n')

            p <- regexpr(pattern = "_", file_name)
            pp <- p - 1
            x <- substr(file_name, 1, pp)

        } else
        {
            # cat('no match:',file_name,'\n')
            x <- file_name
        }

        x
    })

    xxx <- unique(unlist(xx))
    res2 <- unlist(res)

    cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
    cmd8 = "-p 4 -o"
    cmd9 = "mv"

    for (i in 1:length(xxx))
    {
        sample.name <- xxx[i]
        sample.name.out.dir <- file.path(output.dir, sample.name)

        if (!dir.exists(sample.name.out.dir))
        {
            dir.create(sample.name.out.dir, recursive = TRUE)
        }

        y <- res2[grep(xxx[i], res2)]

        # print(y) print(length(y))

        if (length(y) == 2)
        {
            yy1 <- basename(y[1])
            yy2 <- basename(y[2])

            p1 <- regexpr(pattern = "_", yy1)
            pp1 <- p1 + 1
            x1 <- substr(yy1, pp1, pp1)

            p2 <- regexpr(pattern = "_", yy2)
            pp2 <- p2 + 1
            x2 <- substr(yy2, pp2, pp2)


            sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                x1, x2))


            if (!dir.exists(sample.name.out.dir.3))
            {
                dir.create(sample.name.out.dir.3, recursive = TRUE)
            }

            cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                genome.index, y[1], y[2], sep = " ")
            # cmd15= paste(cmd.input,cmd14)

            cmd15 = paste(cmd4, cmd2, cmd14)

            system(cmd15)

        } else
        {
            sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")

            if (!dir.exists(sample.name.out.dir.8))
            {
                dir.create(sample.name.out.dir.8, recursive = TRUE)
            }
            cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                genome.index, y[1], sep = " ")
            # cmd25= paste(cmd.input,cmd24)
            cmd25 = paste(cmd4, cmd2, cmd24)
            system(cmd25)
        }

    }

}

#' @examples
#' R -e 'library(ThreeUTR);ThreeUTR:::subsetBam('/scratch/projects/bbc/aiminy_project/DoGs/Bam_split/','$HOME/1833_common_gene.bed','/scratch/projects/bbc/aiminy_project/DoGs/BAMSubSet',BigMem=FALSE)'
#'
subsetBam <- function(input.bam.file.dir, region.bed.file, output.bw.file.dir,
    BigMem = FALSE)
    {
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, res, region.bed, BigMem,
        output.bw.file.dir)
        {
        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            if (BigMem == TRUE)
            {
                cmd0 = "72:00 -n 16 -q bigmem -R 'rusage[mem=36864] span[ptile=8]' -u aimin.yan@med.miami.edu"
            } else
            {
                cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            }

            job.name = paste0("subsetBam.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))
            cmd2 = paste("samtools view -b -L", region.bed.file, res[[u]], "\\>",
                file.path(output.bw.file.dir, paste0(file_name, "_subset.bam")),
                sep = " ")

            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            cmd3 = paste("samtools view -b -L", region.bed.file, res[[u]], ">",
                file.path(output.bw.file.dir, paste0(file_name, "_subset.bam")),
                sep = " ")
        }

        cmd <- cmd3
        system(cmd)
        cmd
    }, m.id, res, region.bed, BigMem, output.bw.file.dir)

}
