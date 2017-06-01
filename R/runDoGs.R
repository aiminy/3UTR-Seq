#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runDoGs("SRP058633",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline")'
#'
#'
#' sample.info.file=file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt")
#' gene.gtf="/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
#' genome.index="/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
#' processed.gene.gtf="/projects/ctsi/bbc/aimin/annotation/"
#' output.dir="/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline"
#'
#'
runDoGs <- function(sra.accession.number,sample.info.file,gene.gtf,genome.index,processed.gene.gtf,output.dir) {

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

# re <- ThreeUTR:::useWget2Download(sra.accession.number,file.path(output.dir,"SRAFiles"))
#
# # This setting works
# Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useFastqDumpConvertSra2Fastq('
# input=file.path(output.dir,"SRAFiles")
# output=file.path(output.dir,"Fastqfiles")
# Rfun2 <- ',wait.job.name = "wgetDownload")'
# Rfun3 <- paste0(Rfun1,'\\"',input,'\\"',',\\"',output,'\\"',Rfun2)
#
# run1 <- createBubRfun(Rfun3,"sra2fastq","wgetDownload")
# system(run1)

Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useTophat4Alignment2('
input=file.path(output.dir,"Fastqfiles")
output=file.path(output.dir,"Alignment")
gene.gtf=gene.gtf
genome.index=genome.index
wait.job.name = 'wait.job.name = NULL'
Rfun2 <- ')'

Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\",','\\"',gene.gtf,'\\",','\\"',genome.index,'\\",',wait.job.name)
Rfun <-paste0(Rfun1,Rinput,Rfun2)

test <- createBubRfun(Rfun,"Alignment",NULL)
system(test)

# useTophat4Alignment2(file.path(output.dir,"FastqFiles"),file.path(output.dir,"Alignment"),gene.gtf,genome.index,"parallel",wait.job.name="sra2fastq")

# processBamFiles(file.path(output.dir,"Alignment"),file.path(output.dir,"ProcessedBam"),wait.job.name="sra2fastq")
#
# convertbam2bed(file.path(output.dir,"ProcessedBam"),file.path(output.dir,"Bed"),wait.job.name="sra2fastq")
#
# removeReadsOnExonIntron(file.path(output.dir,"Bed"),processed.gene.gtf,file.path(output.dir,"BedRmEandI"),wait.job.name="sra2fastq")
#
# getCount4Downstream(file.path(output.dir,"BedRmEandI"),processed.gene.gtf,file.path(output.dir,"Counts"),wait.job.name="sra2fastq")
#
# res <- convertCountFile2Table(file.path(output.dir,"Counts"),"*.txt",wait.job.name="sra2fastq")
#
# res.new <- ThreeUTR:::matchAndDE(res,sample.info.file,group.comparision = c("condition","Untreated","Treated"),wait.job.name="sra2fastq")

}

#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runDoGsOnCluster("SRP058633",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline")'

runDoGsOnCluster <- function(sra.accession.number,sample.info.file,gene.gtf,genome.index,processed.gene.gtf,output.dir) {

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  re <- ThreeUTR:::useWget2Download(sra.accession.number,file.path(output.dir,"SRAFiles"))

  # This setting works
  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertSra2FastqUseJobArray('
  input=file.path(output.dir,"SRAFiles")
  output=file.path(output.dir,"Fastqfiles")
  Rfun2 <- ')'
  Rfun3 <- paste0(Rfun1,'\\"',input,'\\"',',\\"',output,'\\"',Rfun2)

  sra2fastq <- createBsubJobArrayRfun(Rfun3,"sra2fastq[1-8]","wgetDownload")
  system(sra2fastq)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::alignmentUseJobArray('
  input=file.path(output.dir,"Fastqfiles")
  output=file.path(output.dir,"Alignment")
  gene.gtf=gene.gtf
  genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\",','\\"',gene.gtf,'\\",','\\"',genome.index,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  alignment <- createBsubJobArrayRfun(Rfun,"Alignment[1-8]","sra2fastq")
  system(alignment)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::processBamFilesUseJobArray('
  input=file.path(output.dir,"Alignment")
  output=file.path(output.dir,"Bam")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  processbam <- createBsubJobArrayRfun(Rfun,"ProcessBam[1-8]","Alignment")
  system(processbam)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertBam2bedUsingJobArray('
  input=file.path(output.dir,"Bam")
  output=file.path(output.dir,"BedFromBam")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  bam2bed <- createBsubJobArrayRfun(Rfun,"Bam2Bed[1-8]","ProcessBam")
  system(bam2bed)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::removeReadsOnExonIntronUsingJobArray('
  input=file.path(output.dir,"BedFromBam")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"BedRmExonIntron")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  rm.exon.intron <- createBsubJobArrayRfun(Rfun,"RmExonIntron[1-8]","Bam2Bed")

  #rm.exon.intron <- createBsubJobArrayRfun(Rfun,"RmExonIntron[1-8]",NULL)
  system(rm.exon.intron)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::getCount4DownstreamUsingJobArray('
  input=file.path(output.dir,"BedRmExonIntron")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"Counts")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  counting <- createBsubJobArrayRfun(Rfun,"Count[1-8]","RmExonIntron")
  #counting <- createBsubJobArrayRfun(Rfun,"Count[1-8]",NULL)
  system(counting)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);library(org.Hs.eg.db);re <- ThreeUTR:::CountAndDE('
  input=file.path(output.dir,"Counts")
  #processed.gene.gtf=processed.gene.gtf
  sample.info.file=sample.info.file
  output=file.path(output.dir,"Results")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',sample.info.file,'\\",','\\"',output,'\\"')
  Rfun <- paste0(Rfun1,Rinput,Rfun2)

  get.DE <- createBsubJobArrayRfun(Rfun,"Summary[1]","Count")
 #get.DE <- createBsubJobArrayRfun(Rfun,"Summary[1]",NULL)
  system(get.DE)

}

#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runSpliceJunction("/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline2")'

runSpliceJunction <- function(output.dir,wait.job=NULL) {

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::processSpliceJunctionFilesUseJobArray('
  input=file.path(output.dir,"Alignment")
  output=file.path(output.dir,"SpliceBed")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  processbam <- createBsubJobArrayRfun(Rfun,"ProcessSpliceJunction[1-8]",wait.job.name=wait.job)

  system(processbam)

}

#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runRmExonAndIntronFromSplicingJunctions("/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline2","/projects/ctsi/bbc/aimin/annotation/")'
#'
runRmExonAndIntronFromSplicingJunctions <- function(output.dir, processed.gene.gtf,wait.job=NULL) {
  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::removeReadsOnExonIntronUsingJobArray('
  input=file.path(output.dir,"SpliceBed")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"SpliceBedRmExonIntron")
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  rm.exon.intron <- createBsubJobArrayRfun(Rfun,"RmExonIntron[1-8]",wait.job.name=wait.job)

  system(rm.exon.intron)
}

#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runGenerateSubSetBam("/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline2","/projects/ctsi/bbc/aimin/annotation/")'
#'
runGenerateSubSetBam <- function(output.dir, processed.gene.gtf,wait.job=NULL) {
  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::generateSubSetBam('
  input=file.path(output.dir,"Bam")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"BamExample")
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  rm.exon.intron <- createBsubJobArrayRfun(Rfun,"GeneRateBamExample[1-8]",wait.job.name=wait.job)

  system(rm.exon.intron)
}

#' If you already have BAM files available, you can use this function to perform analysis
#'
#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runDoGsOnClusterStartFromBam(file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/Example")'

runDoGsOnClusterStartFromBam <- function(sample.info.file,gene.gtf,genome.index,processed.gene.gtf,output.dir) {

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertBam2bedUsingJobArray('
  input=file.path(output.dir,"Bam")
  output=file.path(output.dir,"BedFromBam")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  #bam2bed <- createBsubJobArrayRfun(Rfun,"Bam2Bed[1-8]","ProcessBam")

  bam2bed <- createBsubJobArrayRfun(Rfun,"Bam2Bed[1-8]",NULL)

  system(bam2bed)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::removeReadsOnExonIntronUsingJobArray('
  input=file.path(output.dir,"BedFromBam")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"BedRmExonIntron")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  rm.exon.intron <- createBsubJobArrayRfun(Rfun,"RmExonIntron[1-8]","Bam2Bed")

  #rm.exon.intron <- createBsubJobArrayRfun(Rfun,"RmExonIntron[1-8]",NULL)
  system(rm.exon.intron)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::getCount4DownstreamUsingJobArray('
  input=file.path(output.dir,"BedRmExonIntron")
  processed.gene.gtf=processed.gene.gtf
  output=file.path(output.dir,"Counts")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',processed.gene.gtf,'\\",','\\"',output,'\\"')
  Rfun <-paste0(Rfun1,Rinput,Rfun2)

  counting <- createBsubJobArrayRfun(Rfun,"Count[1-8]","RmExonIntron")
  #counting <- createBsubJobArrayRfun(Rfun,"Count[1-8]",NULL)
  system(counting)

  Rfun1 <- 'library(ChipSeq);library(ThreeUTR);library(org.Hs.eg.db);re <- ThreeUTR:::CountAndDE('
  input=file.path(output.dir,"Counts")
  #processed.gene.gtf=processed.gene.gtf
  sample.info.file=sample.info.file
  output=file.path(output.dir,"Results")
  #gene.gtf=gene.gtf
  #genome.index=genome.index
  #wait.job.name = 'wait.job.name = "sra2fastq"'
  Rfun2 <- ')'

  Rinput <- paste0('\\"',input,'\\",','\\"',sample.info.file,'\\",','\\"',output,'\\"')
  Rfun <- paste0(Rfun1,Rinput,Rfun2)

  get.DE <- createBsubJobArrayRfun(Rfun,"Summary[1]","Count")
  #get.DE <- createBsubJobArrayRfun(Rfun,"Summary[1]",NULL)
  system(get.DE)

}
