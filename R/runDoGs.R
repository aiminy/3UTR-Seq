#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::runDoGs("SRP058633",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/TestPipeline")'
#'
runDoGs <- function(sra.accession.number,sample.info.file,gene.gtf,genome.index,processed.gene.gtf,output.dir) {

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

re <- ThreeUTR:::useWget2Download(sra.accession.number,file.path(output.dir,"SRAFiles"))

Rfun1 <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useFastqDumpConvertSra2Fastq('
input=file.path(output.dir,"SRAFiles")
output=file.path(output.dir,"Fastqfiles")
Rfun2 <- ',wait.job.name = "wgetDownload")'
Rfun <- paste0(Rfun1,'"',input,'"',',"',output,'"',Rfun2)

test <- createBubRfun(Rfun,"sra2fastq","wgetDownload")
system(Rfun)

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
