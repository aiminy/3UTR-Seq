# 3UTR-Seq
An R package for processing and analyzing 3UTR,DOGs and 5UTR data

#### To Install

##### In R console

library(devtools)

install_github("aiminy/3UTR-Seq",ref = '3UTR')

##### In pegasus terminal 

R -e 'library(devtools);install_github("aiminy/3UTR-Seq",ref = "3UTR")'

##### To start analysis
```{r}
Rscript ~/R/lib64/R/library/ThreeUTR/bin/rscript/utr_interactive.r
```
You will see a list of job options, then you can choose the job you want to do
```{}
Your operating system is:  linux 
Choose analysis: 
	Avaliable analysis: 

    1. Download_SRA 

    2. SRA2Fastq 

    3. subset_Fastq_file 

    4. testAlignment 

    5. Alignment 

    6. QC 

    7. Bam2Bw 

    8. Counts 

    9. DE_analysis 

    10. All
```
You choose 1, you will see the following:
```{}
You choose to download SRA files, please define the following setting parameters: 
 sra.accession.number (ex: SRP058633 Use an example from DoGs paper)
 output.dir (ex:DoGsExample)
```
You choose 2, you will see the following:
```{}
You choose to convert SRA files to Fastq files, please define the following setting parameters: 
 sra.file.dir (ex: /nethome/axy148/DoGsExample)
 output.dir (ex:DoGsFastq)
```
You choose 5, you will see the fllowing:
```
You choose to perform alignment, please define the following setting parameters: 
 fastq.file.dir (ex:/scratch/projects/bbc/aiminy_project/DoGsFastq)
 output.dir (ex:/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBam)
 gene.model.file(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf or /projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/refGene.txt
 genome.index(ex:/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome)
```
