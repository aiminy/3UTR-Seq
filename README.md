**An pipeline for processing and analyzing DOGs data**

**Abstract**

Vilborg et al. discovered a new class of long chromatin-associated RNAs, named as a ‘downstream of gene’-containing transcript (DoG). They found that that DoGs is upregulated in response to osmotic stress. Their bioinformatic analysis identified KCl-induced DoGs downstream of more than 10% of all human protein-coding genes. However there is a lack of streamlined bioinforamtics pipleline for processing and analyzing this type of RNA sequence data. This paper aims to develep a pipeline to fill in these gaps for further research on DoGs. 

**Method**

+ To Install

In R console

library(devtools)

install_github("aiminy/3UTR-Seq",ref = '3UTR')

Install update ThreeUTR without restarting R

detach("package:ThreeUTR", unload=TRUE)

library(ThreeUTR)

In pegasus terminal 

R -e 'library(devtools);install_github("aiminy/3UTR-Seq",ref = "3UTR")'

To start analysis
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
You choose 7, you will see the following:
```
# You choose to process bam files
# Firstly, we remove reads for which the insert size is < 100 or > 200 by using the follwoing procedure
R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::splitBam('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs/Bam_split')'

#Secondly, we convert bam files to strand-specific BigWig files by using the following procedure
R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2StrandBw2("/scratch/projects/bbc/aiminy_project/DoGs/Bam_split/","/scratch/projects/bbc/aiminy_project/DoGs/BW_split")'
```
After this, you can load strand-specific BigWig to IGV to visualize DoGs. The following is an example to show DoGs
![Image of DoGs](inst/extdata/DoGs.png)
The above Figure shows two samples under two conditions(KCl-treated and untreated). There are 4 panels, the top two panels show the coverage of forward and reverse strand for a sample treated by KCl;the botton two panels show the coverage of forward and reverse strand for an untreated sample. It is clear that we can observe a DoG in the downstream of CXXC4 gene in KCl-treated sample(labeled by red arrow).

Since we are interested in intergenic reads instead of reads overlapping with exons and introns, so we perform the following procedure:
```{r}
# Step1: Convert the aligned bam files to bed files

R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertbam2bed('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs')'

# Step2: Remove reads overlappping with exons and intron firstly

R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::removeReadsOnExonIntron("/scratch/projects/bbc/aiminy_project/DoGs/BedFileFromBam","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/BedRmExonIntron")'

# Step3: get counts of intergenic reads with 45kb downstream of transcripts 

R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::getCount4Downstream(""/scratch/projects/bbc/aiminy_project/DoGs/BedRmExonIntron","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/Counts45KB")'

# Step4: convert count files to count table, and perform differential DoGs analysis

res21 <- ThreeUTR:::sumCount4Downstream("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),c(2,1))
```
The following Figure shows differentila DoGs analysis results for two transcripts for CXXC4(shown on the above Figure) under two conditions(KCl-treated and untreated). It is clear that there are differential DoGs in the downstream of CXXC4 gene between two conditions, which is consistent with the above visualization
![Image of DeDoGs](inst/extdata/De.png)

**Summary**

Vilborg identifyied a new class of long chromatin-associated RNA,named as DoGs. RNA sequencing data for discoverying DoGs is a new type of data, and the streamlined bioinformatics pipeline for processing and analyzing this new type of RNA sequencing data is not available yet. Here we developed a streamlined pipeline to proceesing DoGs RNA-Seq data. We believe that this bioinforamtics can enhance the research related to DoGs.
