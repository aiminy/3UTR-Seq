#!/bin/bash
#declare -A chro

sort -k1,1 <(tail -n +69 mgp.v3.indels.rsIDdbSNPv137.addchr.vcf) | sort --stable -n -k2,2 > mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.vcf

chro=('chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chrM' 'chrX' 'chrY')

for each in ${chro[@]}; do 
echo $each;
cat <(grep $each mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.vcf) >> mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.2.vcf
done

cat <(head -n68 mgp.v3.indels.rsIDdbSNPv137.addchr.vcf) mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.2.vcf > mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.with.header.vcf

/media/DATA/NGS_tools/gatk_pipeline/resources/usr/bin/vcfsorter.pl mm10.dict mgp.v3.indels.rsIDdbSNPv137.addchr.vcf > mgp.v3.indels.rsIDdbSNPv137.by.sorter.vcf 2>STDERR

#Generate a vcf having indel only
vcftools --vcf mgp.v3.indels.rsIDdbSNPv137.by.sorter.vcf --keep-only-indels --out mgp.v3.indels.rsIDdbSNPv137.by.sorter.indel.only.vcf --recode --recode-INFO-all

java -Xmx1g -jar /path/to/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator \
  -R /path/to/reference.fasta \
  -o /path/to/output.intervals \
  -known /path/to/indel_calls.vcf

#IndelRealigner to be tested
java -Xmx1g -jar $GATK -T RealignerTargetCreator -R /home/aiminyan/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -o ./output.intervals -known ./
java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir jar $GATK -I <lane-level.bam> -R <ref.fasta> -T IndelRealigner -targetIntervals output.intervals -o <realignedBam.bam> -known --consensusDeterminationModel KNOWNS_ONLY -LOD 0.4

#BaseRecalibrator OK
java -Xmx4g -jar $GATK -T BaseRecalibrator -R /home/aiminyan/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -knownSites /media/DATA/Box/zLi/zLi/myWork/mgp.v3.indels.rsIDdbSNPv137.by.sorter.vcf -I /media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_.bam -o /media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_Recal.grp

#PrintReads OK
java -jar $GATK -T PrintReads -R /home/aiminyan/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -I /media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_.bam -BQSR /media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_Recal.grp -o /media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_recali.bam

#Downsampling ok
java -jar ~/workspace/rnaseq/tools/picard-tools-1.140/picard.jar DownsampleSam I=/media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_.bam O=/media/DATA/Box/zLi/zLi/myWork/zLi_MouseExome_201519492-01_S_3_008.bam R=1234 P=0.08

#INPUT=File
#I=File                        The input SAM or BAM file to downsample.  Required. 

#OUTPUT=File
#O=File                        The output, downsampled, SAM or BAM file to write.  Required. 

#RANDOM_SEED=Long
#R=Long                        Random seed to use if reproducibilty is desired.  Setting to null will cause multiple 
#                              invocations to produce different results.  Default value: 1. This option can be set to 
#                              'null' to clear the default value. 

#PROBABILITY=Double
#P=Double         

#cat <(grep chr10 mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.vcf) <(grep chr11 mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.vcf) <> mgp.v3.indels.rsIDdbSNPv137.addchr.sorted.2.vcf
#chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr1 chr2 chr4 chr5 chr6 chr7 chr#8 chr9 chrM chrX chrY