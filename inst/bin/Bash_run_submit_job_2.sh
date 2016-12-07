#!/bin/sh

TOPHAT_BINARY=/path/to/tophat/executable/tophat
GENE_REFERENCE=/ProjectName/reference/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf
BOWTIE_INDEX=/ProjectName/reference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
P=4 

for SAMPLE_ID in {1..N}
do
cat > /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh <<EOF
$TOPHAT_BINARY -G $GENE_REFERENCE -p $P -o /ProjectName/alignments/sample${SAMPLE_ID} $BOWTIE_INDEX /ProjectName/data/sample${SAMPLE_ID}_1.fastq /ProjectName/data/sample${SAMPLE_ID}_2.fastq
mv /ProjectName/alignments/sample${SAMPLE_ID}/accepted_hits.bam /ProjectName/alignments/sample${SAMPLE_ID}.bam
EOF

qsub -l mf=20G,h_vmem=5G -pe local $P -m e -M myemail@email.com /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh

done