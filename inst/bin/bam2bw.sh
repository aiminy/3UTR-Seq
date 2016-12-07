#! /bin/bash

for bam in *bam
do 
echo $bam 
genomeCoverageBed -ibam $bam -bg -g hg19.genome.info > $(basename $bam .bam).bdg
done

for bdg in *bdg
do 
echo $bdg
bdg2bw $bdg hg19.genome.info 
done