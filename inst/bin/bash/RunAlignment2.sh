#! /bin/bash

#awk '{print "bowtie -p 8 --chunkmbs 2000 --sam mm10_index/genome -1 <(gunzip " $1")","-2 <(gunzip "$2")"}' Fq_24.txt

#sh RunAlignment2.sh Fq_24.txt  
while read line; do

f1=`echo "$line" | awk -F"\t" '{print $1}'`
f2=`echo "$line" | awk -F"\t" '{print $2}'`

#f=`echo "$line"`

echo "$f1\t$f2"
sample_name=`echo "$f1" | awk -F"_" '{print $3}'`

#echo "$sample_name"_aligned_reads.sam

#echo "bowtie -p 8 --chunkmbs 2000 --sam mm10_index/genome -1 <(gunzip "$f1") -2 <(gunzip "$f2") "$sample_name"_aligned_reads.sam&" >> Run_Fq_24.sh
#echo "tophat mm10_index_bt2/genome ~/RNAseqData/"$f1" ~/RNAseqData/"$f2" -o "$sample_name"_tophat_out" > Run_"$sample_name"_tophat.sh 


cat > Run_"$sample_name"_tophat.sh <<EOF
tophat -G genes.gtf -p 4 -o "$sample_name"_tophat_out mm10_index_bt2/genome ~/RNAseqData/"$f1" ~/RNAseqData/"$f2"
mv "$sample_name"_tophat_out/accepted_hits.bam "$sample_name".bam
EOF

qsub -l mf=20G,h_vmem=5G -pe local 4 -m e -M jlfmssm@gmail.com Run_"$sample_name"_tophat.sh

done < $1


#for SAMPLE_ID in {1..N}
#do
#cat > /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh <<EOF
#$TOPHAT_BINARY -G $GENE_REFERENCE -p $P -o /ProjectName/alignments/sample${SAMPLE_ID} $BOWTIE_INDEX /ProjectName/data/sample${SAMPLE_ID}_1.fastq /ProjectName/data/sample${SAMPLE_ID}_2.fastq
#mv /ProjectName/alignments/sample${SAMPLE_ID}/accepted_hits.bam /ProjectName/alignments/sample${SAMPLE_ID}.bam
#EOF

#qsub -l mf=20G,h_vmem=5G -pe local $P -m e -M myemail@email.com /ProjectName/alignments/scripts/tophat_${SAMPLE_ID}.sh

#done