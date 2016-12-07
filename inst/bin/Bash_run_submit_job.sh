#! /bin/bash

#sh Bash_run_submit_job.sh  
while read line; do

#f1=`echo "$line" | awk -F"\t" '{print $1}'`
#f2=`echo "$line" | awk -F"\t" '{print $2}'`
f=`echo "$line"`

echo "$f"
#echo "$f1\t$f2"
sample_name=`echo "$f" | awk -F"_" '{print $2}'`

#echo "$sample_name"_aligned_reads.sam

#echo "bowtie -p 8 --chunkmbs 2000 --sam mm10_index/genome -1 <(gunzip "$f1") -2 <(gunzip "$f2") "$sample_name"_aligned_reads.sam&" >> Run_Fq_24.sh

#echo "tophat mm10_index_bt2/genome ~/RNAseqData/"$f1" ~/RNAseqData/"$f2" -o "$sample_name"_tophat_out" > Run_"$sample_name"_tophat.sh 

bsub -P "$sample_name"_run < "$f"

done < $1
