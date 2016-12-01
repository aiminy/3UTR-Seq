#! /bin/bash

#sh Bash_run_submit_job_3.sh FileBam.txt  
while read line; do

f=`echo "$line"`

echo "$f"

sample_name=`echo "$f" | awk -F"." '{print $1}'`

echo "$sample_name"

cat > Run_"$sample_name"_htseq_count.sh <<EOF



htseq-count -f bam "$f" genes.gtf > "$sample_name"_raw_count_0.txt
EOF

bsub -P Bioinformatics4count < Run_"$sample_name"_htseq_count.sh

done < $1