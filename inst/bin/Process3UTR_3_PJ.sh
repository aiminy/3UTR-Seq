#!/bin/bash

#Usage: sh /home/aiminyan/Code/Process3UTR_3.sh input_files_dir
#Example: sh /home/aiminyan/Code/Process3UTR_3_PJ.sh /media/H_driver/PJ/Alignment/75_S1_STAR/

DIR="$1*"

results_dir="$1Results/"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi


for file in $(ls $DIR.bam)

do

f=`echo "$file"`

echo "$f"

#bam_file="$f"

dir_name=$(dirname "$f")

file_name=$(basename "$f")

echo "$dir_name"
echo "$file_name"

sample_name=$(echo "$file_name" | awk -F"." '{print $1}')

echo "$sample_name"

cat > ~/Code/Run_"$sample_name"_3UTR.sh <<EOF

#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"$sample_name".alignments.bam|cut -f 9 > "$results_dir""$sample_name".InsertSizeMetrics.txt
#samtools stats "$DIR""$sample_name".alignments.bam > "$results_dir""$sample_name".summary.txt

/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in "$DIR""$sample_name".sorted.bam --min 100 --max 200 --out "$results_dir""$sample_name".sorted.split.bam  

EOF

sh ~/Code/Run_"$sample_name"_3UTR.sh& 

done
