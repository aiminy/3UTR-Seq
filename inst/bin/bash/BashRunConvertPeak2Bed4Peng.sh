#!/bin/bash

#Usage: sh /home/aiminyan/Code/BashRunCovertPeak2Bed4Peng.sh input_files_dir output_dir_name
#Example: sh /home/aiminyan/Code/BashRunCovertPeak2Bed4Peng.sh /media/H_driver/2016/Yang/MACS/MACS/ /media/H_driver/2016/Yang/BedFromPeng/

DIR="$1"

results_dir="$2"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi

#for file in $(ls $DIR*bam_mm_summits.bed)
#for file in $(ls $DIR*_bam_mm_1.00e-07_peaks.narrowPeak)
for file in $(ls $DIR*.bed)
do

f=`echo "$file"`

echo "$f"

dir_name=$(dirname "$f")

file_name=$(basename "$f")

echo "$dir_name"
echo "$file_name"

sample_name=$(echo "$file_name" | awk -F"." '{print $1}')

echo "$sample_name"

cut -f1,2,3,4 "$f" > "$results_dir""$sample_name"_from_PeakCall4Yang.bed

awk -F"\t" '{print $1,"_",$2,"_",$3}' OFS="" "$f" > "$results_dir""$sample_name"_from_PeakCall4Yang_combine_chr_st_end.bed

awk -F"\t" '{print $1,"_",$2}' OFS="" "$f" > "$results_dir""$sample_name"_from_PeakCall4Yang_combine_chr_st.bed

done
