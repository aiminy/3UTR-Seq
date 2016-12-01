#!/bin/bash

#Usage: sh /home/aiminyan/Code/BashRunCovertPeak2Bed.sh input_files_dir output_dir_name
#Example: sh /home/aiminyan/Code/BashRunCovertPeak2Bed.sh /media/H_driver/2016/Yang/PeakCall/ /media/H_driver/2016/Yang/BedFromPeakCall/

DIR="$1"

results_dir="$2"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi

#for file in $(ls $DIR*bam_mm_summits.bed)
for file in $(ls $DIR*_bam_mm_1.00e-07_peaks.narrowPeak)
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

done
