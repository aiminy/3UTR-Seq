#!/bin/bash   
#sh ~/Script_bash/BashRunConvertBamtoBigwig.sh FileBam.txt                                                                                                                                                 

#mkdir /scratch/projects/bbc/Peak_chip_seq

while read line; do

#f1=`echo "$line" | awk -F"\t" '{print $1}'`
#f2=`echo "$line" | awk -F"\t" '{print $2}'`
f=`echo "$line"`

#echo "$f"

dir_name=$(dirname "$f")
file_name=$(basename "$f")

#echo "$dir_name"

#echo "$file_name"

sample_name=`echo "$file_name" | awk -F"." '{print $1}'`

#echo "$sample_name"

echo "$f"
#echo "$dir_name"

#echo axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam
#echo axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam.bai
#echo axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bw

#scp -i ~/.ssh/id_rsa axy148@pegasus.ccs.miami.edu:"$f" /media/DATA/metaseq/metaseq-Lluis-data/

#sshpass -p "$password" scp axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam /media/DATA/metaseq/metaseq-Lluis-data/
#sshpass -p "$password" scp axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam.bai /media/DATA/metaseq/metaseq-Lluis-data/
#sshpass -p "$password" scp axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bw /media/DATA/metaseq/metaseq-Lluis-data/

scp -i ~/.ssh/id_rsa axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam /media/DATA/metaseq/metaseq-Lluis-data/
scp -i ~/.ssh/id_rsa axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bam.bai /media/DATA/metaseq/metaseq-Lluis-data/
scp -i ~/.ssh/id_rsa axy148@pegasus.ccs.miami.edu:"$dir_name"/"$sample_name".bw /media/DATA/metaseq/metaseq-Lluis-data/


done < $1
