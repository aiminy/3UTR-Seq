#Usage:sh /home/aiminyan/Code/Run_JunctionSeq_2.sh /media/H_driver/2015/Nimer_Cheng/Bam_files_6.txt

#Use starnd inforamtion and specify to the second strand

while read line; do

f=`echo "$line"`

echo "$f"

dir_name=$(dirname "$f")
file_name=$(basename "$f")

echo "$dir_name"
echo "$file_name"

sample_name=$(echo "$file_name" | awk -F"." '{print $1}')
echo "$sample_name"

#output_name=$("$dir_name"/"$sample_name"_second_strand_based)
#sample_name=`echo "$f" | awk -F"." '{print $1}'`
#echo "$output_name"

java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded --stranded --stranded_fr_secondstrand --noGzipOutput --maxReadLength 1000 --keepMultiMapped "$f" /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf "$dir_name"/"$sample_name"_second_strand_based

done < $1
