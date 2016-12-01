#Usage:sh /home/aiminyan/Code/Run_use_gene_counts_in_Junction_Seq.sh /media/H_driver/2015/Nimer_Cheng/6.samples.QC.geneCounts.formatted.for.DESeq.txt

#Use gene counts

while read line; do

f=`echo "$line"`

#echo "$f"

dir_name=$(dirname "$f")
file_name=$(basename "$f")

echo "$dir_name"
#echo "$file_name"

sample_name=$(echo "$dir_name" | awk -F"/" '{print $6}')


#echo "$f ""$sample_name"

num_line=`wc -l "$f"`

num_line_2=$(echo "$num_line" | awk -F" " '{print $1}')

#echo "$num_line_2"
num_line_3=$((num_line_2-5))

#echo "$num_line_3"

sed -n 1,"$num_line_3"p "$f" > "$dir_name"_gene_counts_only.txt  

dir_name_2=$(dirname "$dir_name"_gene_counts_only.txt)

echo "$dir_name_2"

echo "$dir_name"_gene_counts_only.txt "$sample_name" >> "$dir_name_2""/"combine_6.txt

#output_name=$("$dir_name"/"$sample_name"_second_strand_based)
#sample_name=`echo "$f" | awk -F"." '{print $1}'`
#echo "$output_name"

#java -Xmx4000M -jar /home/aiminyan/QoRTs/QoRTsFullExampleData/QoRTsRelease/QoRTs.jar QC --singleEnded --stranded --stranded_fr_secondstrand --noGzipOutput --maxReadLength 1000 --keepMultiMapped "$f" /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf "$dir_name"/"$sample_name"_second_strand_based

done < $1
