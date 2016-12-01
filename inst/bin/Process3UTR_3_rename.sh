#!/bin/bash

#Usage: sh /home/aiminyan/Code/Process3UTR_3.sh input_files_dir
#Example: sh /home/aiminyan/Code/Process3UTR_3.sh /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/

DIR="$1"

results_dir="$1Results4NewData/"

if [ -e $results_dir ]; then echo "already exists";
else
    mkdir $results_dir

fi


for file in $(ls $DIR*.bam)

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

#cat > ~/Code/Run_"$sample_name"_3UTR.sh <<EOF

#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"$sample_name".alignments.bam|cut -f 9 > "$results_dir""$sample_name".InsertSizeMetrics.txt
#samtools stats "$DIR""$sample_name".alignments.bam > "$results_dir""$sample_name".summary.txt
#bedtools bamtobed -i "$DIR""$sample_name".alignments.bam > "$results_dir""$sample_name".bed

#bedtools intersect -v -a "$results_dir""$sample_name".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_intron.bed > "$results_dir""$sample_name".rm.exon.intron.hg19.bed

#bedtools window -l 4500 -r 0 -a "$results_dir""$sample_name".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > "$results_dir""$sample_name".gene.downstream.hg19.bed

#bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.downstream.hg19.strand.based.bed

#bedtools window -l 0 -r 4500 -a /media/H_driver/2016/Ramin_azhang/Annotation/knowngene.bed -b "$results_dir""$sample_name".rm.exon.intron.bed > "$results_dir""$sample_name".gene.downstream.a.as.gene.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".rm.exon.intron.hg19.bed -l 0 -r 4500 -sw   > "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed

#awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.downstream.count.hg19.strand.based.txt  

#/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"$sample_name".alignments.bam --min 100 --max 200 --out /media/H_driver/2016/Ramin_azhang/Results/"$sample_name".alignments.split.bam  

#EOF


cat > ~/Code/Run_"$sample_name"_3UTR_bam2bed.sh <<EOF
#Convert bam file to bed file
bedtools bamtobed -i "$DIR""$sample_name".alignments.bam > "$results_dir""$sample_name".bed
EOF

cat > ~/Code/Run_"$sample_name"_3UTR_rm_exon_intron.sh <<EOF
#Remove reads in exons and introns
bedtools intersect -v -a "$results_dir""$sample_name".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_intron.bed > "$results_dir""$sample_name".rm.exon.intron.hg19.bed
EOF

cat > ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads.sh <<EOF
#For each gene, find the reads that fall in the 4500bp region of downstream of each gene
#bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.downstream.hg19.strand.based.bed

#For each gene, find the reads that fall in the 4500bp region of downstream of each gene
#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".rm.exon.intron.hg19.bed -l 0 -r 4500 -sw   > "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed

#For each gene, find the reads that fall in the 4500bp region of upstream of each gene
#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".rm.exon.intron.hg19.bed -l 4500 -r 0 -sw   > "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.bed

#For each gene, find the reads that fall in the 4500bp region of upstream and downstream of each gene
#bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.downstream.hg19.strand.based.bed

EOF

cat >  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count.sh <<EOF

#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
#awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.downstream.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
awk -F '\t' '{print \$10}' "$results_dir""$sample_name".gene.downstream.hg19.strand.based.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.and.downstream.count.hg19.strand.based.txt  

EOF

#sh ~/Code/Run_"$sample_name"_3UTR.sh&
#wait
#sh ~/Code/Run_"$sample_name"_3UTR_strand_based.sh&

#sh ~/Code/Run_"$sample_name"_3UTR_bam2bed.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_rm_exon_intron.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads.sh&
#wait

sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count.sh&

done
