#!/bin/bash

#Usage: sh /home/aiminyan/Code/Process3UTR_3.sh input_files_dir output_dir_name
#Example: sh /home/aiminyan/Code/Process3UTR_3.sh /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/ Results4Check2

DIR="$1"

#results_dir="$1Results4Feature/"

results_dir="$1$2/"

#results_dir="$2/"

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

cat > ~/Code/Run_check_"$sample_name".sh <<EOF

#infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -i "$DIR""$sample_name".alignments.bam > "$results_dir""$sample_name"_check_report.txt

#echo "$sample_name" >> "$results_dir"check_summary_report.txt
#infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -i "$DIR""$sample_name".alignments.bam >> "$results_dir"check_summary_report.txt
 
#/media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed 

echo "$sample_name" >> "$results_dir"check_summary_report_2.txt
infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed  -i "$DIR""$sample_name".alignments.bam >> "$results_dir"check_summary_report_2.txt

EOF

cat > ~/Code/Run_process_hg19_gene.sh <<EOF

#nl /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed

#awk '\$7=="+"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed

#awk '\$7=="-"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed

awk '{\$8=\$8-\$3;\$9=\$9-\$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed | awk -v OFS="\t" '\$1=\$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed
awk '{\$8=\$8-\$3;\$9=\$9-\$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed | awk -v OFS="\t" '\$1=\$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed

awk '{\$3 = \$4; \$4=\$4+4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_downstream.bed

awk '{\$4 = \$3; \$3=\$3-4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_downstream.bed

cat /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_downstream.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_downstream.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream.bed
 
sort -nk1 /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted.bed

awk -v OFS="\t" '\$1=\$1' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_2.bed

cut -f2- /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_3.bed

awk -F"\t" '{ \$2 = (\$2 <0 ? 0 : \$2) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_3.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_33.bed  

awk '{\$7=\$2+\$7;\$8=\$2+\$8;\$10=1;\$11=\$3-\$2",";\$12="0,";print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_33.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed  

#awk -F"\t" '{ \$2 = (\$2 <0 ? 0 : \$2) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.bed

EOF


cat > ~/Code/Run_"$sample_name"_3UTR_count_plus_minus.sh <<EOF
#FPKM_count_2.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_3.bed -d '+-,-+' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.minus.and.minus.plus.2.txt

FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed -d '++,--' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.plus.and.minus.minus.4.DoGs

FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed -d '+-,-+' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.minus.and.minus.plus.4.DoGs

FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_7.bed -d '++,--' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.plus.and.minus.minus.4.gene.plus.DoGs

FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_7.bed -d '+-,-+' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.minus.and.minus.plus.4.gene.plus.DoGs

#FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.bed -d '++,--' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.plus.and.minus.minus.4

#FPKM_count.py -i "$DIR""$sample_name".alignments.bam -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.bed -d '+-,-+' -q 0 -s 1 -o "$results_dir""$sample_name".gene.downstream.count.plus.minus.and.minus.plus.4

EOF



cat > ~/Code/Run_process_hg19_gene_4_gene_plus_DoGs.sh <<EOF

#nl /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed

#awk '\$7=="+"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed

#awk '\$7=="-"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed

awk '{\$8=\$8-\$3;\$9=\$9-\$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed | awk -v OFS="\t" '\$1=\$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed

awk '{\$8=\$8-\$3;\$9=\$9-\$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed | awk -v OFS="\t" '\$1=\$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed

awk '{\$4=\$4+4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_plus_downstream.bed

awk '{\$3=\$3-4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_plus_downstream.bed

cat /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_plus_downstream.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_plus_downstream.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream.bed
 
sort -nk1 /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted.bed

awk -v OFS="\t" '\$1=\$1' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_2.bed

cut -f2- /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_3.bed

awk '{\$10=1;\$11=\$3-\$2",";\$12="0,";print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_3.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_4.bed  

awk -F"\t" '{ \$2 = (\$2 <0 ? 0 : \$2) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_4.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_5.bed

awk '{\$11=\$3-\$2",";print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_5.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_6.bed  

awk '{\$7=\$2+\$7;\$8=\$2+\$8;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_6.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_7.bed

EOF


cat > ~/Code/Run_process_hg19_gene_to_SAF.sh <<EOF
awk '{print \$4,\$1,\$2,\$3,\$6}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.SAF

awk -F"\t" '{ \$3 = (\$3 <0 ? 0 : \$3) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.SAF > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF

EOF

cat > ~/Code/Run_"$sample_name"_3UTR_count_use_feature_count.sh <<EOF
#/home/aiminyan/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF -F SAF -s 1 -o "$results_dir""$sample_name".gene.downstream.count.fw.FC.txt "$DIR""$sample_name".alignments.bam

#/home/aiminyan/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF -F SAF -s 2 -o "$results_dir""$sample_name".gene.downstream.count.rv.FC.txt "$DIR""$sample_name".alignments.bam

/home/aiminyan/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF -F SAF -s 1 -M --fraction -o "$results_dir""$sample_name".gene.downstream.count.fw.m.FC.txt "$DIR""$sample_name".alignments.bam

/home/aiminyan/subread-1.5.0-p3-Linux-x86_64/bin/featureCounts -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF -F SAF -s 2 -M --fraction -o "$results_dir""$sample_name".gene.downstream.count.rv.m.FC.txt "$DIR""$sample_name".alignments.bam

EOF

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
bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".rm.exon.intron.hg19.bed -l 0 -r 4500 -sw   > "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed

#For each gene, find the reads that fall in the 4500bp region of upstream of each gene
bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".rm.exon.intron.hg19.bed -l 4500 -r 0 -sw   > "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.bed

#For each gene, find the reads that fall in the 4500bp region of upstream and downstream of each gene
bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.bed

EOF

cat >  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count.sh <<EOF

#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.downstream.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
awk -F '\t' '{print \$10}' "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.and.downstream.count.hg19.strand.based.txt  

EOF

cat > ~/Code/Run_"$sample_name"_3UTR_get_DOG_use_total_reads.sh <<EOF
#For each gene, find the reads that fall in the 4500bp region of downstream of each gene
bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".bed -l 0 -r 4500 -sw   > "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.total2.bed

#For each gene, find the reads that fall in the 4500bp region of upstream of each gene
bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".bed -l 4500 -r 0 -sw   > "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.total2.bed

#For each gene, find the reads that fall in the 4500bp region of upstream and downstream of each gene
bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.total2.bed

EOF


cat >  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_total.sh <<EOF

#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.total.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.downstream.count.hg19.strand.based.total.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '{print \$4}' "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.total.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.count.hg19.strand.based.total.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
awk -F '\t' '{print \$10}' "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.total.bed | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.and.downstream.count.hg19.strand.based.total.txt  

EOF


cat >  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_4_Complementarity.sh <<EOF

#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
awk -F '\t' '\$6!=\$18' "$results_dir""$sample_name".gene.downstream.a.as.gene.strand.based.hg19.bed | awk '{print \$4}' | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.downstream.com.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '\$6!=\$18' "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.bed | awk '{print \$4}' | sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.com.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
#awk -F '\t' '$6!=$18'  "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.bed | awk '{print \$10}'| sort | uniq -c | sort -nr > "$results_dir""$sample_name".gene.upstream.and.downstream.com.count.hg19.strand.based.txt  

EOF

cat >  ~/Code/Run_process_bed_4_using_coverageBed.sh <<EOF
awk -F"\t" '{OFS="\t";print \$1,\$2,\$3,\$4,\$5,\$6}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed | awk -v OFS="\t" '\$1=\$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4_coverageBed.bed
EOF

cat >  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_using_coverageBed.sh <<EOF
intersectBed -a "$DIR""$sample_name".alignments.bam -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4_coverageBed.bed -wb -s -bed > "$results_dir""$sample_name".gene.DoGs.using.coverageBed.2.txt
EOF

cat > ~/Code/Check_if_DoGs_overlap_with_transcript.sh <<EOF

#awk -F"\t" '{print \$1,\$2,\$3,\$4,\$5,\$6}' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed

#awk -F"\t" '{print \$2,\$3,\$4,\$1,\$4-\$3,\$5}' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs.bed

#awk '\$3>\$2' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed 

#/home/aiminyan/bedops/applications/bed/sort-bed/bin/sort-bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body_sorted.bed

#/home/aiminyan/bedops/applications/bed/sort-bed/bin/sort-bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_sorted.bed

#/home/aiminyan/bedops/applications/bed/bedops/bin/bedops --difference /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_sorted.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body_sorted.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_on_DoGs_but_not_on_transcripts.bed

#/home/aiminyan/bedops/applications/bed/bedops/bin/bedops --not-element-of 1 /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_sorted.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body_sorted.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_not_overlap_with_any_transcript.bed


#/home/aiminyan/bedops/applications/bed/bedextract/bin/bedextract/ /media/H_driver/2016/Ramin_azhang/Annotation/hg19_on_DoGs_but_not_on_transcripts.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_sorted.bed

#bedtools intersect -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_1.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -sw > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_3.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_4.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -u > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_u.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -c > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_c.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -v > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_v.bed

bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed  -l 0 -r 0 > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_GoGs_DoGs_overlap.bed


#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -sw -sm > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_4_same_st.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed  -l 0 -r 0 -sw -Sm > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_4_diff_st.bed

#bedtools intersect -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs.bed >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_DoGs_overlap_2.bed

EOF


cat > ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_use_DoGs_only.sh <<EOF
#For each gene, find the reads that fall in the 4500bp region of downstream of each gene
bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed -b "$results_dir""$sample_name".bed -l 0 -r 0 -sw   > "$results_dir""$sample_name".DoGs.only.strand.based.hg19.total2.bed

#For each gene, find the reads that fall in the 4500bp region of upstream of each gene
#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "$results_dir""$sample_name".bed -l 4500 -r 0 -sw   > "$results_dir""$sample_name".gene.upstream.a.as.gene.strand.based.hg19.total2.bed

#For each gene, find the reads that fall in the 4500bp region of upstream and downstream of each gene
#bedtools window -l 0 -r 4500 -a "$results_dir""$sample_name".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "$results_dir""$sample_name".gene.upstream.and.downstream.hg19.strand.based.total2.bed

EOF


#sh  ~/Code/Run_check_"$sample_name".sh&
#wait

#sh  ~/Code/Run_process_hg19_gene.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_count_plus_minus.sh&
#wait

#sh ~/Code/Run_process_hg19_gene_to_SAF.sh&
#wait

#sh ~/Code/Run_process_hg19_gene_4_gene_plus_DoGs.sh
#wait

#sh  ~/Code/Run_"$sample_name"_3UTR_count_use_feature_count.sh&

#sh ~/Code/Run_"$sample_name"_3UTR.sh&
#wait
#sh ~/Code/Run_"$sample_name"_3UTR_strand_based.sh&

#sh ~/Code/Run_"$sample_name"_3UTR_bam2bed.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_rm_exon_intron.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count.sh&

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_use_total_reads.sh&
#wait

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_total.sh&

#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_4_Complementarity.sh&

#sh  ~/Code/Run_process_bed_4_using_coverageBed.sh&
#sh ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_count_using_coverageBed.sh&

#sh ~/Code/Check_if_DoGs_overlap_with_transcript.sh&

sh  ~/Code/Run_"$sample_name"_3UTR_get_DOG_reads_use_DoGs_only.sh

done
