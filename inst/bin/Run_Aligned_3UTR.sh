
#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"Aligned".alignments.bam|cut -f 9 > "Results/""Aligned".InsertSizeMetrics.txt
#samtools stats """Aligned".alignments.bam > "Results/""Aligned".summary.txt
#bedtools bamtobed -i """Aligned".alignments.bam > "Results/""Aligned".bed
bedtools intersect -v -a "Results/""Aligned".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/intron.bed
#/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"Aligned".alignments.bam --min 100 --max 200 --out /media/H_driver/2016/Ramin_azhang/Results/"Aligned".alignments.split.bam  

