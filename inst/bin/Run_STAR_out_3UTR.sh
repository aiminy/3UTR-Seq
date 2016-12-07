
#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"STAR_out".alignments.bam|cut -f 9 > "/media/H_driver/PJ/Alignment/75_S1_STAR/Results/""STAR_out".InsertSizeMetrics.txt
#samtools stats "/media/H_driver/PJ/Alignment/75_S1_STAR/*""STAR_out".alignments.bam > "/media/H_driver/PJ/Alignment/75_S1_STAR/Results/""STAR_out".summary.txt

/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in /media/H_driver/PJ/Alignment/75_S1_STAR/STAR_out.sorted.bam --min 100 --max 200 --out /media/H_driver/PJ/Alignment/75_S1_STAR/Results/STAR_out.sorted.split.bam
