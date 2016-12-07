
#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"ATP5A1".alignments.bam|cut -f 9 > "/media/H_driver/PJ/Alignment/75_S1_STAR/Results/""ATP5A1".InsertSizeMetrics.txt
#samtools stats "/media/H_driver/PJ/Alignment/75_S1_STAR/*""ATP5A1".alignments.bam > "/media/H_driver/PJ/Alignment/75_S1_STAR/Results/""ATP5A1".summary.txt

/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in "/media/H_driver/PJ/Alignment/75_S1_STAR/*""ATP5A1".sorted.bam --min 100 --max 200 --out "/media/H_driver/PJ/Alignment/75_S1_STAR/Results/""ATP5A1".sorted.split.bam  

