
#/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in
#/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/2016-02-10-emp1_Dox.alignments.bam
#--min 100 --max 200
#--out /media/H_driver/2016/Ramin_azhang/Results/2016-02-10-emp1_Dox.alignments.split.bam

bam2gff_bed.pl
--in /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/2016-02-10-emp1_Dox.alignments.bam
--bed
--out /media/H_driver/2016/Ramin_azhang/Results/RNA_seq/2016-02-10-emp1_Dox.alignments.bed
