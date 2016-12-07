
#infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -i "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/""hela_wt1".alignments.bam > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1"_check_report.txt

#echo "hela_wt1" >> "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"check_summary_report.txt
#infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -i "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/""hela_wt1".alignments.bam >> "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"check_summary_report.txt
 
#/media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed 

echo "hela_wt1" >> "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"check_summary_report_2.txt
infer_experiment.py -r /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed  -i "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/""hela_wt1".alignments.bam >> "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"check_summary_report_2.txt

