awk '{print $4,$1,$2,$3,$6}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.bed | awk -v OFS="\t" '$1=$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.SAF

awk -F"\t" '{ $3 = ($3 <0 ? 0 : $3) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_4.SAF > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF

