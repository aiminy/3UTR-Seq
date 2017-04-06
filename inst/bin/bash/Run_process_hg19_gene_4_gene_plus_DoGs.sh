
#nl /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed

#awk '$7=="+"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed

#awk '$7=="-"' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_2.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed

awk '{$8=$8-$3;$9=$9-$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos.bed | awk -v OFS="\t" '$1=$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed

awk '{$8=$8-$3;$9=$9-$3;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg.bed | awk -v OFS="\t" '$1=$1' >  /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed

awk '{$4=$4+4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_plus_downstream.bed

awk '{$3=$3-4500;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_plus_downstream.bed

cat /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_plus_downstream.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_neg_plus_downstream.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream.bed
 
sort -nk1 /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream.bed  > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted.bed

awk -v OFS="\t" '$1=$1' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_2.bed

cut -f2- /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_2.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_3.bed

awk '{$10=1;$11=$3-$2",";$12="0,";print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_3.bed | awk -v OFS="\t" '$1=$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_4.bed  

awk -F"\t" '{ $2 = ($2 <0 ? 0 : $2) } 1' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_4.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_5.bed

awk '{$11=$3-$2",";print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_5.bed | awk -v OFS="\t" '$1=$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_6.bed  

awk '{$7=$2+$7;$8=$2+$8;print}' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_6.bed | awk -v OFS="\t" '$1=$1' > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_plus_downstream_sorted_7.bed

