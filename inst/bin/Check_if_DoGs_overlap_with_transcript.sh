
#awk -F"\t" '{print $1,$2,$3,$4,$5,$6}' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_transcript_body.bed

#awk -F"\t" '{print $2,$3,$4,$1,$4-$3,$5}' OFS="\t" /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene_pos_neg_downstream_sorted_5.SAF > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs.bed

#awk '$3>$2' /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs.bed > /media/H_driver/2016/Ramin_azhang/Annotation/hg19_DoGs_2.bed 

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

