
#samtools view -f66 /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"2016-02-10-emp2_Dox".alignments.bam|cut -f 9 > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".InsertSizeMetrics.txt
#samtools stats "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/""2016-02-10-emp2_Dox".alignments.bam > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".summary.txt
#bedtools bamtobed -i "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/""2016-02-10-emp2_Dox".alignments.bam > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".bed

#bedtools intersect -v -a "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_intron.bed > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".rm.exon.intron.hg19.bed

#bedtools window -l 4500 -r 0 -a "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.hg19.bed

bedtools window -l 0 -r 4500 -a "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".rm.exon.intron.hg19.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -sw > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.hg19.strand.based.bed

#bedtools window -l 0 -r 4500 -a /media/H_driver/2016/Ramin_azhang/Annotation/knowngene.bed -b "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".rm.exon.intron.bed > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.a.as.gene.bed

#bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".rm.exon.intron.hg19.bed -l 0 -r 4500 -sw   > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.a.as.gene.strand.based.hg19.bed

#awk -F '\t' '{print $4}' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.a.as.gene.strand.based.hg19.bed | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/""2016-02-10-emp2_Dox".gene.downstream.count.hg19.strand.based.txt  

#/media/aiminyan/DATA/biotoolbox-extra/scripts/split_bam_by_isize.pl --in /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"2016-02-10-emp2_Dox".alignments.bam --min 100 --max 200 --out /media/H_driver/2016/Ramin_azhang/Results/"2016-02-10-emp2_Dox".alignments.split.bam  

