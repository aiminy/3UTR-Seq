#Remove reads in exons and introns
bedtools intersect -v -a "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_dox2".bed -b /media/H_driver/2016/Ramin_azhang/Annotation/hg19_exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/hg19_intron.bed > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_dox2".rm.exon.intron.hg19.bed
