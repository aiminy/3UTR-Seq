
#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
awk -F '\t' '{print $4}' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.downstream.a.as.gene.strand.based.hg19.total.bed | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.downstream.count.hg19.strand.based.total.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '{print $4}' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.upstream.a.as.gene.strand.based.hg19.total.bed | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.upstream.count.hg19.strand.based.total.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
awk -F '\t' '{print $10}' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.upstream.and.downstream.hg19.strand.based.total.bed | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""hela_wt1".gene.upstream.and.downstream.count.hg19.strand.based.total.txt  

