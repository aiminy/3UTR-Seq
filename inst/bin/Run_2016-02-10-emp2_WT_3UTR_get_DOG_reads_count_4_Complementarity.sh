
#For each gene, get count of the reads that fall in the 4500bp region of downstream of each gene
awk -F '\t' '$6!=$18' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.downstream.a.as.gene.strand.based.hg19.bed | awk '{print $4}' | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.downstream.com.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of upstream of each gene
awk -F '\t' '$6!=$18' "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.upstream.a.as.gene.strand.based.hg19.bed | awk '{print $4}' | sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.upstream.com.count.hg19.strand.based.txt

#For each gene, get count of the reads that fall in the 4500bp region of downstream and upstream of each gene 
#awk -F '\t' '!='  "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.upstream.and.downstream.hg19.strand.based.bed | awk '{print $10}'| sort | uniq -c | sort -nr > "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp2_WT".gene.upstream.and.downstream.com.count.hg19.strand.based.txt  

