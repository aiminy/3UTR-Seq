tophat -G ~/genes.gtf -p 4 -o "trimmed_2015-11-20-WT10_S5_combined4f"_tophat_out ~/mm10_index_bt2/genome "trimmed_2015-11-20-WT10_S5_combined4f.fastq"
mv "trimmed_2015-11-20-WT10_S5_combined4f"_tophat_out/accepted_hits.bam "trimmed_2015-11-20-WT10_S5_combined4f".bam
