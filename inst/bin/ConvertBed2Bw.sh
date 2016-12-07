sort -k1,1n -k2,2n GSE19013_MCF7_E2-ethl_chr7_17_norm.bed > GSE19013_MCF7_E2-ethl_chr7_17_norm_sorted.bed 
awk '{printf "%s\t%d\t%d\t%2.3f\n" , "chr"$1,$2,$3,$4}' GSE19013_MCF7_E2-ethl_chr7_17_norm_sorted.bed > GSE19013_MCF7_E2-ethl_chr7_17_norm_sorted.bedgraph
/home/aiminyan/kentUtils/bin/linux.x86_64/bedGraphToBigWig GSE19013_MCF7_E2-ethl_chr7_17_norm_sorted.bedgraph ../Data_Chip-Seq/hg19.genome GSE19013_MCF7_E2-ethl_chr7_17_norm_sorted.bw
