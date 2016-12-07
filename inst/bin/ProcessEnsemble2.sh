#awk -F "\t" '{OFS="\t"; $1 = "chr"$1; print}' /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.gtf | awk -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.gtf
#sort -k1,1 -k4,4n /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.gtf > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf

#awk -F "\t" '{OFS="\t"; $1 = "chr"$1; print}' /media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.gtf | awk -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' > /media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.gtf
sort -k1,1 -k4,4n /media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.gtf > /media/H_driver/Aimin_project/GTF_Files/Homo_sapiens.GRCh38.84.processed.sorted.2.gtf



