#awk -F "\t" '{OFS="\t"; $1 = "chr"$1; print}' /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.gtf | awk -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.gtf
#sort -k1,1 -k4,4n /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.gtf > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.DEXSeq.cleaned.sorted.gtf

awk -F "\t" '{OFS="\t"; $1 = "chr"$1; print}' /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.gtf | awk -F"\t" '{OFS="\t"; if($1=="chrMT") $1="chrM"; print}' > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.gtf
sort -k1,1 -k4,4n /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.gtf > /media/DATA/mus_musculus/Mus_musculus.GRCm38.83.processed.sorted.gtf



