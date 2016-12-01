awk ' \
    BEGIN \
    { \
        binSize = 1000000; \
        binIdx = 0; \
    } \
    { \
        chr = $1; \
        start = $2; \
        stop = $3; \
        for (binStart = start; binStart < (stop - binSize); binStart += binSize) { \
            print chr"\t"binStart"\t"(binStart + binSize)"\tbin-"binIdx; \
            binIdx++; \
        } \
    }' /media/aiminyan/DATA/Yang_ChipSeq_RNA_Seq/Bin_data/chrList.bed \
q> /media/aiminyan/DATA/Yang_ChipSeq_RNA_Seq/Bin_data/myBins.bed
