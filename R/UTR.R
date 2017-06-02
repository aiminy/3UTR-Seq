#' AdjustBatch
#'
#' @param re.rMAT.com.total.test.2
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/'
#'
#' re.3UTR<-AdjustBatch(re.rMAT.com.total.test.2,dir.name,'By_P','By_Gene','Interested_gene_only','Interested_gene_only_sig')
#'
#' re.3UTR.one.side<-AdjustBatch(re.rMAT.com.total.test.2,dir.name,'By_P_one_side','By_Gene_one_side',
#' 'Interested_gene_only_one_side','Interested_gene_only_sig_one_side')
#'
#'
#' re.3UTR.one.side.norm<-AdjustBatch(re.rMAT.com.total.test.2,DE.norm.with.rpkm.norm.2,dir.name,'By_P_one_side_norm_2','By_Gene_one_side_norm_2',
#' 'Interested_gene_only_one_side_norm_2','Interested_gene_only_sig_one_side_norm_2')
#'
AdjustBatch <- function(re.rMAT.com.total.test.2, DE.norm.with.rpkm.norm.2,
    dir.name, out_by_p, out_by_gene, out_gene_interested, out_gene_interested_sig)
    {

    # names(re.rMAT.com.total.test.2)

    df.NT <- counts(re.rMAT.com.total.test.2$dds)
    colData <- colData(re.rMAT.com.total.test.2$dds)

    cell <- factor(rep(c("emp", "hela"), c(2, 2)))
    cell = rep(cell, 2)

    colData.1 <- data.frame(condition = factor(rep(c("Dox", "WT"), c(4, 4))),
        cell = cell)

    # Use original design, not adjusted by any batch factor
    dds <- DESeqDataSetFromMatrix(df.NT, colData.1, formula(~condition))
    colData(dds)

    print(design(dds))
    dds <- DESeq(dds)
    res.dds <- results(dds)
    head(res.dds)
    print(sizeFactors(estimateSizeFactors(dds)))
    print(head(counts(dds, normalized = TRUE)))

    # Check hidden batch effects dds <- estimateSizeFactors(dds)
    dds <- DESeqDataSetFromMatrix(df.NT, colData.1, formula(~condition))
    dat <- counts(dds)

    # idx <- rowMeans(dat) > 1 dat <- dat[idx,]

    mod <- model.matrix(~condition, colData(dds))
    mod0 <- model.matrix(~1, colData(dds))
    svseq.NT <- svaseq(dat, mod, mod0)

    par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
    stripchart(svseq.NT$sv[, 1] ~ dds$cell, vertical = TRUE, main = "SV1")
    abline(h = 0)
    stripchart(svseq.NT$sv[, 2] ~ dds$cell, vertical = TRUE, main = "SV2")
    abline(h = 0)

    # Adjust by surrogate variables 1 and 2
    ddssva <- dds
    ddssva$SV1 <- svseq.NT$sv[, 1]
    ddssva$SV2 <- svseq.NT$sv[, 2]
    design(ddssva) <- ~SV1 + SV2 + condition
    ddssva <- DESeq(ddssva)
    print(sizeFactors(estimateSizeFactors(ddssva)))
    print(head(counts(ddssva, normalized = TRUE)))
    res.ddssva.2.sva.adjusted <- results(ddssva)
    head(res.ddssva.2.sva.adjusted)
    summary(res.ddssva.2.sva.adjusted)

    # Adjusted by cell type

    dds.cell <- DESeqDataSetFromMatrix(df.NT, colData.1, design = ~cell + condition)

    # design(dds.cell) <- ~ cell + condition
    print(design(dds.cell))

    dds.cell <- DESeq(dds.cell)
    res.dds.cell <- results(dds.cell)

    re.condition <- results(dds.cell, contrast = c("condition", "WT", "Dox"))
    head(re.condition)
    summary(re.condition)

    re.cell <- results(dds.cell, contrast = c("cell", "emp", "hela"))
    summary(re.cell)

    head(res.dds.cell)
    print(sizeFactors(estimateSizeFactors(dds.cell)))
    print(head(counts(dds.cell, normalized = TRUE)))

    # dds.cell.2 <- DESeq(dds.cell.2) res.dds.cell.2<- results(dds.cell.2)
    # head(res.dds.cell.2) summary(res.dds.cell.2) summary(res.dds.cell) #Use
    # limma based method
    # #dat.adjusted<-cbind(dat[,c(3,4,7,8)]*(2^re.cell$log2FoldChange),dat[,c(3,4,7,8)])
    # dds.ad <- DESeqDataSetFromMatrix(df.NT,colData.1,design=~ cell+condition)
    # dds.ad <- estimateSizeFactors(dds.ad) norm.counts <- counts(dds.ad,
    # normalized=TRUE) log.norm.counts <- log2(norm.counts + 1)
    # y<-log.norm.counts batch <-
    # c('emp','emp','hela','hela','emp','emp','hela','hela')
    # y2<-removeBatchEffect(y, batch) #par(mfrow=c(1,2))
    # boxplot(as.data.frame(y),main='Original')
    # boxplot(as.data.frame(y2),main='Batch corrected') f.st134.frma
    # <-factor(colData.1$condition) design.st134.frma <-
    # model.matrix(~0+f.st134.frma) colnames(design.st134.frma) <-
    # levels(f.st134.frma) fit.st134.all.probes.frma <- lmFit(y2,
    # design.st134.frma) cont.matrix.st134.frma <-
    # makeContrasts(stWTDox='WT-Dox',levels=design.st134.frma)
    # fit2.st134.all.probes.frma <- contrasts.fit(fit.st134.all.probes.frma,
    # cont.matrix.st134.frma) fit2.st134.all.probes.frma <-
    # eBayes(fit2.st134.all.probes.frma)
    # TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(y2)[1],sort.by='p')
    # cutoff <- 0.3 wtResCont <- decideTests(fit2.st134.all.probes.frma, p.value
    # = cutoff, method = 'global') summary(wtResCont)
    # length(which(TopTableSt34.all.probes.frma$adj.P.Val<0.3&TopTableSt34.all.probes.frma$logFC<0))
    # length(which(TopTableSt34.all.probes.frma$adj.P.Val<0.3&TopTableSt34.all.probes.frma$logFC>0))
    # #matrix of SVs svseq.NT$sv #number of SVs svseq.NT$n.sv save(svseq.NT,
    # file='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/svseq.NT.RData')
    # colData.sva <- cbind(colData,svseq.NT$sv) colData.sva <-
    # data.frame(condition=factor(rep(c('Dox', 'WT'), c(4, 4)),svseq.NT$sv)) #
    # preprocess data colData.sva <- data.frame(condition=factor(rep(c('Dox',
    # 'WT'), c(4, 4))), svseq.NT$sv) colnames(colData.sva)[-1] = paste0('SV',
    # 1:svseq.NT$sv) ncol <- dim(svseq.NT$sv)[2] SV1<-paste(paste0('SV',
    # 1:svseq.NT$sv), collapse = '+') design.sva.1 <- formula(paste0('~ SV1
    # +','condition')) design.sva.2 <- formula(paste0('~ condition +','SV1'))
    # #design.sva <- formula(paste0('~ paste(paste0('SV', 1:svseq.NT$sv),
    # collapse = '+') +', condition )) colData.sva.1<-colData.sva[,c(2,1)]
    # ddssva.NT <- DESeqDataSetFromMatrix(df.NT,colData.sva, design.sva.1)
    # rownames(colData(ddssva.NT)) <- colnames(df.NT) ddssva.NT <-
    # DESeq(ddssva.NT) res.sva.NT <- results(ddssva.NT) head(res.sva.NT)
    # ddssva.NT.1 <- DESeqDataSetFromMatrix(df.NT,colData.sva.1, design.sva.1)
    # rownames(colData(ddssva.NT.1)) <- colnames(df.NT) ddssva.NT.1 <-
    # DESeq(ddssva.NT.1) res.sva.NT.1 <- results(ddssva.NT.1) head(res.sva.NT.1)
    # ddssva.NT.2 <- DESeqDataSetFromMatrix(df.NT,colData.sva, design.sva)
    # rownames(colData(ddssva.NT.2)) <- colnames(df.NT) ddssva.NT.2 <-
    # DESeq(ddssva.NT.2) res.sva.NT.2 <- results(ddssva.NT.2) head(res.sva.NT)

    MatchGene <- function(res.sva.NT, re.rMAT.com.total.test.2)
    {
        temp <- re.rMAT.com.total.test.2$DE[, c(1, 10:24)]
        rownames(temp) <- temp[, 1]
        temp2 <- merge(as.data.frame(res.sva.NT), temp, by = 0, sort = FALSE)
        re.FC <- temp2
        re.FC.sorted <- re.FC[order(re.FC$pvalue), ]
        return(re.FC.sorted)
    }


    re.FC.sorted <- MatchGene(res.ddssva.2.sva.adjusted, re.rMAT.com.total.test.2)

    DE.norm.with.rpkm.norm.3 <- apply(DE.norm.with.rpkm.norm.2, 2, as.numeric)

    rownames(DE.norm.with.rpkm.norm.3) <- rownames(DE.norm.with.rpkm.norm.2)

    DE.rpkm.norm <- cbind(DE.norm.with.rpkm.norm.3, DE.norm.with.rpkm.norm.3[,
        9:16] * (1000/4500))

    rownames(re.FC.sorted) = re.FC.sorted$Row.names

    re.FC.sorted <- merge(re.FC.sorted[, -1], DE.rpkm.norm, by = 0, sort = FALSE)

    p.one.sided <- unlist(lapply(re.FC.sorted$stat, function(x)
    {
        convert.z.score(x, two.side = F)
    }))

    p.one.sided.adjust <- p.adjust(p.one.sided, method = "BH")

    re.FC.sorted <- cbind(re.FC.sorted, p.one.sided, p.one.sided.adjust)

    re.FC.sorted.sorted.by.gene <- re.FC.sorted[order(re.FC.sorted$gene), ]

    write.csv(re.FC.sorted, file = paste0(dir.name, "3UTR_DE_", out_by_p, ".csv"),
        row.names = FALSE, quote = FALSE)

    write.csv(re.FC.sorted.sorted.by.gene, file = paste0(dir.name, "3UTR_DE_",
        out_by_gene, ".csv"), row.names = FALSE, quote = FALSE)

    target <- re.FC.sorted[which(re.FC.sorted$pvalue < 0.05 & re.FC.sorted$log2FoldChange <
        0), ]$gene

    # which(target %in% genes.interested[,1])

    re.genes.interested <- re.FC.sorted.interested <- re.FC.sorted[which(re.FC.sorted$gene %in%
        genes.interested[, 1]), ]

    re.genes.interested.sig <- re.FC.sorted.interested[which(re.FC.sorted.interested$pvalue <
        0.05 & re.FC.sorted.interested$log2FoldChange < 0), ]

    write.csv(re.genes.interested, file = paste0(dir.name, "3UTR_DE_", out_gene_interested,
        ".csv"), row.names = FALSE, quote = FALSE)

    write.csv(re.genes.interested.sig, file = paste0(dir.name, "3UTR_DE_", out_gene_interested_sig,
        ".csv"), row.names = FALSE, quote = FALSE)

    # re.com.total.pos<-re.FC.sorted[which(re.FC.sorted$`FoldChange(WT/Dox)`<1),]
    # re.com.total.neg<-re.FC.sorted[which(re.FC.sorted$`FoldChange(WT/Dox)`>=1),]

    re <- list(re.FC.sorted = re.FC.sorted, re.FC.sorted.sorted.by.gene = re.FC.sorted.sorted.by.gene,
        re.genes.interested = re.genes.interested, re.genes.interested.sig = re.genes.interested.sig)

    return(re)
}
#' CheckBatch
#'
#' @param df.NT
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' df.NT<-re.DESeq.DoGs.plus.minus.gene.interested[,c(1,9:16)]
#' rownames(df.NT)<-df.NT[,1]
#' df.NT.2<- df.NT[,-1]
#'
#' CheckBatch(dir.name,'Btach_using_DoGs',df.NT.2)
#'
CheckBatch <- function(dir.name, output_bacth_check, df.NT)
{

    cell <- factor(rep(c("emp", "hela"), c(2, 2)))
    cell = rep(cell, 2)

    colData.1 <- data.frame(condition = factor(rep(c("Dox", "WT"), c(4, 4))),
        cell = cell)

    # Check hidden batch effects dds <- estimateSizeFactors(dds)
    dds <- DESeqDataSetFromMatrix(df.NT, colData.1, formula(~condition))
    dat <- counts(dds)

    # idx <- rowMeans(dat) > 1 dat <- dat[idx,]

    mod <- model.matrix(~condition, colData(dds))
    mod0 <- model.matrix(~1, colData(dds))
    svseq.NT <- svaseq(dat, mod, mod0)

    png(paste0(dir.name, output_bacth_check, ".png"))
    par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
    stripchart(svseq.NT$sv[, 1] ~ dds$cell, vertical = TRUE, main = "SV1")
    abline(h = 0)
    stripchart(svseq.NT$sv[, 2] ~ dds$cell, vertical = TRUE, main = "SV2")
    abline(h = 0)
    dev.off()

}
#' Compare2Pos
#'
#' @param in1
#' @param in2
#'
#' @return
#' @export Compare2Pos
#'
#' @examples
#'
#' #real data
#' in1='/Volumes/Bioinformatics$/2016/Ramin_azhang/ResultsDoGsOnlyRmOne/RealData/3UTR_DE_DoGs_adjust_by_batch_all.csv'
#'
#' #permutated data
#'
#' in1='/Volumes/Bioinformatics$/2016/Ramin_azhang/ResultsDoGsOnlyRmOne/Counts4DoGsOnlyRmOnePermutation/3UTR_DE_DoGs_adjust_by_batch_all.csv'
#' threshold=0.05
#' in2='/Volumes/Bioinformatics$/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/final_list.csv'
#'
#' Compare2Pos(in1,threshold,in2)
#'
Compare2Pos <- function(in1, threshold, in2)
{

    a <- read.table(in1, header = T, sep = ",")
    b <- read.table(in2, header = F, sep = ",")

    print(head(a))

    print(head(b))

    g1 <- unique(as.character(a[which(a$padj <= threshold), 22]))
    g2 <- unique(as.character(b[, 1]))

    g11 <- g1[-which(is.na(g1))]

    g22 <- g2

    print(g11)

    print(g22)

    re <- list(R = g11, pos = g22, a = a, b = b)

    names(re)[1] <- paste0("R", "_", threshold)

    venn.plot <- venn.diagram(x = re[c(1, 2)], filename = paste0(dirname(in1),
        "/", names(re)[1], names(re)[2], "_overlap_venn.tiff"), height = 3000,
        width = 3500, resolution = 1000, col = "black", lty = "dotted", lwd = 1,
        fill = c("red", "blue"), alpha = 0.5, label.col = c(rep("black", 3)),
        cex = 0.5, fontfamily = "serif", fontface = "bold", cat.col = c("red",
            "blue"), cat.cex = 0.5, cat.pos = 0.5, cat.dist = 0.05, cat.fontfamily = "serif")

    return(re)

}
CompareToFeature.R <- function(Re.unadjusted.adjusted)
{

    length(unique(as.character(re.DoGs.adjusted.by.batch$gene)))

    grep("-", unique(as.character(re.DoGs.adjusted.by.batch$gene)))

    emp1.WT <- read.table("/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4Feature/2016-02-10-emp1_WT.gene.downstream.count.rv.m.FC.txt",
        header = T)

    emp1.WT.feature <- emp1.WT[, c(1, 7)]

    colnames(emp1.WT.feature) = c("Row.names", "feature.Count")

    emp1.WT.inhouse <- Re.unadjusted.adjusted$DE[, c(1, 13)]

    Re.unadjusted.adjusted$DE[which(Re.unadjusted.adjusted$DE[, 1] %in% c("uc001aai.1")),
        ]

    emp1.inhouse.feature <- merge(emp1.WT.feature, emp1.WT.inhouse)

    colnames(emp1.inhouse.feature) <- c("Transcript", "Counting_based_featureCounts",
        "Counting_based_bbc_count")

    plot(emp1.inhouse.feature[, 3], emp1.inhouse.feature[, 2], ylab = "featureCounts_based",
        xlab = "bbc_inhouse_developed_program_based")

    head(emp1.inhouse.feature)
}
#' DEAnalysisAdjustByBatch
#'
#' @param df.NT
#'
#' @return
#' @export
#'
#' @examples
#'
#' df.NT<-re.DESeq.DoGs.plus.minus.gene.interested[,c(1,9:16)]
#'
#' df.NT<-Re.unadjusted.adjusted$DE[,c(1,9:16)]
#'
#' rownames(df.NT)<-df.NT[,1]
#' df.NT<- df.NT[,-1]
#'
#' re.DoGs.adjusted.by.batch<-DEAnalysisAdjustByBatch(df.NT.2)
#'

DEAnalysisAdjustByBatch <- function(df.NT, Nbatch)
{


    # print(pkg.env$sample)

    cell <- factor(rep(c("emp", "hela"), Nbatch))
    cell = rep(cell, 2)

    colData <- data.frame(condition = factor(rep(c("Dox", "WT"), c(5, 5))),
        cell = cell)

    dds <- DESeqDataSetFromMatrix(df.NT, colData, formula(~condition))
    dat <- counts(dds)

    mod <- model.matrix(~condition, colData(dds))
    mod0 <- model.matrix(~1, colData(dds))
    svseq.NT <- svaseq(dat, mod, mod0)

    # How many significant surrogate variables is?

    # ncol <- dim(svseq.NT$sv)[2]

    ncol <- svseq.NT$n.sv

    cat(ncol, "\n")

    # ncol<-2 Adjust by surrogate variables 1 and 2
    ddssva <- dds

    if (ncol == 1)
    {
        colData(ddssva) <- cbind2(colData(ddssva), svseq.NT$sv)
        colnames(colData(ddssva))[dim(colData(ddssva))[2]] = paste0("SV", 1)
    } else
    {
        for (i in 1:ncol)
        {
            colData(ddssva) <- cbind2(colData(ddssva), svseq.NT$sv[, i])
            colnames(colData(ddssva))[dim(colData(ddssva))[2]] = paste0("SV",
                i)
        }
    }


    # ddssva$SV1 <- svseq.NT$sv[,1] ddssva$SV2 <- svseq.NT$sv[,2]

    # design(ddssva) <- ~ SV1 + SV2 + condition

    cat("ddssva")
    print(head(ddssva))
    cat("ddsva_done")

    design(ddssva) <- formula(paste("~ ", paste(paste0("SV", 1:ncol), collapse = "+"),
        "+", "condition"))

    # formula(paste('~ ',paste(paste0('SV', 1:2), collapse = '+'),'+',
    # 'condition' ))

    ddssva <- DESeq(ddssva)
    re.DESeq <- results(ddssva)

    re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
    colnames(re.FC)[7] = "FoldChange(WT/Dox)"

    txs.gene <- ReformatTxsGene()

    re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

    re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

    return(re.FC.sorted)

}
#' FilterTranscripts
#'
#' @return
#' @export
#'
#' @examples
#'
#' Re.txs<-FilterTranscripts('/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/','final_list.csv')
#'
#'
#'
FilterTranscripts <- function(input.dir, input.file.pattern)
{

    file.name = paste0(input.dir, dir(input.dir, recursive = TRUE, pattern = input.file.pattern))

    # print(file.name)

    genes.interested <- read.csv(file.name, header = F)

    # print(genes.interested)

    # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

    library("Homo.sapiens")

    txs <- as.data.frame(transcripts(Homo.sapiens, columns = c("TXNAME", "GENEID",
        "SYMBOL")))

    txs.4.genes.interested <- txs[which(txs$SYMBOL %in% genes.interested[, 1]),
        ]

    cat(dim(txs.4.genes.interested))

    # print(head(txs.4.genes.interested))

    # columns(TxDb) library(TxDb.Hsapiens.UCSC.hg19.knownGene) txdb <-
    # TxDb.Hsapiens.UCSC.hg19.knownGene txs1<-transcripts(txdb) columns(txdb)
    # keys <- c('129787') keys <- keys(hgu95av2.db, 'ENTREZID') select(txdb,
    # keys = keys, columns=c('TXCHROM','TXSTART','TXEND','TXSTRAND','TXNAME'),
    # keytype='GENEID') txs2<-select(txdb,keys =
    # keys,columns=c('GENEID','TXCHROM','TXSTART','TXEND','TXSTRAND','TXNAME'),keytype='GENEID')
    # dim(txs2) txs2 txdb
    GRList <- transcriptsBy(txdb, by = "gene")
    # length(GRList) which(names(GRList) %in% c('129787')) GRList[[4237]]
    # as.data.frame(GRList[[4237]])
    num.of.txs.4.gene <- unlist2(lapply(GRList, length))

    boxplot(num.of.txs.4.gene)

    hist(num.of.txs.4.gene)

    # which.max(num.of.txs.4.gene)
    # index<-which(names(unlist(num.of.txs.4.gene))=='4237')
    # GRList.data.frame<-lapply(GRList,as.data.frame) GRList.data.frame[[5922]]

    data2 <- txs
    data3 <- cbind(data2, unlist(data2$SYMBOL))
    colnames(data3)[9] = "Gene"

    txs.num <- count(data3, "Gene")
    hist(txs.num)

    txs.min <- do.call(rbind, by(data3, data3$Gene, function(x) x[which.min(x$width),
        ]))
    txs.max <- do.call(rbind, by(data3, data3$Gene, function(x) x[which.max(x$width),
        ]))

    txs.min.start <- do.call(rbind, by(data3, data3$Gene, function(x) x[which.min(x$start),
        ]))
    txs.max.end <- do.call(rbind, by(data3, data3$Gene, function(x) x[which.max(x$end),
        ]))

    temp <- merge(txs.min.start, txs.max.end, by = "Gene", suffixes = c(".min.start",
        ".max.end"))
    temp2 <- merge(temp, txs.min, by = "Gene")
    temp3 <- merge(temp2, txs.max, by = "Gene", suffixes = c(".shortest", ".longest"))
    temp4 <- merge(temp3, txs.num, by = "Gene")

    # temp3<-temp[order(temp3$Gene),]


    txs.num.4.genes.interested <- txs.num[which(txs.num$Gene %in% genes.interested[,
        1]), ]
    txs.min.4.genes.interested <- txs.min[which(txs.min$Gene %in% genes.interested[,
        1]), ]
    txs.max.4.genes.interested <- txs.max[which(txs.max$Gene %in% genes.interested[,
        1]), ]

    re <- list(txs = txs, txs.num = txs.num, txs.min = txs.min, txs.max = txs.max,
        txs.min.start.max.end = temp4, txs.4.genes.interested = txs.4.genes.interested,
        txs.num.4.genes.interested = txs.num.4.genes.interested, txs.min.4.genes.interested = txs.min.4.genes.interested,
        txs.max.4.genes.interested = txs.max.4.genes.interested)

    return(re)

}

# index<-which(Re.txs$txs.min.start.max.end$seqnames.min.start!=Re.txs$txs.min.start.max.end$seqnames.max.end)
# dim(Re.txs$txs.min.start.max.end[index,])
# names(Re.txs$txs.min.start.max.end)
# head(Re.txs$txs.min.start.max.end[index,c(1,2,10)],275)

# data2<-head(Re.txs$txs.all,100) data2

# do.call(rbind, by(data2, unlist(data2$SYMBOL), function(x)
# x[which.min(x$width),]))

# do.call(rbind,tapply(1:nrow(data2), unlist(data2$SYMBOL), function(x)
# data2[x,][which.min(data2$width[x]),])) dara



cmd2 = "wc -l /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/*.rm.exon.intron.hg19.bed"

Re.out <- system(cmd2, intern = TRUE, ignore.stderr = TRUE)

re <- strsplit(Re.out, split = "\\/")

re2 <- cbind(trimws(do.call("rbind", lapply(re[1:4], "[[", 1))), do.call("rbind",
    lapply(re[1:4], "[[", 9)))
#' InstallRequiredPackage
#'
#' @return
#' @export
#'
#' @examples
#' InstallRequiredPackage()
#'
InstallRequiredPackage <- function()
{

    source("https://bioconductor.org/biocLite.R")
    biocLite("Homo.sapiens")
    biocLite("hgu95av2.db")
    biocLite("rhdf5")
    install.packages("plyr")
    biocLite("rtracklayer")
    install.packages("devtools")
    devtools::install_github("pachterlab/sleuth")
    biocLite("biomaRt")
    install.packages("corrplot")
    biocLite("sva")
    source("https://bioconductor.org/biocLite.R")
    biocLite("ShortRead")
    install.packages("doMC")

    library(Homo.sapiens)
    library(plyr)
    library(rtracklayer)
    library(annotate)

}
#' Loadlibrary
#'
#' @return
#' @export
#'
#' @examples
#'
#' Loadlibrary()
#'
Loadlibrary <- function()
{
    suppressPackageStartupMessages(library(Homo.sapiens))
    suppressMessages(library(hgu95av2.db))
    suppressPackageStartupMessages(library(plyr))
    suppressPackageStartupMessages(library(rtracklayer))
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    suppressPackageStartupMessages(library(annotate))
    suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
    suppressPackageStartupMessages(library("sleuth"))
    suppressPackageStartupMessages(library(corrplot))
    suppressPackageStartupMessages(library(sva))
    suppressPackageStartupMessages(library(DESeq2))
}
#' MatchTxs2Gene
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' dir.name='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/'
#' output.file='InterestedGeneTxsSorted.csv'
#'
#'
#' Re.interested.gene<-MatchTxs2Gene.R(re.rMAT,Re.txs,dir.name,output.file)
#'
#'
#' Re.interested.gene[which(Re.interested.gene$SYMBOL %in% c('XYLB')),]
#'
#'
#'
MatchTxs2Gene.R <- function(re.rMAT, Re.txs, dir.name, out.file)
{
    txs.DE <- re.rMAT$DE
    txs.all <- Re.txs$txs
    head(txs.DE)
    head(txs.all)
    rownames(txs.all) <- txs.all$TXNAME
    txs.DE.all <- merge(txs.DE, txs.all, by = 0)

    txs.4.certain.gene <- txs.DE.all[which(txs.DE.all$SYMBOL %in% Re.txs$txs.4.genes.interested$SYMBOL),
        ]

    gene <- unlist2(txs.4.certain.gene$SYMBOL)

    txs.4.certain.gene.2 <- cbind(txs.4.certain.gene, gene)

    txs.4.certain.gene.3 <- txs.4.certain.gene.2[order(txs.4.certain.gene.2$gene),
        ]

    write.csv(txs.4.certain.gene.3, file = paste0(dir.name, out.file))

    return(txs.4.certain.gene.2)

}
#' Title
#'
#' @param re.out
#' @param dir.name
#' @param by_what
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' dir.name=/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/'
#'
#' by_what='*.rm.exon.intron.hg19.bed'
#'
#' NormalizedCount(re.out,dir.name,by_what)

NormalizedCount <- function(re.out, dir.name, by_what)
{


    pattern_check = grep(glob2rx("*.rm.exon.intron.hg19.bed"), by_what)

    if (length(pattern_check) != 0)
    {
        cmd2 = paste0("wc -l ", dir.name, by_what)
    } else
    {
        cmd2 = paste0("wc -l $(ls ", dir.name, by_what, " | grep -v hg19)")
    }

    num.intergenic.reads <- system(cmd2, intern = TRUE, ignore.stderr = TRUE)

    num.intergenic.reads.2 <- strsplit(num.intergenic.reads, split = "\\/")

    # cat('OK') print(num.intergenic.reads.2) cat('OK')

    n = length(num.intergenic.reads.2) - 1
    num.intergenic.reads.3 <- cbind(trimws(do.call("rbind", lapply(num.intergenic.reads.2[1:n],
        "[[", 1))), do.call("rbind", lapply(num.intergenic.reads.2[1:n], "[[",
        9)))

    name.sample <- unlist(lapply(strsplit(num.intergenic.reads.3[, 2], split = "\\."),
        "[[", 1))

    print(names(re.out))

    print(name.sample)

    normalized.count.gene.matched.final <- list()

    for (i in 1:8)
    {
        count.gene.matched <- re.out[(grep(name.sample[i], names(re.out)))]
        num.intergenic.reads.matched <- as.numeric(num.intergenic.reads.3[grep(name.sample[i],
            num.intergenic.reads.3[, 2]), 1])

        # print(count.gene.matched) print(num.intergenic.reads.matched)

        normalized.count.gene.matched <- (count.gene.matched[[1]][, 1]/num.intergenic.reads.matched) *
            (10^6)
        count.gene.matched.2 <- data.frame(cbind(as.character(count.gene.matched[[1]][,
            2]), count.gene.matched[[1]][, 1], normalized.count.gene.matched))

        colnames(count.gene.matched.2) = c("Gene", paste0("Counts.", name.sample[i]),
            paste0("Normalized.Counts.", name.sample[i]))

        normalized.count.gene.matched.final[[i]] <- count.gene.matched.2
    }

    merge.all <- function(x, y)
    {
        merge(x, y, all = TRUE, by = "Gene")
    }

    output <- Reduce(merge.all, normalized.count.gene.matched.final)

    row.has.na <- apply(output, 1, function(x)
    {
        any(is.na(x))
    })
    print(sum(row.has.na))

    final.filtered <- output[!row.has.na, ]

    re <- list(final.filtered = final.filtered, num.intergenic.reads.2 = num.intergenic.reads.2)

    return(re)

}
#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' raw.count<-re.rMAT.com.total$DE[,10:17]
#'
#' raw.count.dds<-counts(re.rMAT.com.total.test.2$dds)
#'
#' raw.count.DE<-re.rMAT.com.total.test.2$DE[,c(1,10:17)]
#'
#' raw.count.DE.2<-raw.count.DE[match(rownames(raw.count.dds),raw.count.DE$tx_name),]
#'
#' head(raw.count.DE.2)
#'
#' raw.count.2<-apply(raw.count,2,as.numeric)
#'
#' head(raw.count)
#'
#'
#' nor.fac<-norm_factors(raw.count.2)
#'
#' estimateSizeFactorsForMatrix(raw.count.2)
#'
#'
norm_factors <- function(mat)
{
    # nz <- apply(mat, 1, function(row) !any(round(row) == 0)) mat_nz <-
    # mat[nz,] p <- ncol(mat) geo_means <- exp(apply(mat_nz, 1, function(row)
    # mean(log(row)) ))

    loggeomeans <- rowMeans(log(mat))

    counts <- mat

    sf <- apply(counts, 2, function(cnts)
    {
        exp(median(log(cnts) - loggeomeans))
    })

    # <- sweep(mat_nz, 1, geo_means, `/`)

    # sf <- apply(s, 2, median) scaling <- exp( (-1 / p) * sum(log(sf)))

    # sf * scaling
    sf
}
#' @title ProcessOutputFilesFrom3UTR
#'
#' @description  read the mapping results from the Dogs of each gene
#'
#' @param dir.name: the path for input files
#' @param input.file.pattern: input file pattern
#'
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/'
#'
#' #new data
#' dir.name='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/'
#' input.file.pattern='*downstream.count.hg19.strand.based.txt'
#'
#' com.input.file.pattern='*downstream.com.count.hg19.strand.based.txt'
#'
#' total.com.input.file.pattern='*.gene.downstream.count.hg19.strand.based.total2.txt'
#'
#' #use intergenic reads as normal factor
#' normal.factor='*.rm.exon.intron.hg19.bed'
#'
#' #use total reads as normal factor
#' normal.factor='*.bed'
#'
#' sink('testcode.txt')
#' re.rMAT<-ProcessOutputFilesFrom3UTR(dir.name,input.file.pattern,normal.factor,'intergenic')
#' re.rMAT.com<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,'com_intergenic')
#'
#' re.rMAT.com.total<-ProcessOutputFilesFrom3UTR(dir.name,total.com.input.file.pattern,normal.factor,'com_total')
#'
#' re.rMAT.com.total.interested<-re.rMAT.com.total$DE[which(re.rMAT.com.total$DE$gene %in% genes.interested[,1]),]
#'
#' re.com.total.pos<-re.rMAT.com.total.interested[which(re.rMAT.com.total.interested$`FoldChange(WT/Dox)`<1),]
#' re.com.total.neg<-re.rMAT.com.total.interested[which(re.rMAT.com.total.interested$`FoldChange(WT/Dox)`>=1),]
#'
#'
#' re.com.total.pos[which(re.com.total.pos$gene %in% c('TPCN2')),]
#'
#' re.rMAT.com.interested<-re.rMAT.com$DE[which(re.rMAT.com$DE$gene %in% genes.interested[,1]),]
#' re.pos<-re.rMAT.com.interested[which(re.rMAT.com.interested$`FoldChange(WT/Dox)`<1),]
#' re.neg<-re.rMAT.com.interested[which(re.rMAT.com.interested$`FoldChange(WT/Dox)`>=1),]
#'
#'
#' re.rMAT.total<-ProcessOutputFilesFrom3UTR(dir.name,input.file.pattern,normal.factor,'total')
#' re.rMAT.com.total<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,'com_total')
#'
#' sink()
#'
#' sink('test_com_total.txt')
#' re.rMAT.com.total.test<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,'com_total_test')
#' sink()
#'
#' re.rMAT.com.total.test.2<-ProcessOutputFilesFrom3UTR(dir.name,total.com.input.file.pattern,normal.factor,'com_total_test_3')
#'
#' geoMeans <- exp(rowMeans(log(counts(re.rMAT.com.total.test$dds))))
#'
#' temp.dds<-re.rMAT.com.total.test.2$dds
#'
#' normalizationFactors(temp.dds)

#' size..factor<-estimateSizeFactors(re.rMAT.com.total.test$dds,geoMeans=geoMeans)
#' sizeFactors(size..factor)
#'
#' temp.dds <- estimateSizeFactors(temp.dds)
#'
#' count.norm<-counts(temp.dds,normalized = TRUE)
#'
#' head(count.norm)
#'
#' DE.norm.with.rpkm.norm<-cbind(count.norm,re.rMAT.com.total.test.2$MergedNorma[,c(1,3,7,11,13,5,9,15,17)])
#'
#' DE.norm.with.rpkm.norm.2<-DE.norm.with.rpkm.norm[,-9]
#'
#' col.name<-c(paste0('D.',sapply(strsplit(colnames(DE.norm.with.rpkm.norm.2)[1:8],'\\.'),'[[',2)),
#' paste0('R.',sapply(strsplit(colnames(DE.norm.with.rpkm.norm.2)[9:16],'\\.'),'[[',3)))
#'
#' col.name.2<-gsub('2016-02-10-','',col.name)
#'
#' cbind(colnames(DE.norm.with.rpkm.norm.2),col.name.2)
#'
#' colnames(DE.norm.with.rpkm.norm.2)<-col.name.2
#'
#'
#' DE.norm.with.rpkm.norm.3<-apply(DE.norm.with.rpkm.norm.2,2,as.numeric)
#'
#'par(mfrow = c(4, 2))  # 3 rows and 2 columns
#'for (i in 1:8) {
#' plot(DE.norm.with.rpkm.norm.3[,c(i,i+8)])
#' model <- lm(DE.norm.with.rpkm.norm.3[,i+8] ~ DE.norm.with.rpkm.norm.3[,i], data = as.data.frame(DE.norm.with.rpkm.norm.3))
#' abline(model, col = 'red')
#'}
#'
#' plot(apply(DE.norm.with.rpkm.norm.2,2,as.numeric)[,1:3])
#'
#' M <- cor(DE.norm.with.rpkm.norm.3[,1:8])
#'
#' corrplot.mixed(M)
#'
#' corrplot(M,method='ellipse',order = 'hclust', addrect = 2)
#'
#' corrplot(M,method='number',order = 'hclust', addrect = 2)
#'
#' corrplot(M,method='number',order = 'FPC', addrect = 2)
#'
#
#' corrplot(M,method='ellipse')
#'
#' corrplot(M,method='ellipse',type='upper',)
#' corrplot(M, method='number')
#'
#' corrplot(M,method='ellipse',order='hclust', addrect=4)

#'
#' size.factor.2<-estimateSizeFactors(re.rMAT.com.total.test.2$dds)
#' sizeFactors(size.factor.2)

# corrplot(M, order='hclust', addrect=2, col='wb', bg='gold2')


ProcessOutputFilesFrom3UTR <- function(dir.name, input.file.pattern, normal.factor,
    out)
    {

    file.name = paste0(dir.name, dir(dir.name, recursive = TRUE, pattern = input.file.pattern))
    file.name.2 <- as.list(file.name)

    names(file.name.2) = sapply(strsplit(file.name, split = "\\/"), "[[", 9)

    print(file.name.2)

    re.out <- lapply(file.name.2, function(u)
    {
        re = read.table(u, header = F)
        colnames(re) = c("Count", "GeneName")
        re
    })

    temp.name <- strsplit(names(file.name.2), split = "\\.")

    temp.name.2 <- trimws(do.call("rbind", lapply(temp.name, "[[", 1)))

    final.filtered.norm <- NormalizedCount(re.out, dir.name, normal.factor)

    final.filtered <- final.filtered.norm$final.filtered
    num.intergenic.reads.2 <- final.filtered.norm$num.intergenic.reads.2

    n = length(num.intergenic.reads.2) - 1
    # suppressPackageStartupMessages(library(DESeq2))

    print(head(final.filtered))

    countData <- apply(final.filtered[, c(2, 6, 10, 12, 4, 8, 14, 16)], 2, as.numeric)

    rownames(countData) <- final.filtered[, 1]

    colData <- data.frame(condition = factor(c(rep("Dox", 4), rep("WT", 4))))
    dds <- DESeqDataSetFromMatrix(countData, colData, formula(~condition))

    # size.factor<-estimateSizeFactors(dds)

    # print(size.factor)

    re.DESeq <- results(DESeq(dds))

    print(re.DESeq)

    re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
    colnames(re.FC)[7] = "FoldChange(WT/Dox)"

    txs.gene <- ReformatTxsGene()

    re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

    # print(head(final.filtered))

    colnames(final.filtered)[1] <- "tx_name"

    final.filtered.2 <- cbind(final.filtered[, grep(glob2rx("tx_name"), colnames(final.filtered))],
        final.filtered[, grep(glob2rx("Normalized*emp*Dox*"), colnames(final.filtered))],
        final.filtered[, grep(glob2rx("Normalized*hela*dox*"), colnames(final.filtered))],
        final.filtered[, grep(glob2rx("Normalized*emp*WT*"), colnames(final.filtered))],
        final.filtered[, grep(glob2rx("Normalized*hela*wt*"), colnames(final.filtered))])

    colnames(final.filtered.2)[1] <- "tx_name"

    re.FC <- merge(re.FC, final.filtered.2, by = "tx_name")

    re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

    input.file.pattern.2 <- sub("\\*", "_", input.file.pattern)

    write.csv(re.FC.sorted, file = paste0(dir.name, "3UTR_DE_", out, input.file.pattern.2,
        ".csv"))

    re2 <- cbind(trimws(do.call("rbind", lapply(num.intergenic.reads.2[1:n],
        "[[", 1))), do.call("rbind", lapply(num.intergenic.reads.2[1:n], "[[",
        9)))

    re.out.3 <- list(ReadDownstream45kb = re.out, intergenic_reads = re2, MergedNorma = final.filtered,
        DE = re.FC.sorted, dds = dds)

    return(re.out.3)

}
#' @title ProcessOutputFilesFrom3UTR_recount
#'
#' @description read the mapping results from the Dogs of each gene
#'
#' @param dir.name: the path for input files
#' @param input.file.pattern: input file pattern
#'
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#' #new data
#' dir.name='/media/H_driver/2016/Ramin_azhang/'
#'
#' input.file.pattern='*downstream.count.plus.minus.and.minus.plus.3.FPKM.xls'
#'
#' input.file.pattern.4.interested.gene='final_list.csv'
#'
#' out_put='for_bioinfo_core/RNA_seq/Results4Check/'
#'
#' out_by_p=paste0(out_put,'DE_by_p.csv')
#' out_by_gene=paste0(out_put,'DE_by_gene.csv')
#' out_gene_interested=paste0(out_put,'DE_by_interested_gene.csv')
#' out_gene_interested_sig=paste0(out_put,'DE_by_interested_gene_sig.csv')
#'
#' re<-ProcessOutputFilesFrom3UTR_recount(dir.name,input.file.pattern,input.file.pattern.4.interested.gene,out_by_p,out_by_gene,out_gene_interested,out_gene_interested_sig)
#'

ProcessOutputFilesFrom3UTR_recount <- function(dir.name, input.file.pattern,
    input.file.pattern.4.interested.gene, out_by_p, out_by_gene, out_gene_interested,
    out_gene_interested_sig)
    {

    file.name = paste0(dir.name, dir(dir.name, recursive = TRUE, pattern = input.file.pattern))
    file.name.2 <- as.list(file.name)

    names(file.name.2) = sapply(strsplit(file.name, split = "\\/"), "[[", 9)

    # print(file.name.2)

    re.out <- lapply(file.name.2, function(u)
    {
        re = read.table(u, header = F)
        colnames(re) = c("chrom", "st", "end", "accession", "mRNA_size", "gene_strand",
            "Frag_count", "FPM", "FPKM")
        re
    })

    temp.name <- strsplit(names(file.name.2), split = "\\.")
    temp.name.2 <- trimws(do.call("rbind", lapply(temp.name, "[[", 1)))

    re.out.2 <- do.call("cbind", lapply(re.out, function(x)
    {
        x
    }))

    re.out.2.name <- gsub("gene.downstream.count.plus.minus.and.minus.plus.3.FPKM.xls.",
        "", colnames(re.out.2))

    colnames(re.out.2) <- re.out.2.name

    chrom <- as.data.frame(apply(re.out.2[, grep("chrom", colnames(re.out.2))],
        1, unique))
    st <- as.data.frame(apply(re.out.2[, grep("\\.st", colnames(re.out.2))],
        1, unique))
    end <- as.data.frame(apply(re.out.2[, grep("end", colnames(re.out.2))],
        1, unique))
    accession <- as.data.frame(apply(re.out.2[, grep("accession", colnames(re.out.2))],
        1, unique))
    mRNA_size <- as.data.frame(apply(re.out.2[, grep("mRNA_size", colnames(re.out.2))],
        1, unique))
    gene_strand <- as.data.frame(apply(re.out.2[, grep("gene_strand", colnames(re.out.2))],
        1, unique))

    colnames(chrom) <- "chrom"
    colnames(st) <- "st"
    colnames(end) <- "end"
    colnames(accession) <- "accession"
    colnames(mRNA_size) <- "mRNA_size"
    colnames(gene_strand) <- "gene_strand"

    temp <- cbind.data.frame(chrom, st, end, accession, mRNA_size, gene_strand)

    Frag_count <- re.out.2[, grep("Frag_count", colnames(re.out.2))]
    FPM <- re.out.2[, grep("FPM", colnames(re.out.2))]
    FPKM <- re.out.2[, grep("FPKM", colnames(re.out.2))]

    re.out.3 <- cbind(temp, Frag_count, FPM, FPKM)

    countData <- apply(re.out.3[, c(7, 9, 11, 12, 8, 10, 13, 14)], 2, as.numeric)

    rownames(countData) <- re.out.3[, 4]

    dat <- countData
    idx <- rowMeans(dat) > 1
    dat <- dat[idx, ]


    cell <- factor(rep(c("emp", "hela"), c(2, 2)))
    cell = rep(cell, 2)

    colData <- data.frame(condition = factor(rep(c("Dox", "WT"), c(4, 4))),
        cell = cell)

    dds <- DESeqDataSetFromMatrix(dat, colData, formula(~condition))

    # counts(dds)


    mod <- model.matrix(~condition, colData(dds))
    mod0 <- model.matrix(~1, colData(dds))
    svseq.NT <- svaseq(counts(dds), mod, mod0)

    par(mfrow = c(2, 1), mar = c(3, 5, 3, 1))
    stripchart(svseq.NT$sv[, 1] ~ dds$cell, vertical = TRUE, main = "SV1")
    abline(h = 0)
    stripchart(svseq.NT$sv[, 2] ~ dds$cell, vertical = TRUE, main = "SV2")
    abline(h = 0)


    # Adjust by surrogate variables 1 and 2
    ddssva <- dds
    ddssva$SV1 <- svseq.NT$sv[, 1]
    ddssva$SV2 <- svseq.NT$sv[, 2]
    design(ddssva) <- ~SV1 + SV2 + condition
    ddssva <- DESeq(ddssva)

    # print(sizeFactors(estimateSizeFactors(ddssva)))
    # print(head(counts(ddssva,normalized = TRUE)))

    res.ddssva.2.sva.adjusted <- results(ddssva)
    # head(res.ddssva.2.sva.adjusted) summary(res.ddssva.2.sva.adjusted)

    re.FC <- cbind(as.data.frame(res.ddssva.2.sva.adjusted), 2^res.ddssva.2.sva.adjusted[,
        2])
    colnames(re.FC)[7] = "FoldChange(WT/Dox)"

    MatchGene <- function(res.sva.NT, re.rMAT.com.total.test.2)
    {
        temp <- re.rMAT.com.total.test.2$DE[, c(1, 23, 24)]
        rownames(temp) <- temp[, 1]
        # temp<-temp[,-1]
        temp2 <- merge(as.data.frame(res.sva.NT), temp, by = 0, sort = FALSE)
        re.FC <- temp2
        re.FC.sorted <- re.FC[order(re.FC$pvalue), ]
        return(re.FC.sorted)
    }

    re.FC.sorted <- MatchGene(re.FC, re.rMAT.com.total.test.2)

    rownames(re.out.3) = re.out.3$accession
    re.FC.sorted <- re.FC.sorted[, -1]

    rownames(re.FC.sorted) = re.FC.sorted$tx_name

    re.out.3.temp <- re.out.3[, c(1:6, 7, 9, 11, 12, 8, 10, 13, 14, 15, 17,
        19, 20, 16, 18, 21, 22, 23, 25, 27, 28, 24, 26, 29, 30)]

    re.FC.sorted <- merge(re.FC.sorted, re.out.3.temp, by = 0, sort = FALSE)

    p.one.sided <- unlist(lapply(re.FC.sorted$stat, function(x)
    {
        convert.z.score(x, two.side = F)
    }))

    p.one.sided.adjust <- p.adjust(p.one.sided, method = "BH")

    re.FC.sorted <- cbind(re.FC.sorted, p.one.sided, p.one.sided.adjust)

    re.FC.sorted.sorted.by.gene <- re.FC.sorted[order(re.FC.sorted$gene), ]

    write.csv(re.FC.sorted, file = paste0(dir.name, out_by_p), row.names = FALSE,
        quote = FALSE)

    write.csv(re.FC.sorted.sorted.by.gene, file = paste0(dir.name, out_by_gene),
        row.names = FALSE, quote = FALSE)

    file.name = paste0(dir.name, dir(dir.name, recursive = TRUE, pattern = input.file.pattern.4.interested.gene))

    genes.interested <- read.csv(file.name, header = F)

    re.genes.interested <- re.FC.sorted[which(re.FC.sorted$gene %in% genes.interested[,
        1]), ]

    re.genes.interested.sig <- re.genes.interested[which(re.genes.interested$pvalue <
        0.05 & re.genes.interested$log2FoldChange < 0), ]

    write.csv(re.genes.interested, file = paste0(dir.name, out_gene_interested),
        row.names = FALSE, quote = FALSE)

    write.csv(re.genes.interested.sig, file = paste0(dir.name, out_gene_interested_sig),
        row.names = FALSE, quote = FALSE)

    re <- list(re.FC.sorted = re.FC.sorted, re.FC.sorted.sorted.by.gene = re.FC.sorted.sorted.by.gene,
        re.genes.interested = re.genes.interested, re.genes.interested.sig = re.genes.interested.sig)

    return(re)

}
#' @title ProcessOutputFilesFrom5Cases
#'
#' @description  read the mapping results from the Dogs of each gene
#'
#' @param dir.name: the path for input files
#' @param input.file.pattern: input file pattern
#'
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts5CasesEachSample/'
#' input.file.pattern='*count.2.txt'
#'
#' dir.name.gene.list='/media/H_driver/2016/Ramin_azhang/'
#' pattern.4.gene.list='final_list.csv'
#'

#' out.dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts5CasesEachSample/'
#'
#' out.file.pattern.interested='DoGs_adjust_by_batch_interested_gene'
#' out.file.pattern.positive.gene='DoGs_adjust_by_batch_positive'
#' out.file.pattern.negative.gene='DoGs_adjust_by_batch_negative'
#' out.file.pattern.all= 'DoGs_adjust_by_batch_all'
#'
#' Re.unadjusted.adjusted<-ProcessOutputFilesFrom5Cases(dir.name,input.file.pattern,out.dir.name,out.file.pattern.interested,
#' out.file.pattern.positive.gene,
#' out.file.pattern.negative.gene,
#' out.file.pattern.all,
#' dir.name.gene.list,
#' pattern.4.gene.list,
#' adjust_by_batch='YES')
#'
#' save.image(file=paste0(out.dir.name,'re_save_2.RData'))
#' savehistory(file=paste0(out.dir.name,'re_save_2.Rhistory'))
#'

ProcessOutputFilesFrom5Cases <- function(dir.name, input.file.pattern, out.dir.name,
    out.file.pattern.interested, out.file.pattern.positive.gene, out.file.pattern.negative.gene,
    out.file.pattern.all, dir.name.gene.list, pattern.4.gene.list, adjust_by_batch)
    {

    file.name = paste0(dir.name, dir(dir.name, recursive = TRUE, pattern = input.file.pattern))
    file.name.2 <- as.list(file.name)

    names(file.name.2) = sapply(strsplit(file.name, split = "\\/"), "[[", 7)

    print(file.name.2)

    re.out <- lapply(file.name.2, function(u)
    {
        re = read.table(u, header = F)
        colnames(re) = c("Count", "GeneName")
        re
    })

    temp.name <- strsplit(names(file.name.2), split = "\\.")

    temp.name.2 <- trimws(do.call("rbind", lapply(temp.name, "[[", 1)))

    temp.name.3 <- unique(temp.name.2[, 1])

    temp.name.4 <- as.list(temp.name.3)
    names(temp.name.4) <- temp.name.3

    txs.name <- lapply(re.out, "[[", 2)
    txs.name <- unlist(txs.name)
    txs.name <- unique(as.character(txs.name))
    txs.name.count <- data.frame(Count = rep(0, length(txs.name)), GeneName = txs.name)

    # dim(txs.name.count)

    Process4EachSample <- function(one_sample, txs.name.count, re.out)
    {

        count.gene.matched <- re.out[(grep(one_sample, names(re.out)))]

        count.gene.matched.minus.gene <- count.gene.matched[(grep("minus.gene",
            names(count.gene.matched)))]
        count.gene.matched.plus.gene <- count.gene.matched[(grep("plus.gene",
            names(count.gene.matched)))]

        ReformatCount <- function(count.gene.matched.minus.gene, txs.name.count)
        {

            reformat.count.gene.matched <- lapply(count.gene.matched.minus.gene,
                function(x, txs.name.count)
                {
                  y <- txs.name.count
                  if (dim(x[which(x$GeneName %in% y$GeneName), ])[1] != 0)
                  {
                    y[match(x$GeneName, y$GeneName), ]$Count <- x$Count
                  }
                  y
                }, txs.name.count)


            reformat.count.gene.matched.2 <- cbind(apply(do.call("cbind", lapply(reformat.count.gene.matched,
                "[", 2)), 1, unique), do.call("cbind", lapply(reformat.count.gene.matched,
                "[", 1)))

            colnames(reformat.count.gene.matched.2) = c("GeneName", names(count.gene.matched.minus.gene))

            gene <- as.character(reformat.count.gene.matched.2$GeneName)

            count.DoGs.plus.read <- apply(reformat.count.gene.matched.2[, c(grep("plus.read.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.start.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.end.count.2.txt",
                colnames(reformat.count.gene.matched.2)))], 1, sum)

            count.DoGs.minus.read <- apply(reformat.count.gene.matched.2[, c(grep("minus.read.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.start.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.end.count.2.txt",
                colnames(reformat.count.gene.matched.2)))], 1, sum)

            count.gene.plus.read <- apply(reformat.count.gene.matched.2[, c(grep("plus.read.gene.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.start.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("plus.read.overlap.Gene.and.DoGs.end.count.2.txt",
                colnames(reformat.count.gene.matched.2)))], 1, sum)

            count.gene.minus.read <- apply(reformat.count.gene.matched.2[, c(grep("minus.read.gene.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.start.count.2.txt",
                colnames(reformat.count.gene.matched.2)), grep("minus.read.overlap.Gene.and.DoGs.end.count.2.txt",
                colnames(reformat.count.gene.matched.2)))], 1, sum)


            reformat.count.gene.matched.3 <- as.data.frame(cbind(gene, count.gene.plus.read,
                count.gene.minus.read, count.DoGs.plus.read, count.DoGs.minus.read))

            return(reformat.count.gene.matched.3)
        }

        re1 <- ReformatCount(count.gene.matched.minus.gene, txs.name.count)
        re2 <- ReformatCount(count.gene.matched.plus.gene, txs.name.count)
        re <- rbind(re1, re2)

        re <- re[order(re$gene), ]

        return(re)
    }

    re.8.samples <- lapply(temp.name.4, function(u, txs.name.count, re.out)
    {
        YY <- Process4EachSample(u, txs.name.count, re.out)
        YY
    }, txs.name.count, re.out)

    DEAnalysis <- function(countData)
    {
        colData <- data.frame(condition = factor(c(rep("Dox", 4), rep("WT",
            4))))

        dds <- DESeqDataSetFromMatrix(countData, colData, formula(~condition))

        # size.factor<-estimateSizeFactors(dds)

        # print(size.factor)

        re.DESeq <- results(DESeq(dds))

        re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
        colnames(re.FC)[7] = "FoldChange(WT/Dox)"

        txs.gene <- ReformatTxsGene()

        re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

        re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

        return(re.FC.sorted)

    }

    ProcessEachCorner <- function(re.8.samples, i)
    {

        gene <- apply(do.call("cbind", lapply(re.8.samples, "[", 1)), 1, unique)
        count.gene.plus.read.8.samples <- cbind(gene, do.call("cbind", lapply(re.8.samples,
            "[", i)))

        colnames(count.gene.plus.read.8.samples) = c("gene", names(re.8.samples))

        df <- count.gene.plus.read.8.samples
        df <- apply(df[, -1], 2, as.numeric)
        index <- rowSums(df[, -1]) > 0
        dff <- count.gene.plus.read.8.samples[index, ]
        rownames(dff) <- dff$gene
        dff <- dff[, -1]

        countData <- apply(dff[, c(1, 3, 5, 6, 2, 4, 7, 8)], 2, as.numeric)

        rownames(countData) <- rownames(dff)

        re.FC <- DEAnalysis(countData)

        # re.DoGs.adjusted.by.batch<-DEAnalysisAdjustByBatch(countData)

        # re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
        # colnames(re.FC)[7]='FoldChange(WT/Dox)' txs.gene<-ReformatTxsGene()
        # re.FC<-list(re.FC=re.FC,re.DoGs.adjusted.by.batch=re.DoGs.adjusted.by.batch)

        return(re.FC)

    }

    re.BL <- ProcessEachCorner(re.8.samples, 2)  #BL
    re.TL <- ProcessEachCorner(re.8.samples, 3)  #TL
    re.BR <- ProcessEachCorner(re.8.samples, 4)  #BR
    re.TR <- ProcessEachCorner(re.8.samples, 5)  #TR

    # head(re.BL) head(re.BR)

    # head(re.TL) head(re.TR)

    # re.TR[which(re.TR$Row.names =='uc003tzi.4'),]

    # Check for uc003tzi.4:gene(+) re.BL[[1]][which(re.BL[[1]]$Row.names
    # =='uc003tzi.4'),] #BL for gene(+) re.TL[[1]][which(re.TL[[1]]$Row.names
    # =='uc003tzi.4'),] #TL for gene(+) gene
    # re.BR[[1]][which(re.BR[[1]]$Row.names =='uc003tzi.4'),] #BR for gene(+)
    # re.TR[[1]][which(re.TR[[1]]$Row.names =='uc003tzi.4'),] #TR for gene(+)
    # DoGs

    # Check for uc002szf.1:gene(-) re.BL[[1]][which(re.BL[[1]]$Row.names
    # =='uc002szf.1'),] #BR for gene(-) re.TL[[1]][which(re.TL[[1]]$Row.names
    # =='uc002szf.1'),] # re.BR[[1]][which(re.BR[[1]]$Row.names
    # =='uc002szf.1'),] # re.TR[[1]][which(re.TR[[1]]$Row.names
    # =='uc002szf.1'),] #

    # Check for uc001aac.4:gene(-) re.BL[[1]][which(re.BL[[1]]$Row.names
    # =='uc001aac.4'),]#BR for gene(-) Gene
    # re.TL[[1]][which(re.TL[[1]]$Row.names =='uc001aac.4'),]#TR for gene(-)
    # re.TR[[1]][which(re.BR[[1]]$Row.names =='uc001aac.4'),]#TL for gene(-)
    # re.BR[[1]][which(re.TR[[1]]$Row.names =='uc001aac.4'),]#BL for gene(-)
    # DoGs

    # if(adjust_by_batch=='NO'){
    re.BL.4.plus.gene.BR.4.minus.gene <- re.BL
    re.TL.4.plus.gene.TR.4.minus.gene <- re.TL
    re.BR.4.plus.gene.BL.4.minus.gene <- re.BR
    re.TR.4.plus.gene.TL.4.minus.gene <- re.TR
    # }else{ re.BL.4.plus.gene.BR.4.minus.gene<-re.BL[[2]]
    # re.TL.4.plus.gene.TR.4.minus.gene<-re.TL[[2]]
    # re.BR.4.plus.gene.BL.4.minus.gene<-re.BR[[2]]
    # re.TR.4.plus.gene.TL.4.minus.gene<-re.TR[[2]] }

    # Get the counts for DoGs of plus and minus gene
    DoGs.4.plus.Gene <- re.TR.4.plus.gene.TL.4.minus.gene[which(re.TR.4.plus.gene.TL.4.minus.gene$strand ==
        "+"), ]
    DoGs.4.minus.Gene <- re.BR.4.plus.gene.BL.4.minus.gene[which(re.BR.4.plus.gene.BL.4.minus.gene$strand ==
        "-"), ]

    DoGs.4.plus.minus.Gene <- rbind(DoGs.4.plus.Gene, DoGs.4.minus.Gene)

    # DoGs.4.plus.Gene[which(DoGs.4.plus.Gene$Row.names =='uc003tzi.4'),] #TR
    # for gene(+) DoGs DoGs.4.minus.Gene[which( DoGs.4.minus.Gene$Row.names
    # =='uc001aac.4'),]#BL for gene(-) DoGs

    GeneTypeBasedDE <- function(DoGs.4.plus.Gene)
    {
        Count.DoGs.4.plus.Gene <- DoGs.4.plus.Gene[, c(1, 9:16)]
        rownames(Count.DoGs.4.plus.Gene) <- Count.DoGs.4.plus.Gene$Row.names
        Count.DoGs.4.plus.Gene.2 <- Count.DoGs.4.plus.Gene[, -1]
        head(Count.DoGs.4.plus.Gene.2)

        if (adjust_by_batch == "NO")
        {
            re.DESeq.DoGs.plus.gene <- DEAnalysis(Count.DoGs.4.plus.Gene.2)
        } else
        {
            re.DESeq.DoGs.plus.gene <- DEAnalysisAdjustByBatch(Count.DoGs.4.plus.Gene.2)
        }

        return(re.DESeq.DoGs.plus.gene)
    }

    re.DESeq.DoGs.plus.gene <- GeneTypeBasedDE(DoGs.4.plus.Gene)
    # head(re.DESeq.DoGs.plus.gene[which(re.DESeq.DoGs.plus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.gene$pvalue<0.05),])

    re.DESeq.DoGs.minus.gene <- GeneTypeBasedDE(DoGs.4.minus.Gene)
    # head(re.DESeq.DoGs.minus.gene[which(re.DESeq.DoGs.minus.gene$log2FoldChange<0&re.DESeq.DoGs.minus.gene$pvalue<0.05),])

    re.DESeq.DoGs.plus.minus.gene <- GeneTypeBasedDE(DoGs.4.plus.minus.Gene)
    # head(re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.minus.gene$pvalue<0.05),])

    gene.interested <- ReadGeneList(dir.name.gene.list, pattern.4.gene.list)

    re.DESeq.DoGs.plus.minus.gene.interested <- re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$gene %in%
        gene.interested[, 1]), ]

    # re.DESeq.DoGs.plus.minus.gene.interested[which(re.DESeq.DoGs.plus.minus.gene.interested$gene
    # %in% c('TPCN2')),] re.FC.sorted<-re.FC[order(re.FC$pvalue),]


    write.csv(re.DESeq.DoGs.plus.gene, file = paste0(out.dir.name, "3UTR_DE_",
        out.file.pattern.positive.gene, ".csv"), row.names = FALSE, quote = FALSE)

    write.csv(re.DESeq.DoGs.minus.gene, file = paste0(out.dir.name, "3UTR_DE_",
        out.file.pattern.negative.gene, ".csv"), row.names = FALSE, quote = FALSE)

    write.csv(re.DESeq.DoGs.plus.minus.gene, file = paste0(out.dir.name, "3UTR_DE_",
        out.file.pattern.all, ".csv"), row.names = FALSE, quote = FALSE)
    # Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)

    write.csv(re.DESeq.DoGs.plus.minus.gene.interested, file = paste0(out.dir.name,
        "3UTR_DE_", out.file.pattern.interested, ".csv"))


    # write.csv(re.DESeq.DoGs.plus.minus.gene,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE) Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)


    re.out.3 <- list(re.out = re.out, re.8.samples = re.8.samples, DE_positive_gene = re.DESeq.DoGs.plus.gene,
        DE_negative_gene = re.DESeq.DoGs.minus.gene, DE = re.DESeq.DoGs.plus.minus.gene,
        DE_interested = re.DESeq.DoGs.plus.minus.gene.interested)

    return(re.out.3)

}
DEAnalysis <- function(countData, wt.index, dox.index)
{

    a <- length(wt.index)
    b <- length(dox.index)

    colData <- data.frame(condition = factor(c(rep("Dox", b), rep("WT", a))))

    dds <- DESeqDataSetFromMatrix(countData, colData, formula(~condition))

    # size.factor<-estimateSizeFactors(dds)

    # print(size.factor)

    re.DESeq <- results(DESeq(dds))

    re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
    colnames(re.FC)[7] = "FoldChange(WT/Dox)"

    txs.gene <- ReformatTxsGene()

    re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

    re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

    return(re.FC.sorted)

}

#' getutrcount
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name <- '/Volumes/Bioinformatics$/Aimin_project/UTR/NewCounts'
#' input.file.pattern <- 'count.txt'
#'
#' res <- ThreeUTR:::getutrcount(dir.name, input.file.pattern,file.path(system.file('extdata',package = 'ThreeUTR'),'sample_infor.txt'))
#'
getutrcount <- function(dir.name, input.file.pattern, sample.infor.file, group.comparision = c(1,
    2))
    {
    file.name = file.path(dir.name, dir(dir.name, recursive = TRUE, pattern = input.file.pattern))
    file.name.2 <- as.list(file.name)

    names(file.name.2) = basename(file.name)


    file.name.3 <- lapply(file.name.2, function(u)
    {
        if (!file.size(u) == 0)
        {
            re = u
        } else
        {
            re = NULL
        }
        re
    })

    filterByRmNull <- function(re.peaks.only.bed)
    {
        re.peaks.only.bed.2 <- re.peaks.only.bed[lapply(re.peaks.only.bed, length) >
            0]

        # names(re.peaks.only.bed.2)=unlist(lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2){
        # tmp=re.peaks.only.bed.2 x=tmp[[u]] path_name=dirname(x)
        # file_name=basename(x) file_name },re.peaks.only.bed.2))
        return(re.peaks.only.bed.2)
    }

    file.name.4 <- filterByRmNull(file.name.3)

    sample.infor <- read.table(sample.infor.file, header = TRUE)

    names(file.name.4) = unlist(lapply(1:length(file.name.4), function(u, file.name.4,
        sample.infor)
        {

        tmp = file.name.4
        x = tmp[[u]]
        path_name = dirname(x)
        file_name = basename(x)
        pos = gregexpr("-", file_name)
        a = pos[[1]][1] + 1
        b = pos[[1]][2] - 1
        c = substr(file_name, a, b)
        con = sample.infor[match(c, sample.infor$Sample), ]$Condition
        file_name_con = paste0(con, "-", file_name)

        file_name_con


    }, file.name.4, sample.infor))

    re.out <- lapply(file.name.4, function(u)
    {
        if (!file.size(u) == 0)
        {
            re = read.table(u, header = F)
            colnames(re) = c("Count", "GeneName")
            re
        }
    })

    names(re.out) = names(file.name.4)

    temp.name <- strsplit(names(file.name.4), split = "\\.")
    temp.name.2 <- trimws(do.call("rbind", lapply(temp.name, "[[", 1)))

    temp.name.3 <- unique(temp.name.2[, 1])

    temp.name.4 <- as.list(temp.name.3)
    names(temp.name.4) <- temp.name.3

    txs.name <- lapply(re.out, "[[", 2)
    txs.name <- unlist(txs.name)
    txs.name <- unique(as.character(txs.name))
    txs.name.count <- data.frame(Count = rep(0, length(txs.name)), GeneName = txs.name)

    # dim(txs.name.count)

    Process4EachSample <- function(one_sample, txs.name.count, re.out)
    {

        count.gene.matched <- re.out[(grep(one_sample, names(re.out)))]

        count.gene.matched.minus.gene <- count.gene.matched[(grep("minus.gene",
            names(count.gene.matched)))]
        count.gene.matched.plus.gene <- count.gene.matched[(grep("plus.gene",
            names(count.gene.matched)))]

        ReformatCount <- function(count.gene.matched.minus.gene, txs.name.count)
        {

            reformat.count.gene.matched <- lapply(count.gene.matched.minus.gene,
                function(x, txs.name.count)
                {
                  y <- txs.name.count
                  if (dim(x[which(x$GeneName %in% y$GeneName), ])[1] != 0)
                  {
                    y[match(x$GeneName, y$GeneName), ]$Count <- x$Count
                  }
                  y
                }, txs.name.count)


            reformat.count.gene.matched.2 <- cbind(apply(do.call("cbind", lapply(reformat.count.gene.matched,
                "[", 2)), 1, unique), do.call("cbind", lapply(reformat.count.gene.matched,
                "[", 1)))

            colnames(reformat.count.gene.matched.2) = c("GeneName", names(count.gene.matched.minus.gene))

            gene <- as.character(reformat.count.gene.matched.2$GeneName)

            count.DoGs.plus.read <- apply(as.data.frame(reformat.count.gene.matched.2[,
                c(grep("plus.read.DoGs.count.txt", colnames(reformat.count.gene.matched.2)))]),
                1, sum)
            count.DoGs.minus.read <- apply(as.data.frame(reformat.count.gene.matched.2[,
                c(grep("minus.read.DoGs.count.txt", colnames(reformat.count.gene.matched.2)))]),
                1, sum)


            reformat.count.gene.matched.3 <- as.data.frame(cbind(gene, count.DoGs.plus.read,
                count.DoGs.minus.read))

            return(reformat.count.gene.matched.3)
        }

        re1 <- ReformatCount(count.gene.matched.minus.gene, txs.name.count)
        re2 <- ReformatCount(count.gene.matched.plus.gene, txs.name.count)
        re <- rbind(re1, re2)

        re <- re[order(re$gene), ]

        return(re)
    }

    re.8.samples <- lapply(temp.name.4, function(u, txs.name.count, re.out)
    {
        YY <- Process4EachSample(u, txs.name.count, re.out)
        YY
    }, txs.name.count, re.out)

    ProcessEachCorner <- function(re.8.samples, i, sample.infor, group.comparision)
    {

        gene <- apply(do.call("cbind", lapply(re.8.samples, "[", 1)), 1, unique)
        count.gene.plus.read.8.samples <- cbind(gene, do.call("cbind", lapply(re.8.samples,
            "[", i)))

        colnames(count.gene.plus.read.8.samples) = c("gene", names(re.8.samples))

        df <- count.gene.plus.read.8.samples
        df <- apply(df[, -1], 2, as.numeric)
        index <- rowSums(df) > 0
        dff <- count.gene.plus.read.8.samples[index, ]
        rownames(dff) <- dff$gene
        dff <- dff[, -1]

        X <- toupper(unique(as.character(sample.infor$Condition)))

        wt.index <- grep(X[group.comparision[1]], toupper(colnames(dff)))

        dox.index <- grep(X[group.comparision[2]], toupper(colnames(dff)))

        real.index <- c(dox.index, wt.index)

        n.sample <- length(real.index)

        permutation.index <- real.index

        countData <- apply(dff[, permutation.index], 2, as.numeric)

        rownames(countData) <- rownames(dff)

        txs.gene <- ReformatTxsGene()

        xx <- txs.gene$txs_genes_DF_2

        mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

        results <- getBM(attributes = c("ucsc", "ensembl_transcript_id"), filters = "ensembl_transcript_id",
            values = xx, mart = mart)

        xxx <- results[-which(results$ucsc == ""), ]

        rownames(xxx) <- xxx$ucsc

        re.FC <- merge(countData, xxx, by = 0)

        re.FC <- countData

        return(re.FC)

    }

    re.BR <- ProcessEachCorner(re.8.samples, 2, sample.infor)  #BR
    re.TR <- ProcessEachCorner(re.8.samples, 3, sample.infor)  #TR


    # re.BR.4.plus.gene.BL.4.minus.gene <- re.BR
    # re.TR.4.plus.gene.TL.4.minus.gene <- re.TR # Get the counts for DoGs of
    # plus and minus gene DoGs.4.plus.Gene <-
    # re.TR.4.plus.gene.TL.4.minus.gene[which(re.TR.4.plus.gene.TL.4.minus.gene$strand
    # == '+'), ] DoGs.4.minus.Gene <-
    # re.BR.4.plus.gene.BL.4.minus.gene[which(re.BR.4.plus.gene.BL.4.minus.gene$strand
    # == '-'), ] DoGs.4.plus.minus.Gene <- rbind(DoGs.4.plus.Gene,
    # DoGs.4.minus.Gene) res <- list(DoGs.4.plus.Gene=DoGs.4.minus.Gene,
    # DoGs.4.minus.Gene=DoGs.4.minus.Gene,
    # DoGs.4.plus.minus.Gene=DoGs.4.plus.minus.Gene)

    res <- list(reBR = re.BR, reTR = re.TR, re.8.samples = re.8.samples)


    return(res)

}

#' Title
#'
#' @param res
#'
#' @return
#' @export
#'
#' @examples
#'
#' y <- ThreeUTR:::annotatetranscript(res$reBR)
#'
annotatetranscript <- function(res)
{

    x <- res
    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
        host = "ensembl.org")

    t2g <- getBM(attributes = c("ucsc", "ensembl_transcript_id", "strand", "ensembl_gene_id",
        "external_gene_name"), filters = "ensembl_transcript_id", values = rownames(x),
        mart = mart)

    # xx <- t2g[-which(t2g$ensembl_transcript_id==''),]

    # rownames(xx) <- xx$ensembl_transcript_id

    # xxx <- merge(x,xx,by=0)

    xxx <- t2g

    xxx
}

#' ProcessOutputFilesFromDoGsOnly
#'
#' @param dir.name
#' @param input.file.pattern
#' @param out.dir.name
#' @param out.file.pattern.interested
#' @param out.file.pattern.positive.gene
#' @param out.file.pattern.negative.gene
#' @param out.file.pattern.all
#' @param dir.name.gene.list
#' @param pattern.4.gene.list
#' @param adjust_by_batch
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/'
#'
#' input.file.pattern='*count.2.txt'
#'
#' dir.name.gene.list='/media/H_driver/2016/Ramin_azhang/'
#' pattern.4.gene.list='final_list.csv'
#'
#'
#' out.dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/'
#'
#' for (i in 1:10) {
#'
#' time.string = gsub(':', '-', gsub(' ', '_', as.character(Sys.time())))
#' out.dir.name=paste0('/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOnePermutationAt_',time.string,'/')
#' dir.create(out.dir.name)
#'
#' out.file.pattern.interested='DoGs_adjust_by_batch_interested_gene'
#' out.file.pattern.positive.gene='DoGs_adjust_by_batch_positive'
#' out.file.pattern.negative.gene='DoGs_adjust_by_batch_negative'
#' out.file.pattern.all= 'DoGs_adjust_by_batch_all'
#'
#' Re.unadjusted.adjusted<-ProcessOutputFilesFromDoGsOnly(dir.name,input.file.pattern,out.dir.name,out.file.pattern.interested,
#' out.file.pattern.positive.gene,
#' out.file.pattern.negative.gene,
#' out.file.pattern.all,
#' dir.name.gene.list,
#' pattern.4.gene.list,
#' adjust_by_batch='YES')
#'
#' }
#'
#' save.image(file=paste0(out.dir.name,'re_save_2.RData'))
#' savehistory(file=paste0(out.dir.name,'re_save_2.Rhistory'))

ProcessOutputFilesFromDoGsOnly <- function(dir.name, input.file.pattern, out.dir.name,
    out.file.pattern.interested, out.file.pattern.positive.gene, out.file.pattern.negative.gene,
    out.file.pattern.all, dir.name.gene.list, pattern.4.gene.list, adjust_by_batch,
    permutation.set.up = NULL)
    {

    x <- getutrcount(dir.name, input.file.pattern)

    DoGs.4.plus.Gene <- x[[1]]
    DoGs.4.minus.Gene <- x[[2]]
    DoGs.4.plus.minus.Gene <- x[[3]]

    GeneTypeBasedDE <- function(DoGs.4.plus.Gene)
    {

        # a <- n+9-1

        Count.DoGs.4.plus.Gene <- DoGs.4.plus.Gene
        rownames(Count.DoGs.4.plus.Gene) <- Count.DoGs.4.plus.Gene$Row.names
        Count.DoGs.4.plus.Gene.2 <- Count.DoGs.4.plus.Gene[, -1]

        # head(Count.DoGs.4.plus.Gene.2)

        wt.index <- grep("WT", toupper(colnames(Count.DoGs.4.plus.Gene.2)))

        dox.index <- grep("DOX", toupper(colnames(Count.DoGs.4.plus.Gene.2)))

        real.index <- c(dox.index, wt.index)

        n <- length(real.index)

        if (!is.null(permutation.set.up))
        {
            real.index <- seq(1, n)
            permutation.index <- real.index
            permutation.index = array(sample(real.index))
        } else
        {
            permutation.index <- real.index
        }

        tmp <- Count.DoGs.4.plus.Gene.2

        Count.DoGs.4.plus.Gene.2 <- tmp[, permutation.index]

        cat("get count\n")
        print(head(Count.DoGs.4.plus.Gene.2))
        cat("get count done \n")

        if (adjust_by_batch == "NO")
        {
            re.DESeq.DoGs.plus.gene <- DEAnalysis(Count.DoGs.4.plus.Gene.2,
                wt.index, dox.index)
        } else
        {

            re.DESeq.DoGs.plus.gene <- DEAnalysisAdjustByBatch(Count.DoGs.4.plus.Gene.2,
                Nbatch = c(3, 2))
        }

        return(re.DESeq.DoGs.plus.gene)
    }

    print(head(DoGs.4.plus.Gene))
    re.DESeq.DoGs.plus.gene <- GeneTypeBasedDE(DoGs.4.plus.Gene)
    # head(re.DESeq.DoGs.plus.gene[which(re.DESeq.DoGs.plus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.gene$pvalue<0.05),])

    print(head(DoGs.4.minus.Gene))
    re.DESeq.DoGs.minus.gene <- GeneTypeBasedDE(DoGs.4.minus.Gene)
    # head(re.DESeq.DoGs.minus.gene[which(re.DESeq.DoGs.minus.gene$log2FoldChange<0&re.DESeq.DoGs.minus.gene$pvalue<0.05),])

    print(head(DoGs.4.plus.minus.Gene))
    re.DESeq.DoGs.plus.minus.gene <- GeneTypeBasedDE(DoGs.4.plus.minus.Gene,
        n)
    # head(re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.minus.gene$pvalue<0.05),])

    gene.interested <- ReadGeneList(dir.name.gene.list, pattern.4.gene.list)

    re.DESeq.DoGs.plus.minus.gene.interested <- re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$gene %in%
        gene.interested[, 1]), ]

    # re.DESeq.DoGs.plus.minus.gene.interested[which(re.DESeq.DoGs.plus.minus.gene.interested$gene
    # %in% c('TPCN2')),] re.FC.sorted<-re.FC[order(re.FC$pvalue),]

    ReformatGeneSymbol <- function(re.DoGs.adjusted.by.batch)
    {
        re2 <- re.DoGs.adjusted.by.batch
        re2$gene <- as.character(re2$gene)
        re2$gene <- paste0("'", as.character(re2$gene), "'")
        return(re2)
    }

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.positive.gene,
        ".csv")), row.names = FALSE, quote = FALSE)

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.minus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.negative.gene,
        ".csv")), row.names = FALSE, quote = FALSE)

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.all,
        ".csv")), row.names = FALSE, quote = FALSE)
    # Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)
    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene.interested)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.interested,
        ".csv")))

    # write.csv(re.DESeq.DoGs.plus.minus.gene,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE) Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)


    re.out.3 <- list(dogcount = x, DE_positive_gene = re.DESeq.DoGs.plus.gene,
        DE_negative_gene = re.DESeq.DoGs.minus.gene, DE = re.DESeq.DoGs.plus.minus.gene,
        DE_interested = re.DESeq.DoGs.plus.minus.gene.interested)

    return(re.out.3)

}




idtransform <- function(values_to_be_transforemd)
{

    mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
        host = "ensembl.org")

    t2g <- getBM(attributes = c("ucsc", "ensembl_transcript_id", "strand", "ensembl_gene_id",
        "external_gene_name"), filters = "ucsc", values = values_to_be_transforemd,
        mart = mart)
    t2g


}

#' Title
#'
#' @param input.file
#' @param output.file
#'
#' @return
#' @export
#'
#' @examples
#' input.file <- '/Volumes/Bioinformatics$/Zhen-Gao/UTR/Gencode-3UTR-transcript.tsv'
#' output.file <- '/Volumes/Bioinformatics$/Zhen-Gao/UTR/Gencode-3UTR-transcript-transformed.tsv'
#' ThreeUTR:::gaoidtransform(input.file,output.file)
#'
gaoidtransform <- function(input.file, output.file)
{

    res <- read.table(input.file, header = F)
    colnames(res) <- c("Chr", "Start", "End", "utr", "Symbol", "Strand", "ucsc")

    # head(res)

    y <- ThreeUTR:::idtransform(res$ucsc)

    yy <- merge(res, y, by = "ucsc")
    write.table(yy, file = output.file, append = FALSE, quote = FALSE, sep = "\t",
        eol = "\n", na = "NA", dec = ".", row.names = FALSE)
}
#' ProcessTranscript2
#'
#' @return
#' @export
#'
#' @examples
#'
#' ProcessTranscript2()
ProcessTranscript2 <- function()
{

    # Gene based
    genes_hg19 <- as.data.frame(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

    # Transcripts based
    txs.genes <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)

    gene.symbol <- getSYMBOL(names(txs.genes), data = "org.Hs.eg")

    txs.genes.2 <- txs.genes
    names(txs.genes.2) <- gene.symbol

    gene.symbol.2 <- getSYMBOL(genes_hg19$gene_id, data = "org.Hs.eg")
    gene.hg19.gene.symbol <- cbind(genes_hg19, gene.symbol.2)

    txs.genes.3 <- lapply(txs.genes.2, as.data.frame)

    txs.genes.3.name <- sapply(txs.genes.3, function(x)
    {

        gene <- names(x)
        # })

        # print(x)

        n.txs <- dim(x)[1]
        y <- cbind(x, rep(gene, n.txs))
        y
    }, simplify = FALSE, USE.NAMES = TRUE)

    ListReformat <- function(txs.genes.3)
    {
        txs.genes.3.name <- lapply(seq_along(txs.genes.3), function(i)
        {
            name.gene <- names(txs.genes.3)[[i]]
            n.txs <- dim(txs.genes.3[[i]])[1]
            gene <- rep(name.gene, n.txs)
            y <- cbind(txs.genes.3[[i]], gene)
            y
        })
        names(txs.genes.3.name) <- names(txs.genes.3)
        return(txs.genes.3.name)
    }

    txs.genes.3.name.2 <- ListReformat(txs.genes.3)

    txs.genes.3.name.2.DF <- do.call(rbind.data.frame, txs.genes.3.name.2)

    # txs.genes.3.name.22<-unlist2(txs.genes.3.name.2)

    txs.genes.33 <- do.call(rbind.data.frame, txs.genes.3)
    txs.genes.33.sorted <- txs.genes.33[order(row.names(txs.genes.33)), ]

    # txs.genes.33.sorted.2<-cbind(txs.genes.33.sorted,strsplit(rownames(txs.genes.33.sorted),split='.'))

    # txs.gene.333<-data.frame(Reduce(rbind,txs.genes.3))

    print(txs.genes.3[which(names(txs.genes.3) == "STRA6")])

    txs.genes.4 <- lapply(txs.genes.3, function(x)
    {
        y <- unique(x$seqnames)
        ny <- length(y)
        ny
    })
    txs.genes.5 <- do.call("rbind", lapply(txs.genes.4, "[[", 1))

    head(txs.genes.5)

    gene.not.unique.chr <- unique(names(txs.genes.5[which(txs.genes.5[, 1] !=
        1), ]))

    gene.unique.chr <- unique(names(txs.genes.5[which(txs.genes.5[, 1] == 1),
        ]))


    txs.genes.3.name.not.unique <- txs.genes.3.name.2[which(names(txs.genes.3.name.2) %in%
        gene.not.unique.chr)]

    txs.genes.3.name.unique <- txs.genes.3.name.2[which(names(txs.genes.3.name.2) %in%
        gene.unique.chr)]


    txs.genes.4.strand <- lapply(txs.genes.3.name.unique, function(x)
    {
        y <- unique(x$strand)
        ny <- length(y)
        ny
    })
    txs.genes.5.strand <- do.call("rbind", lapply(txs.genes.4.strand, "[[",
        1))

    temp <- names(txs.genes.5.strand[which(txs.genes.5.strand[, 1] == 2), ])

    txs.genes.3.unique.name.not.strand <- txs.genes.3.name.unique[which(names(txs.genes.3.name.unique) %in%
        temp)]

    length(unique(names(txs.genes.3.unique.name.not.strand)))

    temp2 <- names(txs.genes.5.strand[which(txs.genes.5.strand[, 1] == 1), ])

    txs.genes.3.unique.name.unique.strand <- txs.genes.3.name.unique[which(names(txs.genes.3.name.unique) %in%
        temp2)]

    length(unique(names(txs.genes.3.unique.name.unique.strand)))


    input.dir = "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"
    input.file.pattern = "final_list.csv"
    file.name = paste0(input.dir, dir(input.dir, recursive = TRUE, pattern = input.file.pattern))
    genes.interested <- read.csv(file.name, header = F)

    length(which(names(txs.genes.3.unique.name.unique.strand) %in% genes.interested[,
        1]))

    length(which(names(txs.genes.3.unique.name.not.strand) %in% genes.interested[,
        1]))

    length(which(names(txs.genes.3.name.not.unique) %in% genes.interested[,
        1]))

    txs.genes.3.name.not.unique[which(names(txs.genes.3.name.not.unique) %in%
        genes.interested[, 1])]

    print(gene.hg19.gene.symbol[which(gene.hg19.gene.symbol$gene.symbol.2 ==
        "STRA6"), ])

    re <- list(txs_genes_list = txs.genes.3.name.2, txs_genes_DF = txs.genes.3.name.2.DF)

    return(re)

}
#' Title
#'
#' @param dir.name
#' @param input.file.pattern.4.interested.gene
#'
#' @return
#' @export
#'
#' @examples
#'  dir.name='/media/H_driver/2016/Ramin_azhang/'
#'  input.file.pattern.4.interested.gene='final_list.csv'
#' gene.interested<-ReadGeneList(dir.name, input.file.pattern.4.interested.gene)
#'
ReadGeneList <- function(dir.name, input.file.pattern.4.interested.gene)
{
    file.name = file.path(dir.name, input.file.pattern.4.interested.gene)
    genes.interested <- read.csv(file.name, header = F)
    return(genes.interested)
}
#' ReformatTxsGene
#'
#' @return
#' @export
#'
#' @examples
#'
#' ReformatTxsGene()
#'
ReformatTxsGene <- function()
{

    # Transcripts based
    txs.genes <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
    gene.symbol <- getSYMBOL(names(txs.genes), data = "org.Hs.eg")

    txs.genes.2 <- txs.genes
    names(txs.genes.2) <- gene.symbol

    txs.genes.3 <- lapply(txs.genes.2, as.data.frame)

    ListReformat <- function(txs.genes.3)
    {
        txs.genes.3.name <- lapply(seq_along(txs.genes.3), function(i)
        {
            name.gene <- names(txs.genes.3)[[i]]
            n.txs <- dim(txs.genes.3[[i]])[1]
            gene <- rep(name.gene, n.txs)
            y <- cbind(txs.genes.3[[i]], gene)
            y
        })
        names(txs.genes.3.name) <- names(txs.genes.3)
        return(txs.genes.3.name)
    }

    txs.genes.3.name.2 <- ListReformat(txs.genes.3)

    txs.genes.3.name.2.DF <- do.call(rbind.data.frame, txs.genes.3.name.2)

    txs.genes.3.name.2.DF.2 <- txs.genes.3.name.2.DF

    rownames(txs.genes.3.name.2.DF.2) <- txs.genes.3.name.2.DF$tx_name

    re <- list(txs_genes_list = txs.genes.3.name.2, txs_genes_DF = txs.genes.3.name.2.DF,
        txs_genes_DF_2 = txs.genes.3.name.2.DF.2)

    return(re)

}
Sleuth <- function(dplyr, select, run_accession, condition, mutate, biomaRt,
    useMart, getBM, rename, ensembl_transcript_id, ensembl_gene_id, external_gene_name)
    {
    base_dir <- "~/Downloads/"

    sample_id <- dir(file.path(base_dir, "results"))

    kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results",
        id, "kallisto"))

    s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE,
        stringsAsFactors = FALSE)
    s2c <- dplyr::select(s2c, sample = run_accession, condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    print(s2c)
    so <- sleuth_prep(s2c, ~condition)
    so <- sleuth_fit(so)
    so <- sleuth_fit(so, ~1, "reduced")
    so <- sleuth_lrt(so, "reduced", "full")
    models(so)
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
        host = "ensembl.org")

    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
        "external_gene_name"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id,
        ext_gene = external_gene_name)
    so <- sleuth_prep(s2c, ~condition, target_mapping = t2g)

    so <- sleuth_fit(so)

    so <- sleuth_fit(so, ~1, "reduced")

    so <- sleuth_lrt(so, "reduced", "full")

    sleuth_live(so)
    results_table <- sleuth_results(so, "reduced:full", test_type = "lrt")

    # gene level
    so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = "ens_gene")
}


TxsDE <- function(dplyr, select, run_accession, condition, mutate, biomaRt,
    useMart, getBM, rename, ensembl_transcript_id, ensembl_gene_id, external_gene_name)
    {
    base_dir <- "~/Downloads/"

    sample_id <- dir(file.path(base_dir, "results"))

    kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results",
        id, "kallisto"))

    s2c <- read.table(file.path(base_dir, "hiseq_info.txt"), header = TRUE,
        stringsAsFactors = FALSE)
    s2c <- dplyr::select(s2c, sample = run_accession, condition)
    s2c <- dplyr::mutate(s2c, path = kal_dirs)
    print(s2c)

    # transcript based
    so <- sleuth_prep(s2c, ~condition)
    so <- sleuth_fit(so)
    so <- sleuth_fit(so, ~1, "reduced")
    so <- sleuth_lrt(so, "reduced", "full")
    models(so)

    # Including gene names into transcript-level analysis
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
        host = "ensembl.org")

    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
        "external_gene_name"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id,
        ext_gene = external_gene_name)


    so <- sleuth_prep(s2c, ~condition, target_mapping = t2g)

    so <- sleuth_fit(so)

    so <- sleuth_fit(so, ~1, "reduced")

    so <- sleuth_lrt(so, "reduced", "full")

    sleuth_live(so)
    results_table <- sleuth_results(so, "reduced:full", test_type = "lrt")

    # gene level based
    so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = "ens_gene")
    so <- sleuth_fit(so)
    so <- sleuth_fit(so, ~1, "reduced")
    so <- sleuth_lrt(so, "reduced", "full")
    models(so)

    # get results
    results_table_gene <- sleuth_results(so, "reduced:full", test_type = "lrt")

}
#' TestUseLongestTxsOfGene
#'
#' @param Re.unadjusted.adjusted
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' out.dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts5CasesEachSample/'
#' out.file.pattern='DoGs_using_longest_txs_adjust_by_batch'
#'
#' re.longest.transcript<-TestUseLongestTxsOfGene(Re.unadjusted.adjusted,out.dir.name,out.file.pattern)
#'
TestUseLongestTxsOfGene <- function(Re.unadjusted.adjusted, out.dir.name, out.file.pattern)
{

    tempData <- Re.unadjusted.adjusted$DE

    # which(tempData$gene==NA)

    tempdata.byGSym.2 <- tempData

    rownames(tempdata.byGSym.2) = NULL

    data.byGSym = ddply(tempdata.byGSym.2, c("gene"), function(h)
    {
        # summary = apply(h,2,max)
        y <- as.data.frame(h[which.max(h$width), ])
        y
    })

    df.NT <- data.byGSym[, c(1, 9:16)]
    rownames(df.NT) <- df.NT[, 1]
    df.NT.2 <- df.NT[, -1]

    re.DoGs.adjusted.by.batch <- DEAnalysisAdjustByBatch(df.NT.2)

    re2 <- re.DoGs.adjusted.by.batch
    re2$gene <- as.character(re2$gene)
    re2$gene <- paste0("'", as.character(re2$gene), "'")

    re.DoGs.adjusted.by.batch[which(re.DoGs.adjusted.by.batch$Row.names == "uc010paz.2"),
        ]

    # length(unique(as.character(re.DoGs.adjusted.by.batch$gene)))
    # unique(as.character(re.DoGs.adjusted.by.batch$gene))
    # [grep('uc002jah.2',unique(as.character(re.DoGs.adjusted.by.batch$Row.names)))]

    write.csv(re2, file = paste0(out.dir.name, "3UTR_DE_", out.file.pattern,
        ".csv"), row.names = FALSE, quote = TRUE)

    return(re.DoGs.adjusted.by.batch)

}
installercran <- function(.cran_packages)
{
    .inst <- .cran_packages %in% installed.packages()
    if (any(!.inst))
    {
        install.packages(.cran_packages[!.inst])
    }
    # Load packages into session, and print package version
    sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
}

installerbioc <- function(.bioc_packages)
{
    .inst <- .bioc_packages %in% installed.packages()
    if (any(!.inst))
    {
        source("http://bioconductor.org/biocLite.R")
        biocLite(.bioc_packages[!.inst], ask = F)
    }
    # Load packages into session, and print package version
    sapply(.bioc_packages, require, character.only = TRUE)
}

parsersample <- function()
{
    cell <- factor(rep(c("batch1", "batch2"), c(3, 2)))
    cell = rep(cell, 2)

    sample <- c("R1_Dox.bam", "R2_Dox.bam", "R3_Dox.bam", "R4_Dox.bam", "R5_Dox.bam",
        "R1_WT.bam", "R2_WT.bam", "R3_WT.bam", "R4_WT.bam", "R5_WT.bam")

    colData <- data.frame(sample = sample, condition = factor(rep(c("Dox", "WT"),
        c(5, 5))), cell = cell)

    colData
}

#' parserreadfiles
#'
#' @param input.file.dir
#' @param input.file.type
#' @param output.file.dir
#'
#' @return
#' @export
#'
#' @examples
#' input.file.dir <- '/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq'
#' output.file.dir <- '/Volumes/Bioinformatics$/Aimin_project'
#'
#'
#' res <- parserreadfiles(input.file.dir,'bam',filter.sample='emp1')
#'
parserreadfiles <- function(input.file.dir, input.file.type, sample.group = NULL,
    filter.sample = NULL)
    {
    dir.name = input.file.dir
    dir.name = reformatPath(dir.name)
    file.name = file.path(dir.name, dir(dir.name, recursive = TRUE))
    file.name.2 <- as.list(file.name)
    file.name.3 <- lapply(1:length(file.name.2), function(u, input.file.type,
        file.name.2)
        {
        tmp = file.name.2
        x = tmp[[u]]
        path_name = dirname(x)
        file_name = basename(x)
        n <- length(grep(input.file.type, file_name))
        if (n == 1)
        {
            if (file_ext(file_name) == input.file.type)
            {
                re <- file.path(path_name, file_name)
            } else
            {
                re <- NULL
            }

        } else
        {
            re <- NULL
        }
        re
    }, input.file.type, file.name.2)

    file.name.4 <- file.name.3[lapply(file.name.3, length) > 0]
    names(file.name.4) = unlist(lapply(1:length(file.name.4), function(u, file.name.4)
    {
        tmp = file.name.4
        x = tmp[[u]]
        path_name = dirname(x)
        path_name2 <- basename(path_name)

        file_name = basename(x)

        file_name <- paste0(path_name2, "-", file_name)

        file_name

    }, file.name.4))

    if (!is.null(filter.sample))
    {
        file.name.5 <- file.name.4[-grep(filter.sample, names(file.name.4))]

    } else
    {
        file.name.5 <- file.name.4
    }


    if (!is.null(sample.group))
    {
        x <- unlist(file.name.5)

        file.name.6 <- lapply(1:length(sample.group), function(u,sample.group,
            x)
            {

            xx <- x[grep(toupper(sample.group[[u]]),toupper(x))]
            xx

        }, sample.group, x)

        # g1 <- grep(sample.group[1], toupper(names(res$input)))

        # g2 <- grep(sample.group[2], toupper(names(res$input)))

        # output.dir.name = reformatPath(output.file.dir) temp3 = output.dir.name
        re2 <- list(input = file.name.6, input.file.type = input.file.type,sample.group=sample.group)
    } else
    {
        re2 <- list(input = file.name.5, input.file.type = input.file.type)
    }

    # pkg.env <- new.env(parent = emptyenv()) pkg.env$sample <-
    # ThreeUTR:::parsersample()

    return(re2)


}

#' Example:
#' R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useWget2Download("SRP058633","/nethome/axy148/DoGsExample")'
#'
#'re <- ThreeUTR:::useWget2Download("SRP058633","/nethome/axy148/DoGsExample")
#'
useWget2Download <- function(sra.accession.number, output.dir)
{
    cmd0 <- "wget -c -r -nd -np -L"
    cmd1 <- "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/"

    temp <- file.path(substr(sra.accession.number, 1, 6), sra.accession.number)

    temp2 <- paste0(cmd1, temp)

    cmd2 <- paste(cmd0, temp2, "-P", output.dir, sep = " ")

    m.id <- grep("login", system("hostname", intern = TRUE))

    if(length(m.id) == 1){
    job.name <- "wgetDownload"
    cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
                                  cores = 32, Memory = 25000, span.ptile = 16, job.name)
    #cmd.end.check = "&& echo 'done wget download' > 'Download.txt'"
    cmd3 <- paste(cmd.p,cmd2,sep = " ")
    }else
    {cmd3 <- cmd2}

    cmd3
    system(cmd3,intern = TRUE)

}

useFastqDump <- function(sra.accession.number, output.dir)
{
    cmd0 <- "fastq-dump -I --split-files"

    cmd1 <- paste(cmd0, sra.accession.number, "-O", output.dir, sep = " ")

    system(cmd1)

}

#'Example
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useFastqDumpConvertSra2Fastq("/nethome/axy148/DoGsExample","/scratch/projects/bbc/aiminy_project/DoGsFastq")'
#'
#'On pegasus, you need to run this if you want to start conversion automatically after downloading using wget
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useFastqDumpConvertSra2Fastq("/nethome/axy148/DoGsExample","/scratch/projects/bbc/aiminy_project/DoGsFastq",wait.job.name = "wgetDownload")'
#'
useFastqDumpConvertSra2Fastq <- function(sra.file.dir, output.dir,wait.job.name=NULL)
{

    re <- parserreadfiles(sra.file.dir, "sra")

    res <- re$input

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir,recursive = TRUE)
    }

    m.id <- grep("login", system("hostname", intern = TRUE))

    cmd.l <- lapply(1:length(res), function(u, res,m.id,wait.job.name,output.dir)
    {
        cmd0 <- "fastq-dump --split-3"

        path_name = dirname(res[[u]])

        file_name = file_path_sans_ext(basename(res[[u]]))

        cmd1 <- paste(cmd0, res[[u]], "-O", output.dir, sep = " ")

        if(length(m.id) == 1){
          job.name <-  paste0("sra2fastq.",u)

          if(!is.null(wait.job.name)){
          wait.job.name <- wait.job.name
          cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
                                        cores = 32, Memory = 25000, span.ptile = 16, job.name=job.name,wait.job.name=wait.job.name)
           }else
           {
             cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
                                           cores = 32, Memory = 25000, span.ptile = 16, job.name=job.name)
           }
          cmd2 <- paste(cmd.p,cmd1,sep = "")
        }else
        {cmd2 <- cmd1}

        system(cmd2)
        cmd2
    }, res,m.id,wait.job.name,output.dir)

    re <- list(cmdl = cmd.l, output.dir = output.dir)

    re

}

# infer_experiment.py
useInferExperiment <- function(input.file.dir, ref.gene.bed.file, output.dir)
{
    re <- parserreadfiles(input.file.dir, "bam")

    res <- re$input

    cmd0 <- "infer_experiment.py -i"
    cmd1 <- "-r"
    cmd2 <- ">"

    # output.dir <- file.path(output.dir, 'BamInfo')

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    cmd.l <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd2 <- paste(cmd0, u, cmd1, ref.gene.bed.file, cmd2, file.path(output.dir,
            paste0(file_name, "_infor.txt")), sep = " ")

        system(cmd2)

        cmd2
    }, output.dir)

    re <- list(cmdl = cmd.l, output.dir = output.dir)

    re
}

#' convertbam2bed
#'
#' @param input.bamfile.dir
#' @param output.bedfile.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' input.bamfile.dir <- '/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq'
#' output.bedfile.dir <- '/Volumes/Bioinformatics$/Aimin_project'
#'
#' res <- convertbam2bed(input.bamfile.dir,output.bedfile.dir)
#'
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::convertbam2bed('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs')'
#'
convertbam2bed <- function(input.bamfile.dir, output.bedfile.dir)
{
    res <- parserreadfiles(input.bamfile.dir, "bam")

    res <- res$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    # if (!dir.exists(output.bedfile.dir)) { dir.create(output.bedfile.dir,
    # recursive = TRUE) }

    output.bedfile.dir <- file.path(output.bedfile.dir, "BedFileFromBam")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bedfile.dir)
    {
        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11
        # <- '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13
        # <- paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)

        path_name = dirname(res[[u]])
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        file_name <- paste0(path_name2, "-", file_name)
        if (m.id == 1)
        {
            job.name <- paste0("bam2bed.", u)
            cmd0 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 25000, span.ptile = 16, job.name)

            cmd1 <- "bedtools bamtobed -i"
            cmd2 <- "\\>"

            cmd3 <- paste(cmd0, cmd1, res[[u]], cmd2, file.path(output.bedfile.dir,
                paste0(file_name, ".bed")), sep = " ")
        } else
        {
            cmd3 <- paste(cmd1, res[[u]], cmd2, file.path(output.bedfile.dir,
                paste0(file_name, ".bed")), sep = " ")
        }

        system(cmd3)

        cmd3
    }, m.id, res, output.bedfile.dir)

    re <- list(cmdl = cmd.l, output.bedfile.dir = output.bedfile.dir)

    re

}

#' matchbed2annotation
#'
#' @param input.bedfile.dir
#' @param annotation.bed.file
#' @param ld upstream base
#' @param rd downstream base
#' @param output.matched.bed.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' res <- convertbam2bed(input.bamfile.dir,output.bedfile.dir)
#' input.bedfile.dir <- res$output.bedfile.dir
#' annotation.bed.file <- annotation.bed.file
#' ld <- ld
#' rd <- rd
#' output.matched.bed.file.dir <- output.matched.bed.file.dir
#'
#' res <- matchbed2annotation(input.bedfile.dir,annotation.bed.file,
#' ld,rd,output.matched.bed.file.dir)
#'
matchbed2annotation <- function(input.bedfile.dir, annotation.bed.file, ld,
    rd, output.matched.bed.file.dir)
    {
    res <- parserreadfiles(input.bedfile.dir, "bed")

    res <- res$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    output.bedfile.dir <- file.path(output.matched.bed.file.dir, "MatchedBedFile")

    if (!dir.exists(output.bedfile.dir))
    {
        dir.create(output.bedfile.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, res, m.id, ld, rd, annotation.bed.file,
        output.bedfile.dir)
        {
        # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11
        # <- '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13
        # <- paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)
        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            job.name <- paste0("bedMannot.", u)
            wait.job.name <- paste0("bam2bed.", u)
            cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 25000, span.ptile = 16, job.name, wait.job.name)

            # cmd1 <- 'bedtools bamtobed -i' cmd2 <- '\\>'

            cmd0 <- paste("bedtools window -a", annotation.bed.file, "-b", sep = " ")
            cmd1 <- paste("-l", ld, "-r", rd, "-sw", "\\>", sep = " ")
            cmd2 <- paste(cmd.p, cmd0, res[[u]], cmd1, file.path(output.bedfile.dir,
                paste0(file_name, "_matched.bed")), sep = " ")

        } else
        {
            cmd0 <- paste("bedtools window -a", annotation.bed.file, "-b", sep = " ")
            cmd1 <- paste("-l", ld, "-r", rd, "-sw", ">", sep = " ")
            cmd2 <- paste(cmd.p, cmd0, res[[u]], cmd1, file.path(output.bedfile.dir,
                paste0(file_name, "_matched.bed")), sep = " ")
        }

        system(cmd2)

        cmd2
    }, res, m.id, ld, rd, annotation.bed.file, output.bedfile.dir)

    re <- list(cmdl = cmd.l, output.bedfile.dir = output.bedfile.dir)

    re

}

#' getcountsfromMatchedbed
#'
#' @param input.bedfile.dir
#' @param annotation.bed.file
#' @param ld upstream base
#' @param rd downstream base
#' @param output.matched.bed.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' res <- matchbed2annotation(input.bedfile.dir,annotation.bed.file,
#' ld,rd,output.matched.bed.file.dir)
#'
#' input.bedfile.dir <- res$output.bedfile.dir
#'
#' res <- getcountsfromMatchedbed (input.bedfile.dir,output.count.file.dir)
#'
getcountsfromMatchedbed <- function(input.bedfile.dir, output.count.file.dir,
    filter.sample)
    {

    # bsub -P bbc -J 'count.8' -o %J.count.8.log -e %J.count.8.err -W 72:00 -n
    # 32 -q parallel -R 'rusage[mem= 25000 ] span[ptile= 16 ]' -u
    # aimin.yan@med.miami.edu 'awk -F '\t' '\$6==\'-\' && \$12==\'-\''
    # /scratch/projects/bbc/aiminy_project/DoGs/MatchedBedFile/BAM-SRR2039089-Fs-accepted_hits_matched.bed
    # | awk '\$8 <= \$3 && \$9 > \$3' | awk '{print \$4}' | sort | uniq -c
    # | sort -nr >
    # /scratch/projects/bbc/aiminy_project/DoGs/Counts/BAM-SRR2039089-Fs-accepted_hits_matched.minus.gene.minus.read.below.DoGs.count.txt'

    res <- parserreadfiles(input.bedfile.dir, "bed", filter.sample = filter.sample)

    res <- res$input

    # system('awk -F '\\t' '$6==\'+\' && $12==\'-\''
    # ~/MatchedBedFile/R1_Dox_matched.bed | awk '$8<$2&&$9>=$2' | awk '{print
    # $4}' | sort | uniq -c | sort -nr | head')
    m.id <- grep("login", system("hostname", intern = TRUE))
    cmd0 <- "awk -F '\\t'"

    cmd1 <- "'\\$6==\\\"+\\\" && \\$12==\\\"-\\\"'"
    cmd11 <- "'\\$6==\\\"+\\\" && \\$12==\\\"+\\\"'"
    cmd12 <- "'\\$6==\\\"-\\\" && \\$12==\\\"+\\\"'"
    cmd13 <- "'\\$6==\\\"-\\\" && \\$12==\\\"-\\\"'"

    cmd2 <- "| awk '\\$8 < \\$2 && \\$9 >= \\$2' | awk '{print \\$4}' | sort | uniq -c | sort -nr"  #below
    cmd21 <- "| awk '\\$8 >= \\$2 && \\$9 <= \\$3' | awk '{print \\$4}' | sort | uniq -c | sort -nr"  #DoGs
    cmd22 <- "| awk '\\$8 <= \\$3 && \\$9 > \\$3' | awk '{print \\$4}' | sort | uniq -c | sort -nr"  #over

    cmd3 <- ">"

    cmdtemp <- rbind(cbind(rep(cmd1, 3), c(cmd2, cmd21, cmd22), rep("plus",
        3), rep("minus", 3), c("below.DoGs", "DoGs", "over.DoGs")), cbind(rep(cmd11,
        3), c(cmd2, cmd21, cmd22), rep("plus", 3), rep("plus", 3), c("below.DoGs",
        "DoGs", "over.DoGs")), cbind(rep(cmd12, 3), c(cmd2, cmd21, cmd22), rep("minus",
        3), rep("plus", 3), c("over.DoGs", "DoGs", "below.DoGs")), cbind(rep(cmd13,
        3), c(cmd2, cmd21, cmd22), rep("minus", 3), rep("minus", 3), c("over.DoGs",
        "DoGs", "below.DoGs")))

    output.count.file.dir <- file.path(output.count.file.dir, "Counts")

    if (!dir.exists(output.count.file.dir))
    {
        dir.create(output.count.file.dir, recursive = TRUE)
    }

    counteachcase <- function(res, m.id, cmd0, cmd1, cmd2, cmd3, gene.strand,
        read.strand, location, output.count.file.dir)
        {
        cmd.l <- lapply(1:length(res), function(u, m.id, res, cmd0, cmd1, cmd2,
            cmd3, gene.strand, read.strand, location, output.count.file.dir)
            {
            file_name = file_path_sans_ext(basename(res[[u]]))

            if (m.id == 1)
            {
                job.name <- paste0("count.", u)
                wait.job.name <- paste0("bedMannot.", u)

                cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
                  cores = 32, Memory = 25000, span.ptile = 16, job.name)
                cmd3 <- ">"

                cmd4 <- paste(cmd.p, paste0("\"", cmd0), cmd1, res[[u]], cmd2,
                  cmd3, file.path(output.count.file.dir, paste0(file_name, ".",
                    gene.strand, ".gene.", read.strand, ".read.", location,
                    ".count.txt", "\"")), sep = " ")

            } else
            {
                cmd4 <- paste(cmd0, cmd1, res[[u]], cmd2, cmd3, file.path(output.count.file.dir,
                  paste0(file_name, ".", gene.strand, ".gene.", read.strand,
                    ".read.", location, ".count.txt")), sep = " ")

            }

            cat(cmd4, "\n")

            system(cmd4, intern = TRUE)

            cmd4

        }, m.id, res, cmd0, cmd1, cmd2, cmd3, gene.strand, read.strand, location,
            output.count.file.dir)

        return(cmd.l)

    }

    cmdtempres2 <- apply(cmdtemp, 1, function(u, res, m.id, cmd0, cmd3, output.count.file.dir)
    {
        x <- as.data.frame(t(u))

        cmd1 <- x[, 1]
        cmd2 <- x[, 2]
        gene.strand <- x[, 3]

        read.strand <- x[, 4]

        location <- x[, 5]

        cmdtempres <- counteachcase(res, m.id, cmd0, cmd1, cmd2, cmd3, gene.strand,
            read.strand, location, output.count.file.dir)

        cmdtempres

    }, res, m.id, cmd0, cmd3, output.count.file.dir)

    re <- list(cmdl = cmdtempres2, output.count.file.dir = output.count.file.dir)

    re

}

#' getcounts
#'
#' @param input.bamfile.dir
#' @param annotation.bed.file
#' @param ld
#' @param rd
#' @param output.count.file.dir
#'
#' @return
#' @export
#'
#' @examples
#'
#' getcounts()
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::getcounts('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/projects/ctsi/bbc/aimin/annotation/hg19_DoGs_2.bed',0,0,'/scratch/projects/bbc/aiminy_project/DoGs')'
#'
getcounts <- function(input.bamfile.dir, annotation.bed.file, ld, rd, output.count.file.dir,
    filter.sample)
    {
    res <- convertbam2bed(input.bamfile.dir, output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    annotation.bed.file <- annotation.bed.file
    ld <- ld
    rd <- rd

    res <- matchbed2annotation(input.bedfile.dir, annotation.bed.file, ld, rd,
        output.count.file.dir)

    input.bedfile.dir <- res$output.bedfile.dir

    res <- getcountsfromMatchedbed(input.bedfile.dir, output.count.file.dir,
        filter.sample)

    res
}

#' parserAnnotationFile
#'
#' @return
#' @export
#'
#' @examples
#'
#' input.annotation.file <- '/Volumes/Bioinformatics$/Aimin_project/UTR/3UTR_GFF2.txt'
#'
#' y <- parserAnnotationFile(input.annotation.file)
#'
#'
parserAnnotationFile <- function(input.annotation.file)
{
    dir.name <- dirname(input.annotation.file)

    file.name <- file_path_sans_ext(basename(input.annotation.file))

    x <- read.table(input.annotation.file)

    d <- x[, 5] - x[, 4]

    nn <- substr(x[, 9], 9, 23)

    xx <- cbind.data.frame(x, nn, d)

    counts <- as.data.frame(table(xx$nn))

    colnames(counts) = c("nn", "counts")

    xxx <- merge(xx, counts, by = "nn", sort = FALSE)

    xxxx <- xxx[, c(2, 5, 6, 10, 11, 8)]

    write.table(xxxx, file = file.path(dir.name, paste0(file.name, ".bed")),
        row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

    xxxx
}

convertBam2Bw <- function(input.bam.file.dir, input.chromosome.size.file, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    cmd0 = "samtools sort"

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd <- paste(cmd0, u, file.path(output.bw.file.dir, paste0(file_name,
            "_sorted")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd1 = "samtools index"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        path_name = dirname(u)
        path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        file_name <- paste0(path_name2, "-", file_name)

        cmd <- paste(cmd1, file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam")),
            sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd3 <- "tail -n +2"

    path_name = dirname(input.chromosome.size.file)
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(input.chromosome.size.file))

    file_name <- paste0(path_name2, "-", file_name, ".txt")

    cmd <- paste(cmd3, input.chromosome.size.file, ">", file.path(path_name,
        file_name), sep = " ")

    system(cmd)

    cmd4 <- "genomeCoverageBed -ibam"
    cmd5 <- "-bg -g"
    cmd6 <- ">"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam"))

        cmd <- paste(cmd4, x, cmd5, input.chromosome.size.file, cmd6, file.path(output.bw.file.dir,
            paste0(file_name, ".bdg")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    # re <- list(cmdl = cmd.l, output.dir = output.dir)

    # input.bdg.file.dir <- re$output.dir

    cmd7 <- "LC_COLLATE=C sort -k1,1 -k2,2n"

    # re <- parserreadfiles(input.bdg.file.dir,'bdg')

    # res <- re$input

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".bdg"))

        cmd <- paste(cmd7, x, cmd6, file.path(output.bw.file.dir, paste0(file_name,
            ".sorted_bdg")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

    cmd8 = "bedGraphToBigWig"

    input.chromosome.size.file.m <- file.path(path_name, file_name)

    # re <- parserreadfiles(input.bdg.file.dir,'.sorted_bdg')

    # res <- re$input

    cmd.l <- lapply(res, function(u, input.chromosome.size.file.m, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".sorted_bdg"))

        cmd <- paste(cmd8, x, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, ".bw")), sep = " ")

        system(cmd)

        cmd
    }, input.chromosome.size.file.m, output.bw.file.dir)

    cmd9 = "bigWigToWig"

    cmd.l <- lapply(res, function(u, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, ".bw"))

        cmd <- paste(cmd9, x, file.path(output.bw.file.dir, paste0(file_name,
            ".wig")), sep = " ")

        system(cmd)

        cmd
    }, output.bw.file.dir)

}

#' @examples
#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2StrandBw('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/nethome/axy148/hg19.genome','/scratch/projects/bbc/aiminy_project/DoGs/BW2')'
#'
convertBam2StrandBw <- function(input.bam.file.dir, input.chromosome.size.file,
    output.bw.file.dir)
    {
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    # job.name=paste0('bamSort[',length(res),']')

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamSort.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))

            # job.name=paste0('bamSort[',length(res),']') cmd1 = paste0('bsub -w
            # \'done(\'bamSort[*]\')\'', 'bsub -P bbc -J \'',job.name,paste0('\'
            # -o %J.log '),paste0('-e %J.err -W')) job.name=paste0('Bdg[',u,']') cmd1 =
            # paste0('bsub -w \'done(\'bamIndex[*]\') && done(\'Chrosome\')\'',
            # 'bsub -P bbc -J \'',job.name,paste0('\' -o %J.',job.name,'.log
            # '),paste0('-e %J.',job.name,'.err -W'))

            cmd2 = "samtools sort"
            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            cmd3 = "samtools sort"
        }

        cmd <- paste(cmd3, res[[u]], file.path(output.bw.file.dir, paste0(file_name,
            "_sorted")), sep = " ")

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamIndex.", u)
            wait.job.name = paste0("bamSort.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            # cmd1 = 'bsub -w \'done(\'bamSort\')\' -P bbc -J \'bamIndex\' -o
            # %J.bamIndex.log -e %J.bamIndex.err -W'
            cmd4 = "samtools index"
            cmd5 = paste(cmd1, cmd0, cmd4, sep = " ")
        } else
        {
            cmd5 = "samtools index"
        }

        cmd <- paste(cmd5, file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam")),
            sep = " ")

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

    path_name <- dirname(input.chromosome.size.file)
    path_name2 <- basename(path_name)
    file_name <- file_path_sans_ext(basename(input.chromosome.size.file))
    file_name <- paste0(path_name2, "-", file_name, ".txt")

    if (m.id == 1)
    {
        cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
        cmd1 = "bsub -P bbc -J \"Chrosome\" -o %J.Chrosome.log -e %J.Chrosome.err -W"
        cmd6 <- "tail -n +2"
        cmd7 <- paste(cmd1, cmd0, cmd6, sep = " ")
        cmd8 <- paste(cmd7, input.chromosome.size.file, "\\>", file.path(path_name,
            file_name), sep = " ")
    } else
    {
        cmd7 <- "tail -n +2"
        cmd8 <- paste(cmd7, input.chromosome.size.file, ">", file.path(path_name,
            file_name), sep = " ")
    }

    system(cmd8)

    # cmd11='bsub -w \'done(\'STAR-alignment\')\' -P bbc -J
    # \'samtools-sort\' -o %J.samtools-sort.log -e %J.samtools-sort.err -W'

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_sorted.bam"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

            # u=1
            job.name = paste0("Bdg.", u)
            wait.job.name = paste(paste0("bamIndex.", u, "\")"), "&& done(\"Chrosome\")",
                sep = " ")
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            cmd9 = "genomeCoverageBed -split -strand + -ibam"
            cmd10 = "genomeCoverageBed -split -strand - -ibam"
            cmd11 = paste(cmd1, cmd0, cmd9, sep = " ")
            cmd12 = paste(cmd1, cmd0, cmd10, sep = " ")
            cmd13 <- "\\>"
        } else
        {
            cmd11 = "genomeCoverageBed -split -strand + -ibam"
            cmd12 = "genomeCoverageBed -split -strand - -ibam"
            cmd13 <- ">"
        }

        cmd14 <- "-bg -g"

        cmd.x <- paste(cmd11, x, cmd14, input.chromosome.size.file, cmd13, file.path(output.bw.file.dir,
            paste0(file_name, "_plus.bdg")), sep = " ")

        cmd.y <- paste(cmd12, x, cmd14, input.chromosome.size.file, cmd13, file.path(output.bw.file.dir,
            paste0(file_name, "_minus.bdg")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd

    }, m.id, res, output.bw.file.dir)


    # re <- list(cmdl = cmd.l, output.dir = output.dir)

    # input.bdg.file.dir <- re$output.dir

    # re <- parserreadfiles(input.bdg.file.dir,'bdg')

    # res <- re$input

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.bdg"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.bdg"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("sortBdg.", u)
            wait.job.name = paste(paste0("Bdg.", u, "\")"), "&& done(\"Chrosome\")",
                sep = " ")
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))

            cmd15 <- "LC_COLLATE=C sort -k1,1 -k2,2n"
            cmd16 <- paste(cmd1, cmd0, cmd15, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd16 <- "LC_COLLATE=C sort -k1,1 -k2,2n"
            cmd17 <- ">"
        }

        cmd.x <- paste(cmd16, x, cmd17, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.sorted_bdg")), sep = " ")
        cmd.y <- paste(cmd16, y, cmd17, file.path(output.bw.file.dir, paste0(file_name,
            "_minus.sorted_bdg")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, output.bw.file.dir)


    input.chromosome.size.file.m <- file.path(path_name, file_name)

    # re <- parserreadfiles(input.bdg.file.dir,'.sorted_bdg')

    # res <- re$input

    cmd.l <- lapply(1:length(res), function(u, m.id, res, input.chromosome.size.file.m,
        output.bw.file.dir)
        {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.sorted_bdg"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.sorted_bdg"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

            job.name = paste0("BigWig.", u)
            wait.job.name = paste0("sortBdg.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))
            cmd18 = "bedGraphToBigWig"
            cmd19 <- paste(cmd1, cmd0, cmd18, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd19 <- "bedGraphToBigWig"
        }

        cmd.x <- paste(cmd19, x, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, "_plus.bw")), sep = " ")

        cmd.y <- paste(cmd19, y, input.chromosome.size.file.m, file.path(output.bw.file.dir,
            paste0(file_name, "_minus.bw")), sep = " ")


        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, input.chromosome.size.file.m, output.bw.file.dir)


    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(u) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)

        x = file.path(output.bw.file.dir, paste0(file_name, "_plus.bw"))
        y = file.path(output.bw.file.dir, paste0(file_name, "_minus.bw"))

        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"


            job.name = paste0("Wig.", u)
            wait.job.name = paste0("BigWig.", u)
            cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                job.name, paste0("\" -o %J.", job.name, ".log "), paste0("-e %J.",
                  job.name, ".err -W"))
            cmd20 = "bigWigToWig"
            cmd21 <- paste(cmd1, cmd0, cmd20, sep = " ")
            cmd17 <- "\\>"
        } else
        {
            cmd21 <- "bigWigToWig"
        }

        cmd.x <- paste(cmd21, x, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.wig")), sep = " ")

        cmd.y <- paste(cmd21, y, file.path(output.bw.file.dir, paste0(file_name,
            "_plus.wig")), sep = " ")

        system(cmd.x)
        system(cmd.y)

        cmd <- list(cmd.x = cmd.x, cmd.y = cmd.y)
        cmd
    }, m.id, res, output.bw.file.dir)

}

# HCI-Scripts/BamFile/split_bam_by_isize.pl
#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::splitBam('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs/Bam_split')'
#'
splitBam <- function(input.bam.file.dir, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    # job.name=paste0('bamSort[',length(res),']')

    cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bw.file.dir)
    {
        # path_name = dirname(res[[u]]) path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(res[[u]]))

        # file_name <- paste0(path_name2,'-',file_name)
        u <- 3
        if (m.id == 1)
        {
            cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            job.name = paste0("bamSplit.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))

            # job.name=paste0('bamSort[',length(res),']') cmd1 = paste0('bsub -w
            # \'done(\'bamSort[*]\')\'', 'bsub -P bbc -J \'',job.name,paste0('\'
            # -o %J.log '),paste0('-e %J.err -W')) job.name=paste0('Bdg[',u,']') cmd1 =
            # paste0('bsub -w \'done(\'bamIndex[*]\') && done(\'Chrosome\')\'',
            # 'bsub -P bbc -J \'',job.name,paste0('\' -o %J.',job.name,'.log
            # '),paste0('-e %J.',job.name,'.err -W'))
            if (u <= 6)
            {
                cmd2 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd2 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            }
            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            if (u <= 6)
            {
                cmd3 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd3 = paste("$HOME/HCI-Scripts/BamFile/split_bam_by_isize.pl --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_split")),
                  "--in", res[[u]], sep = " ")
            }
        }

        cmd <- cmd3

        system(cmd)

        cmd
    }, m.id, res, output.bw.file.dir)

}

#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2StrandBw2('/scratch/projects/bbc/aiminy_project/DoGs/Bam_split/','/scratch/projects/bbc/aiminy_project/DoGs/BW_split')'
#'

convertBam2StrandBw2 <- function(input.bam.file.dir, output.bw.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
        span.ptile, res, output.bw.file.dir)
        {

        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            job.name = paste0("bam2wig.", u)
            cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 16000, span.ptile = 16, job.name)

            if (u <= 6)
            {
                cmd2 = paste("bam2wig.pl -pe --pos span --strand --bw  --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd2 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            }
            cmd3 = paste(cmd1, cmd2, sep = " ")
        } else
        {
            if (u <= 6)
            {
                cmd3 = paste("bam2wig.pl -pe --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            } else
            {
                cmd3 = paste("bam2wig.pl --pos span --strand --bw --bwapp $HOME/kentUtils/bin/wigToBigWig --out",
                  file.path(output.bw.file.dir, paste0(file_name, "_2.bw")),
                  "--in", res[[u]], sep = " ")
            }
        }

        cmd <- cmd3

        system(cmd)

        cmd
    }, m.id, Wall.time, cores, Memory, span.ptile, res, output.bw.file.dir)

}

#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::convertBam2bed2('/scratch/projects/bbc/aiminy_project/DoGs/BAM','/scratch/projects/bbc/aiminy_project/DoGs/BED2')'
#'

convertBam2bed2 <- function(input.bam.file.dir, output.bed.file.dir)
{
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bed.file.dir))
    {
        dir.create(output.bed.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
        span.ptile, res, output.bed.file.dir)
        {

        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            job.name = paste0("bam2bed.", u)
            cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 16000, span.ptile = 16, job.name)

            cmd2 = paste("bam2gff_bed.pl -pe --bed --out", file.path(output.bed.file.dir,
                paste0(file_name, ".bed")), "--in", res[[u]], sep = " ")
            cmd3 = paste(cmd1, cmd2, sep = " ")
        } else
        {
            cmd3 = paste("bam2gff_bed.pl -pe --bed --out", file.path(output.bed.file.dir,
                paste0(file_name, ".bed")), "--in", res[[u]], sep = " ")
        }

        cmd <- cmd3

        system(cmd)

        cmd
    }, m.id, Wall.time, cores, Memory, span.ptile, res, output.bed.file.dir)

}

#'bedtools intersect -v -a 'Results/''Aligned'.bed -b /media/H_driver/2016/Ramin_azhang/Annotation/exons.bed /media/H_driver/2016/Ramin_azhang/Annotation/intron.bed

#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::removeReadsOnExonIntron("/scratch/projects/bbc/aiminy_project/DoGs/BED2","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/BedRmExonIntron")'
#'
#'
#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::removeReadsOnExonIntron("/scratch/projects/bbc/aiminy_project/DoGs/BedFileFromBam","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/BedRmExonIntron")'
#'
removeReadsOnExonIntron <- function(input.bed.file.dir, annotation.bed.file.dir,
    output.bed.file.dir)
    {
    re <- parserreadfiles(input.bed.file.dir, "bed")

    res <- re$input

    annotationBed <- parserreadfiles(annotation.bed.file.dir,"bed",sample.group=c("hg19_exons.bed","hg19_intron.bed"))

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bed.file.dir))
    {
        dir.create(output.bed.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
        span.ptile, res, annotationBed, output.bed.file.dir)
        {

        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            job.name = paste0("bedRmExonIntron.", u)
            cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                Memory = 16000, span.ptile = 16, job.name)
            exon.intron <- paste(unlist(annotationBed$input),collapse=" ")
            cmd2 = paste("bedtools intersect -v -a", res[[u]], "-b", exon.intron,
                "\\>", file.path(output.bed.file.dir, paste0(file_name, "_rm_exon_intron.bed")),
                sep = " ")
            cmd3 = paste(cmd1, cmd2, sep = " ")
        } else
        {
            cmd3 = paste("bedtools intersect -v -a", res[[u]], "-b", exon, intron,
                "\\>", file.path(output.bed.file.dir, paste0(file_name, "_rm_exon_intron.bed")),
                sep = " ")
        }

        cmd <- cmd3

        cat(cmd, "\n\n")

        system(cmd)

        cmd
    }, m.id, Wall.time, cores, Memory, span.ptile, res, annotationBed, output.bed.file.dir)

}

#'bedtools window -a /media/H_driver/2016/Ramin_azhang/Annotation/hg19_gene.bed -b "/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/""2016-02-10-emp1_WT".rm.exon.intron.hg19.bed -l 0 -r 4500 -sw   >

#'R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::getCount4Downstream(""/scratch/projects/bbc/aiminy_project/DoGs/BedRmExonIntron","/projects/ctsi/bbc/aimin/annotation/","/scratch/projects/bbc/aiminy_project/DoGs/Counts45KB")'
#'
getCount4Downstream <- function(input.bed.file.dir, annotation.bed.file.dir,
                                      output.count.file.dir)
  {
    re <- parserreadfiles(input.bed.file.dir, "bed")

    res <- re$input

    annotationBed <- parserreadfiles(annotation.bed.file.dir,"bed",sample.group=c("hg19_gene.bed"))

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.count.file.dir))
    {
      dir.create(output.count.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
                                            span.ptile, res, annotationBed, output.count.file.dir)
    {

      file_name = file_path_sans_ext(basename(res[[u]]))

      if (m.id == 1)
      {
        job.name = paste0("bed2count.", u)
        cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
                                     Memory = 16000, span.ptile = 16, job.name)
        exon.intron <- paste(unlist(annotationBed$input),collapse=" ")
        cmd2 = paste("bedtools window -a",exon.intron,"-b",res[[u]],"-l 0 -r 45000 -sw -c",
                     "\\>", file.path(output.count.file.dir, paste0(file_name, "_downstream_count.txt")),
                     sep = " ")
        cmd3 = paste(cmd1, cmd2, sep = " ")
      } else
      {
        cmd3 = paste("bedtools window -a",exon.intron,"-b",res[[u]],"-l 0 -r 45000 -sw -c",
                             ">", file.path(output.count.file.dir, paste0(file_name, "_downstream_count.txt")),
                             sep = " ")
      }

      cmd <- cmd3

      cat(cmd, "\n\n")

      system(cmd)

      cmd
    }, m.id, Wall.time, cores, Memory, span.ptile, res, annotationBed, output.count.file.dir)

  }

#' res <- convertCountFile2Table("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt")
#'
convertCountFile2Table <- function(input.count.file.dir,input.file.pattern)
{
  file.name = file.path(input.count.file.dir, dir(input.count.file.dir, recursive = TRUE, pattern = input.file.pattern))

  res <- lapply(file.name, function(u)
  {
    if (!file.size(u) == 0)
    {
      re <- read.table(u,header=FALSE)

      file_name = basename(u)

      re1 <- as.data.frame(re[,c(4,dim(re)[2])],stringsAsFactors=FALSE)

      colnames(re1) <- c("Gene",file_name)

    }
    re1
  })

  merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
  }

  res1 <- Reduce(merge.all,res)

  if(length(which(is.na(res1))) > 0){
    res1[is.na(res1)] <- 0
  }

  res2 <- res1[,-1]
  row.names(res2) <- res1$Gene

  return(res2)

}

#' res <- convertCountFile2Table("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt")

#' res.new <- ThreeUTR:::matchAndDE(res,file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),group.comparision = c("condition","Untreated","Treated"))
#'
#' res.new <- ThreeUTR:::matchAndDE(res,file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),group.comparision = c("condition","Treated","Untreated"))
#'
matchAndDE <- function(res,sample.infor.file,group.comparision = c("condition","Untreated","Treated"))
{
  file.name <- colnames(res)

  con.name <- group.comparision[1]

  sample.infor <- read.table(sample.infor.file, header = TRUE)
  sample.name <- unique(as.character(sample.infor[,1]))

  file.name.new  <- lapply(file.name, function(u,sample.name)
  {
      x <- u
      y <- sample.name

      xx <- lapply(y,function(u,x){

        if(length(grep(u,x))>0)
        {
        z <- u
        z
        }

      },x)

      xx
  },sample.name)
  colnames(res) <- unlist(file.name.new)

  group <- unique(as.character(sample.infor[,2]))

  a.name <- group[which(group %in% group.comparision[2])]
  b.name <- group[which(group %in% group.comparision[3])]

  group.1.index <- sample.infor[which(sample.infor[,2] %in% a.name),1]

  group.2.index <- sample.infor[which(sample.infor[,2] %in% b.name),1]

  countData <- res[,c(which(colnames(res) %in% group.1.index),which(colnames(res) %in% group.2.index))]

   a <- length(group.1.index)
   b <- length(group.2.index)

   colData <- data.frame(row.names=colnames(countData),x = factor(c(rep(a.name,a),rep(b.name,b))))

   colnames(colData)[1] <-  con.name

   dds <- DESeqDataSetFromMatrix(countData, colData, formula(paste("~",con.name)))

  x <- colnames(colData)[which(colnames(colData) %in% group.comparision[1])]
  c1 <- a.name
  c2 <- b.name
  contrast.set <- c(x,c1,c2)

  dds <- DESeq(dds)

  re.DESeq <- results(dds,contrast = contrast.set)

  re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
  colnames(re.FC)[7] = paste0("FoldChange(",a.name,"-vs-",b.name,")")

  txs.gene <- ReformatTxsGene()

  re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

  re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

  rr <- list(countData = countData,colData=colData,dds=dds, re.DESeq=re.DESeq,re.FC=re.FC,re.FC.sorted=re.FC.sorted)
  return(rr)
}


#' library(org.Hs.eg.db)
#' res12 <- ThreeUTR:::sumCount4Downstream("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),c(1,2))

#'res21 <- ThreeUTR:::sumCount4Downstream("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt",file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),c(2,1))
#'
sumCount4Downstream <- function(input.count.file.dir,input.file.pattern,sample.infor.file,group.comparision = c("condition",1,2))
{
  file.name = file.path(input.count.file.dir, dir(input.count.file.dir, recursive = TRUE, pattern = input.file.pattern))
  #file.name.2 <- as.list(file.name)

  #names(file.name.2) = basename(file.name)

  res <- lapply(file.name, function(u)
  {
    if (!file.size(u) == 0)
    {
      re <- read.table(u,header=FALSE)

      file_name = basename(u)
      pos = gregexpr("-", file_name)
      a = pos[[1]][1] + 1
      b = pos[[1]][2] - 1
      sampleName = substr(file_name, a, b)

      re1 <- as.data.frame(re[,c(4,dim(re)[2])],stringsAsFactors=FALSE)

      colnames(re1) <- c("gene",sampleName)

    }
    re1
  })

  merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="gene")
  }

  res1 <- Reduce(merge.all,res)

  if(length(which(is.na(res1))) > 0){
    res1[is.na(res1)] <- 0
  }

  res2 <- res1[,-1]
  row.names(res2) <- res1$gene


  #colnames(res)[dim(res)[2]]="Sample"

  # filterByRmNull <- function(a.list)
  # {
  #   a.list.2 <- a.list[lapply(a.list, length) >0]
  #   return(a.list.2)
  # }
  #
  # file.name.4 <- filterByRmNull(file.name.3)
  #
   sample.infor <- read.table(sample.infor.file, header = TRUE)

   group <- unique(as.character(sample.infor$Condition))

   group.1.index <- sample.infor[which(sample.infor$Condition %in% group[group.comparision[2]]),]$Sample

   group.2.index <- sample.infor[which(sample.infor$Condition %in% group[group.comparision[3]]),]$Sample

   res3 <- res2[,c(which(colnames(res2) %in% group.1.index),which(colnames(res2) %in% group.2.index))]

   countData <- res3

   a <- length(group.1.index)
   b <- length(group.2.index)

   colData <- data.frame(row.names=colnames(countData),condition = factor(c(rep(group[group.comparision[2]],a),rep(group[group.comparision[3]],b))))

    dds <- DESeqDataSetFromMatrix(countData, colData, formula(~condition))

    x <- colnames(colData)[which(colnames(colData) %in% group.comparision[1])]
    c1 <- group[group.comparision[2]]
    c2 <- group[group.comparision[3]]

    contrast.set <- c(x,c1,c2)

    re.DESeq <- results(DESeq(dds),contrast = contrast.set)

    re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[, 2], counts(dds))
    colnames(re.FC)[7] = paste0("FoldChange(",group[group.comparision[2]],"-vs-",group[group.comparision[3]],")")

    txs.gene <- ReformatTxsGene()

    re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

    re.FC.sorted <- re.FC[order(re.FC$pvalue), ]



   rr <- list(res1=res1,res2=res2,res3=res3,colData=colData,re.DESeq=re.DESeq,re.FC=re.FC,re.FC.sorted=re.FC.sorted)


   #res <- merge(res,sample.infor,by="Sample")

   #res <- cbind.data.frame(res,paste(res[,c(paste0("V",1):paste0("V",12))],collapse="-"),stringsAsFactors=FALSE)



  # names(file.name.4) = unlist(lapply(1:length(file.name.4), function(u, file.name.4,
  #                                                                    sample.infor)
  # {
  #
  #   tmp = file.name.4
  #   x = tmp[[u]]
  #   path_name = dirname(x)
  #   file_name = basename(x)
  #   pos = gregexpr("-", file_name)
  #   a = pos[[1]][1] + 1
  #   b = pos[[1]][2] - 1
  #   c = substr(file_name, a, b)
  #   con = sample.infor[match(c, sample.infor$Sample), ]$Condition
  #   file_name_con = paste0(con, "-", file_name)
  #
  #   file_name_con
  #
  #
  # }, file.name.4, sample.infor))

  return(rr)

}

# volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
#   with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
#   with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
#   with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
#   with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
#   if (labelsig) {
#     require(calibrate)
#     with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
#   }
#   legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
# }
# png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
# volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
# dev.off()


prepareDaPars <- function(input.wig.file.dir, sample.group = c("Dox", "WT"),
    output.dir)
    {
    re <- parserreadfiles(input.wig.file.dir, "wig")

    res <- re$input

    cmd = "tail -n +2"

    cmd.l <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)
        # path_name2 <- basename(path_name)

        file_name = file_path_sans_ext(basename(u))

        # file_name <- paste0(path_name2,'-',file_name)

        x = u

        cmd <- paste(cmd, x, ">", file.path(output.dir, paste0(file_name, "2.wig")),
            sep = " ")

        system(cmd)

        cmd
    }, output.dir)

    g1 <- "1,2,3,4,5"
    g2 <- "a,b,c,d,e"

    cat(c("Annotated_3UTR=hg19_refseq_extracted_3UTR.bed\n", paste0("Group1_Tophat_aligned_Wig=",
        g1, "\n"), paste0("Group2_Tophat_aligned_Wig=", g2, "\n"), "Output_directory=DaPars_Dox_WT_data/\n",
        "Output_result_file=DaPars_Dox_WT_data\n", "#Parameters\n
        Num_least_in_group1=1\n
        Num_least_in_group2=1\n
        Coverage_cutoff=30\n
        FDR_cutoff=0.05\n
        PDUI_cutoff=0.5\n
        Fold_change_cutoff=0.59\n"),
        file = file.path(output.dir, "config_test.txt"), sep = "")
}

select3UTR <- function(genome, tablename)
{
    library(GenomicFeatures)
    library(dplyr)
    refSeq <- makeTxDbFromUCSC(genome = "hg19", tablename = "knownGene")
    threeUTRs <- threeUTRsByTranscript(refSeq, use.names = TRUE)
    length_threeUTRs <- width(ranges(threeUTRs))
    the_lengths <- as.data.frame(length_threeUTRs)
    the_lengths <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
    the_lengths <- unique(the_lengths[, c("group_name", "sum(value)")])
    colnames(the_lengths) <- c("RefSeq Transcript", "3' UTR Length")

}

subsetFastq <- function(input.fastq.files.dir, output.dir, n)
{
    re <- parserreadfiles(input.fastq.files.dir, "fastq")
    res <- re$input

    # cmd0 = paste0('head -',n)
    cmd0 = "seqtk sample -s100"
    # read1.fq 10000

    cmd1 = "\\>"
    # cmd1 = '>'

    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd3 = "bsub -P bbc -J \"subSetFastq\" -o %J.subSetFastq.log -e %J.subSetFastq.err -W"

    # cmd4='bsub -w \'done(\'subSetFastq\')\' -P bbc -J \'tophat\' -o
    # %J.tophat.log -e %J.tophat.err -W'

    xx <- lapply(res, function(u, output.dir)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        sample.name.out = file.path(output.dir, paste0(file_name, "-test-",
            n, ".fastq"))

        cmd = paste(cmd3, cmd2, cmd0, u, n, cmd1, sample.name.out, sep = " ")

        print(cmd)

        cmd

        system(cmd, intern = TRUE)

    }, output.dir)

}

checkStrand <- function(input.alignment.dir)
{
    re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
        pattern = "junctions.bed"))

    cmd0 <- "wc -l"

    y <- lapply(re, function(u)
    {
        cmd <- paste(cmd0, u, sep = " ")

        # system(cmd,intern = TRUE)
        cmd
    })

    yyy <- lapply(y, function(u)
    {
        yy <- system(u, intern = TRUE)
        yy
    })

    yyyy <- unlist(yyy)

    return(yyyy)
}

#' @examples
#' R -e 'library(ThreeUTR);ThreeUTR:::processBamFiles('/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2','/scratch/projects/bbc/aiminy_project/DoGs/BAM')'
#'
processBamFiles <- function(input.alignment.dir, output.dir)
{
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
        pattern = "accepted_hits.bam"))

    cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd1 = "bsub -P bbc -J \"bamProcess\" -o %J.bamProcess.log -e %J.bamProcess.err -W"

    cmd2 <- "mv"

    y <- lapply(re, function(u)
    {
        # /SRR2038198/Fs12/accepted_hits.bam

        dir.name.0 <- dirname(u)
        dir.name.1 <- dirname(dir.name.0)

        x <- basename(dir.name.0)
        y <- basename(dir.name.1)
        file.name <- basename(u)

        sample.name <- paste(y, x, file.name, sep = "-")

        cmd <- paste(cmd1, cmd0, cmd2, u, file.path(output.dir, sample.name),
            sep = " ")

        # system(cmd,intern = TRUE
        cmd
    })

    print(y)

    yyy <- lapply(y, function(u)
    {
        yy <- system(u, intern = TRUE)
        yy
    })

    yyyy <- unlist(yyy)

    return(yyyy)
}

testAlignment <- function(output.dir, gene.model.file, genome.index)
{
    cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

    cmd4 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W"

    re <- parserreadfiles(output.dir, "fastq")
    res <- re$input

    xx <- lapply(res, function(u)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        if (regexpr(pattern = "_", file_name) != -1)
        {
            # cat('match:',file_name,'\n')

            p <- regexpr(pattern = "_", file_name)
            pp <- p - 1
            x <- substr(file_name, 1, pp)

        } else
        {
            # cat('no match:',file_name,'\n')
            x <- file_name
        }

        x
    })


    xxx <- unique(unlist(xx))
    res2 <- unlist(res)

    cmd5 = "tophat --library-type fr-unstranded -g 1 -G"
    cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
    cmd7 = "tophat --library-type fr-secondstrand -g 1 -G"
    cmd8 = "-p 4 -o"
    cmd9 = "mv"

    for (i in 1:length(xxx))
    {
        sample.name <- xxx[i]
        sample.name.out.dir <- file.path(output.dir, sample.name)

        if (!dir.exists(sample.name.out.dir))
        {
            dir.create(sample.name.out.dir, recursive = TRUE)
        }

        y <- res2[grep(xxx[i], res2)]

        # print(y) print(length(y))

        if (length(y) == 2)
        {
            yy1 <- basename(y[1])
            yy2 <- basename(y[2])

            p1 <- regexpr(pattern = "_", yy1)
            pp1 <- p1 + 1
            x1 <- substr(yy1, pp1, pp1)

            p2 <- regexpr(pattern = "_", yy2)
            pp2 <- p2 + 1
            x2 <- substr(yy2, pp2, pp2)


            sample.name.out.dir.1 = file.path(sample.name.out.dir, paste0("Us",
                x1, x2))
            sample.name.out.dir.2 = file.path(sample.name.out.dir, paste0("Us",
                x2, x1))
            sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                x1, x2))
            sample.name.out.dir.4 = file.path(sample.name.out.dir, paste0("Fs",
                x2, x1))
            sample.name.out.dir.5 = file.path(sample.name.out.dir, paste0("Ss",
                x1, x2))
            sample.name.out.dir.6 = file.path(sample.name.out.dir, paste0("Ss",
                x2, x1))

            if (!dir.exists(sample.name.out.dir.1))
            {
                dir.create(sample.name.out.dir.1, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.2))
            {
                dir.create(sample.name.out.dir.2, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.3))
            {
                dir.create(sample.name.out.dir.3, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.4))
            {
                dir.create(sample.name.out.dir.4, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.5))
            {
                dir.create(sample.name.out.dir.5, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.6))
            {
                dir.create(sample.name.out.dir.6, recursive = TRUE)
            }

            cmd10 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.1,
                genome.index, y[1], y[2], sep = " ")
            cmd11 = paste(cmd4, cmd2, cmd10)
            print(cmd10)


            system(cmd11)

            cmd12 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.2,
                genome.index, y[2], y[1], sep = " ")
            cmd13 = paste(cmd4, cmd2, cmd12)
            system(cmd13)

            cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                genome.index, y[1], y[2], sep = " ")
            cmd15 = paste(cmd4, cmd2, cmd14)
            system(cmd15)

            cmd16 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.4,
                genome.index, y[2], y[1], sep = " ")
            cmd17 = paste(cmd4, cmd2, cmd16)
            system(cmd17)

            cmd18 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.5,
                genome.index, y[1], y[2], sep = " ")
            cmd19 = paste(cmd4, cmd2, cmd18)
            system(cmd19)

            cmd20 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.6,
                genome.index, y[2], y[1], sep = " ")
            cmd21 = paste(cmd4, cmd2, cmd20)
            system(cmd21)
        } else
        {
            sample.name.out.dir.7 = file.path(sample.name.out.dir, "Us")
            sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")
            sample.name.out.dir.9 = file.path(sample.name.out.dir, "Ss")

            if (!dir.exists(sample.name.out.dir.7))
            {
                dir.create(sample.name.out.dir.7, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.8))
            {
                dir.create(sample.name.out.dir.8, recursive = TRUE)
            }

            if (!dir.exists(sample.name.out.dir.9))
            {
                dir.create(sample.name.out.dir.9, recursive = TRUE)
            }

            cmd22 = paste(cmd5, gene.model.file, cmd8, sample.name.out.dir.7,
                genome.index, y[1], sep = " ")
            cmd23 = paste(cmd4, cmd2, cmd22)
            # print(cmd23)

            system(cmd23)

            cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                genome.index, y[1], sep = " ")
            cmd25 = paste(cmd4, cmd2, cmd24)
            system(cmd25)

            cmd26 = paste(cmd7, gene.model.file, cmd8, sample.name.out.dir.9,
                genome.index, y[1], sep = " ")
            cmd27 = paste(cmd4, cmd2, cmd26)
            system(cmd27)
        }

    }
}

#'ThreeUTR:::useTophat4Alignment("/scratch/projects/bbc/aiminy_project/DoGsFastq","/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","General")
#'
useTophat4Alignment <- function(input.fastq.files.dir, output.dir, gene.model.file = NULL,
    genome.index, cmd.input)
    {
    if (!dir.exists(output.dir))
    {
        dir.create(output.dir, recursive = TRUE)
    }

    if (cmd.input == "General")
    {
        cmd2 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
        cmd4 = "bsub -P bbc -J \"tophat\" -o %J.tophat.log -e %J.tophat.err -W"
    } else
    {

    }

    re <- parserreadfiles(input.fastq.files.dir, "fastq")
    res <- re$input

    xx <- lapply(res, function(u)
    {
        path_name = dirname(u)

        file_name = file_path_sans_ext(basename(u))

        if (regexpr(pattern = "_", file_name) != -1)
        {
            # cat('match:',file_name,'\n')

            p <- regexpr(pattern = "_", file_name)
            pp <- p - 1
            x <- substr(file_name, 1, pp)

        } else
        {
            # cat('no match:',file_name,'\n')
            x <- file_name
        }

        x
    })

    xxx <- unique(unlist(xx))
    res2 <- unlist(res)

    cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
    cmd8 = "-p 4 -o"
    cmd9 = "mv"

    for (i in 1:length(xxx))
    {
        sample.name <- xxx[i]
        sample.name.out.dir <- file.path(output.dir, sample.name)

        if (!dir.exists(sample.name.out.dir))
        {
            dir.create(sample.name.out.dir, recursive = TRUE)
        }

        y <- res2[grep(xxx[i], res2)]

        # print(y) print(length(y))

        if (length(y) == 2)
        {
            yy1 <- basename(y[1])
            yy2 <- basename(y[2])

            p1 <- regexpr(pattern = "_", yy1)
            pp1 <- p1 + 1
            x1 <- substr(yy1, pp1, pp1)

            p2 <- regexpr(pattern = "_", yy2)
            pp2 <- p2 + 1
            x2 <- substr(yy2, pp2, pp2)


            sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                x1, x2))


            if (!dir.exists(sample.name.out.dir.3))
            {
                dir.create(sample.name.out.dir.3, recursive = TRUE)
            }

            cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                genome.index, y[1], y[2], sep = " ")
            # cmd15= paste(cmd.input,cmd14)

            cmd15 = paste(cmd4, cmd2, cmd14)

            system(cmd15)

        } else
        {
            sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")

            if (!dir.exists(sample.name.out.dir.8))
            {
                dir.create(sample.name.out.dir.8, recursive = TRUE)
            }
            cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                genome.index, y[1], sep = " ")
            # cmd25= paste(cmd.input,cmd24)
            cmd25 = paste(cmd4, cmd2, cmd24)
            system(cmd25)
        }

    }

}

#' @examples
#' R -e 'library(ThreeUTR);ThreeUTR:::subsetBam('/scratch/projects/bbc/aiminy_project/DoGs/Bam_split/','$HOME/1833_common_gene.bed','/scratch/projects/bbc/aiminy_project/DoGs/BAMSubSet',BigMem=FALSE)'
#'
subsetBam <- function(input.bam.file.dir, region.bed.file, output.bw.file.dir,
    BigMem = FALSE)
    {
    re <- parserreadfiles(input.bam.file.dir, "bam")

    res <- re$input

    m.id <- grep("login", system("hostname", intern = TRUE))

    if (!dir.exists(output.bw.file.dir))
    {
        dir.create(output.bw.file.dir, recursive = TRUE)
    }

    cmd.l <- lapply(1:length(res), function(u, m.id, res, region.bed, BigMem,
        output.bw.file.dir)
        {
        file_name = file_path_sans_ext(basename(res[[u]]))

        if (m.id == 1)
        {
            if (BigMem == TRUE)
            {
                cmd0 = "72:00 -n 16 -q bigmem -R 'rusage[mem=36864] span[ptile=8]' -u aimin.yan@med.miami.edu"
            } else
            {
                cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"
            }

            job.name = paste0("subsetBam.", u)
            cmd1 = paste0("bsub -P bbc -J \"", job.name, paste0("\" -o %J.",
                job.name, ".log "), paste0("-e %J.", job.name, ".err -W"))
            cmd2 = paste("samtools view -b -L", region.bed.file, res[[u]], "\\>",
                file.path(output.bw.file.dir, paste0(file_name, "_subset.bam")),
                sep = " ")

            cmd3 = paste(cmd1, cmd0, cmd2, sep = " ")
        } else
        {
            cmd3 = paste("samtools view -b -L", region.bed.file, res[[u]], ">",
                file.path(output.bw.file.dir, paste0(file_name, "_subset.bam")),
                sep = " ")
        }

        cmd <- cmd3
        system(cmd)
        cmd
    }, m.id, res, region.bed, BigMem, output.bw.file.dir)

}
#' Title
#'
#' @param z
#' @param one.sided
#'
#' @return
#' @export
#'
#' @examples
#' convert.z.score(1.454648e-04,two.sided=T)
#' convert.z.score(1.454648e-04,two.sided=F)
#'
#' convert.z.score(-1.428500e-04,two.sided=T)
#' convert.z.score(-1.428500e-04,two.sided=F)
#' lapply(re.3UTR$re.FC.sorted$pvalue,function(x){
#'  y<-convert.z.score(x,two.side=F)
#'  y
#'  }
#'  )
#'
#' one.sided.p<-convert.z.score(re.3UTR$re.FC.sorted$pvalue,two.sided=F)
#'
#'
convert.z.score <- function(z, two.sided = T)
{
    if ((two.sided))
    {
        pval = 2 * pnorm(-abs(z))
    } else if (sign(z) == -1)
    {
        pval = pnorm(-abs(z))
    } else
    {
        pval = 1 - pnorm(-abs(z))
    }
    return(pval)
}
sanitizeColData <- function(object)
{
    if (is.null(mcols(colData(object))))
    {
        mcols(colData(object)) <- DataFrame(type = rep("input", ncol(colData(object))),
            description = character(ncol(colData(object))))
    }
    class(mcols(colData(object))$type) <- "character"
    class(mcols(colData(object))$description) <- "character"
    mcols(colData(object))$type[is.na(mcols(colData(object))$type)] <- ""
    mcols(colData(object))$description[is.na(mcols(colData(object))$description)] <- ""
    object
}

estimateSizeFactors.DESeqDataSet <- function(object, type = c("ratio", "iterate"),
    locfunc = stats::median, geoMeans, controlGenes, normMatrix)
    {
    type <- match.arg(type, c("ratio", "iterate"))
    # Temporary hack for backward compatibility with 'old' DESeqDataSet objects.
    # Remove once all serialized DESeqDataSet objects around have been updated.
    if (!.hasSlot(object, "rowRanges"))
        object <- updateObject(object)
    object <- sanitizeColData(object)
    if (type == "iterate")
    {
        sizeFactors(object) <- estimateSizeFactorsIterate(object)
    } else
    {
        if ("avgTxLength" %in% assayNames(object))
        {
            nm <- assays(object)[["avgTxLength"]]
            nm <- nm/exp(rowMeans(log(nm)))  # divide out the geometric mean
            normalizationFactors(object) <- estimateNormFactors(counts(object),
                normMatrix = nm, locfunc = locfunc, geoMeans = geoMeans, controlGenes = controlGenes)
            message("using 'avgTxLength' from assays(dds), correcting for library size")
        } else if (missing(normMatrix))
        {
            sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object),
                locfunc = locfunc, geoMeans = geoMeans, controlGenes = controlGenes)
        } else
        {
            normalizationFactors(object) <- estimateNormFactors(counts(object),
                normMatrix = normMatrix, locfunc = locfunc, geoMeans = geoMeans,
                controlGenes = controlGenes)
            message("using 'normMatrix', adding normalization factors which correct for library size")
        }
    }
    object
}

filterByRmNull <- function(a.list)
{  a.list.2 <- a.list[lapply(a.list, length) >0]
   return(a.list.2)
}

# Rfun <- 'library(ChipSeq);library(ThreeUTR);re <- ThreeUTR:::useFastqDumpConvertSra2Fastq("/nethome/axy148/DoGsExample","/scratch/projects/bbc/aiminy_project/DoGsFastq",wait.job.name = "wgetDownload")'

createBubRfun <- function(Rfun,job.name,wait.job.name){
  x <- ChipSeq:::usePegasus("parallel","72:00",16,25000,8,job.name,wait.job.name)
  xx <- paste(x,paste0("\"R -e ",paste0("\'",Rfun,"\'"),"\""),sep=" ")
  xx
}
#test <- createBubRfun(Rfun,"sra2fastq[1-8]","wgetDownload")
#system(test)





#' R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::useTophat4Alignment2("/scratch/projects/bbc/aiminy_project/DoGsFastq","/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf","/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome","General","sra2fastq")'

#' On pegasus
#'
#'bsub -P bbc -J "alignment[1-8]" -o %J.alignment.%I -e %J.alignment.%I -W 72:00 -n 16 -q parallel -R 'rusage[mem= 25000 ] span[ptile= 8 ]' -u aimin.yan@med.miami.edu "R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::useTophat4Alignment2(\"/scratch/projects/bbc/aiminy_project/DoGsFastq\",\"/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome\",\"General\")'"
#'
#'bsub -P bbc -J "alignment" -o %J.alignment.log -e %J.alignment.err -W 72:00 -n 16 -q parallel -R 'rusage[mem= 25000 ] span[ptile= 8 ]' -u aimin.yan@med.miami.edu "R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::useTophat4Alignment2(\"/scratch/projects/bbc/aiminy_project/DoGsFastq\",\"/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome\")'"
#'
#'
#'
#'
#'
#'
useTophat4Alignment2 <- function(input.fastq.files.dir, output.dir, gene.model.file = NULL,
                                genome.index, cmd.input="parallel",wait.job.name=NULL)
{

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  re <- parserreadfiles(input.fastq.files.dir, "fastq")
  res <- re$input

  xx <- lapply(res, function(u)
  {
    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    if (regexpr(pattern = "_", file_name) != -1)
    {
      # cat('match:',file_name,'\n')

      p <- regexpr(pattern = "_", file_name)
      pp <- p - 1
      x <- substr(file_name, 1, pp)

    } else
    {
      # cat('no match:',file_name,'\n')
      x <- file_name
    }

    x
  })

  xxx <- unique(unlist(xx))
  res2 <- unlist(res)

  m.id <- grep("login", system("hostname", intern = TRUE))
  cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
  cmd8 = "-p 4 -o"

  for (i in 1:length(xxx))
  {
    sample.name <- xxx[i]
    sample.name.out.dir <- file.path(output.dir, sample.name)

    if (!dir.exists(sample.name.out.dir))
    {
      dir.create(sample.name.out.dir, recursive = TRUE)
    }

    y <- res2[grep(xxx[i], res2)]

    if (length(y) == 2)
    {
      yy1 <- basename(y[1])
      yy2 <- basename(y[2])

      p1 <- regexpr(pattern = "_", yy1)
      pp1 <- p1 + 1
      x1 <- substr(yy1, pp1, pp1)

      p2 <- regexpr(pattern = "_", yy2)
      pp2 <- p2 + 1
      x2 <- substr(yy2, pp2, pp2)


      sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                                                                    x1, x2))

      if (!dir.exists(sample.name.out.dir.3))
      {
        dir.create(sample.name.out.dir.3, recursive = TRUE)
      }

      cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                    genome.index, y[1], y[2], sep = " ")
      #if(length(m.id)==1){

      job.name <- paste0("Alignment.",i)
      if(!is.null(wait.job.name))
      {
      x <- paste0(wait.job.name,".",i)
      cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name,wait.job.name = x)  }else
      {
      cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name)
      }

      cmd15 = paste(cmd.p, cmd14)
      #}else
      #{
      #  cmd15=cmd14
      #}
      system(cmd15,intern = TRUE,ignore.stdout = TRUE)
      cat(cmd15,"\n\n")

    } else
    {
      sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")

      if (!dir.exists(sample.name.out.dir.8))
      {
        dir.create(sample.name.out.dir.8, recursive = TRUE)
      }
      cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                    genome.index, y[1], sep = " ")

      #if(length(m.id)==1){
        job.name <- paste0("Alignment.",i)

        if(!is.null(wait.job.name))
        {
          x <- paste0(wait.job.name,".",i)
          cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name,wait.job.name = x)  }else
          {
            cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name)
          }

        cmd25 = paste(cmd.p, cmd24)

      #}else
      #{
        cmd25=cmd24
      #}
      system(cmd25,intern = TRUE,ignore.stdout = TRUE)
      cat(cmd25,"\n\n")
    }

  }

}

#' bsub -P bbc -J "alignment[1-8]" -o %J.alignment.%I.log -e %J.alignment.%I.err -W 72:00 -n 16 -q parallel -R 'rusage[mem= 25000 ] span[ptile= 8 ]' -u aimin.yan@med.miami.edu "R -e 'library(ChipSeq);library(ThreeUTR);ThreeUTR:::alignmentUseJobArray(\"/scratch/projects/bbc/aiminy_project/DoGsFastq\",\"/scratch/projects/bbc/aiminy_project/DoGs_AlignmentBamTophatGeneral2\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf\",\"/projects/ctsi/bbc/Genome_Ref/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome\")'"

alignmentUseJobArray <- function(input.fastq.files.dir, output.dir, gene.model.file = NULL,
                                 genome.index, cmd.input="parallel",wait.job.name=NULL)
{

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

  re <- parserreadfiles(input.fastq.files.dir, "fastq")
  res <- re$input

  #print(res)
  xx <- lapply(res, function(u)
  {
    path_name = dirname(u)

    file_name = file_path_sans_ext(basename(u))

    if (regexpr(pattern = "_", file_name) != -1)
    {
      # cat('match:',file_name,'\n')

      p <- regexpr(pattern = "_", file_name)
      pp <- p - 1
      x <- substr(file_name, 1, pp)

    } else
    {
      # cat('no match:',file_name,'\n')
      x <- file_name
    }

    x
  })

  xxx <- unique(unlist(xx))
  res2 <- unlist(res)

  print(xxx)
  print(res2)

  #print(xxx[index])
  #print(xxx[as.integer(index)])

  m.id <- grep("login", system("hostname", intern = TRUE))
  cmd6 = "tophat --library-type fr-firststrand -g 1 -G"
  cmd8 = "-p 4 -o"

  #for (i in 1:length(xxx))
  #{
    i <- as.integer(index)
    sample.name <- xxx[i]
    sample.name.out.dir <- file.path(output.dir, sample.name)

    if (!dir.exists(sample.name.out.dir))
    {
      dir.create(sample.name.out.dir, recursive = TRUE)
    }

    y <- res2[grep(xxx[i], res2)]

    if (length(y) == 2)
    {
      yy1 <- basename(y[1])
      yy2 <- basename(y[2])

      p1 <- regexpr(pattern = "_", yy1)
      pp1 <- p1 + 1
      x1 <- substr(yy1, pp1, pp1)

      p2 <- regexpr(pattern = "_", yy2)
      pp2 <- p2 + 1
      x2 <- substr(yy2, pp2, pp2)


      sample.name.out.dir.3 = file.path(sample.name.out.dir, paste0("Fs",
                                                                    x1, x2))

      if (!dir.exists(sample.name.out.dir.3))
      {
        dir.create(sample.name.out.dir.3, recursive = TRUE)
      }

      cmd14 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.3,
                    genome.index, y[1], y[2], sep = " ")
      #if(length(m.id)==1){

      # job.name <- paste0("Alignment.",i)
      # if(!is.null(wait.job.name))
      # {
      #   x <- paste0(wait.job.name,".",i)
      #   cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name,wait.job.name = x)  }else
      #   {
      #     cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name)
      #   }
      #
      # cmd15 = paste(cmd.p, cmd14)
      #}else
      #{
       cmd15=cmd14
      #}
      system(cmd15,intern = TRUE,ignore.stdout = TRUE)
      cat(cmd15,"\n\n")

    } else
    {
      sample.name.out.dir.8 = file.path(sample.name.out.dir, "Fs")

      if (!dir.exists(sample.name.out.dir.8))
      {
        dir.create(sample.name.out.dir.8, recursive = TRUE)
      }
      cmd24 = paste(cmd6, gene.model.file, cmd8, sample.name.out.dir.8,
                    genome.index, y[1], sep = " ")

      #if(length(m.id)==1){
      # job.name <- paste0("Alignment.",i)
      #
      # if(!is.null(wait.job.name))
      # {
      #   x <- paste0(wait.job.name,".",i)
      #   cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name,wait.job.name = x)  }else
      #   {
      #     cmd.p <- ChipSeq:::usePegasus(cmd.input,"72:00",16,25000,8,job.name)
      #   }
      #
      # cmd25 = paste(cmd.p, cmd24)

      #}else
      #{
      cmd25=cmd24
      #}
      system(cmd25,intern = TRUE,ignore.stdout = TRUE)
      cat(cmd25,"\n\n")
    }

  #}

}

convertSra2FastqUseJobArray <- function(sra.file.dir, output.dir)
{

  re <- parserreadfiles(sra.file.dir, "sra")

  res <- re$input

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir,recursive = TRUE)
  }

  m.id <- grep("login", system("hostname", intern = TRUE))

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

#  cmd.l <- lapply(1:length(res), function(u, res,m.id,wait.job.name,output.dir)
#  {
    u <- as.integer(index)
    cmd0 <- "fastq-dump --split-3"

    path_name = dirname(res[[u]])

    file_name = file_path_sans_ext(basename(res[[u]]))

    cmd1 <- paste(cmd0, res[[u]], "-O", output.dir, sep = " ")

    # if(length(m.id) == 1){
    #   job.name <-  paste0("sra2fastq.",u)
    #
    #   if(!is.null(wait.job.name)){
    #     wait.job.name <- wait.job.name
    #     cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
    #                                   cores = 32, Memory = 25000, span.ptile = 16, job.name=job.name,wait.job.name=wait.job.name)
    #   }else
    #   {
    #     cmd.p <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00",
    #                                   cores = 32, Memory = 25000, span.ptile = 16, job.name=job.name)
    #   }
    #   cmd2 <- paste(cmd.p,cmd1,sep = "")
    # }else
    cmd2 <- cmd1

    system(cmd2)
    cat(cmd2,"\n")
 # }, res,m.id,wait.job.name,output.dir)

  #re <- list(cmdl = cmd.l, output.dir = output.dir)

  #re

}

useJobArrayOnPegasus <- function(job.option=c("general","parallel","bigmem")
         , Wall.time, cores, Memory, span.ptile,job.name,wait.job.name=NULL) {

  job.option <- match.arg(job.option)

  #index.job <- regexpr("\\[",job.name)[1]

  #job.name.array <- substr(job.name,1,(index.job-1))

  job.name.array <- job.name

  switch (job.option,
          parallel = {
            cmd0 = paste(Wall.time,"-n",cores,"-q parallel -R 'rusage[mem=",Memory,"] span[ptile=",span.ptile,"]' -u aimin.yan@med.miami.edu",sep = " ")
          },
          bigmem = {
            cmd0 = paste(Wall.time,"-n",cores,"-q bigmem -R 'rusage[mem=",Memory,"] span[ptile=", span.ptile, "]' -u aimin.yan@med.miami.edu",sep = " ")
          },
          general = {
            cmd0 = paste(Wall.time,"-n",cores,"-q general -R 'rusage[mem=",Memory,"] span[ptile=",span.ptile, "]' -u aimin.yan@med.miami.edu",sep = " ")
          }
  )

  if(!is.null(wait.job.name)){
    cmd1 = paste0("bsub -w \"done(\"", wait.job.name, "\")\"", " -P bbc -J \"",
                  job.name, paste0("\" -o %J.", job.name.array, ".log "), paste0("-e %J.",
                                                                           job.name.array, ".err -W"))
  }else{
    cmd1 = paste0("bsub -P bbc -J \"",job.name, paste0("\" -o %J.", job.name.array, ".log "), paste0("-e %J.",job.name.array, ".err -W"))
  }

  cmd = paste(cmd1,cmd0,sep=" ")

  return(cmd)
}

createBsubJobArrayRfun <- function(Rfun,job.name,wait.job.name){
  x <- useJobArrayOnPegasus("parallel","72:00",16,25000,8,job.name,wait.job.name)
  xx <- paste(x,paste0("\"R -e ",paste0("\'",Rfun,"\'"),"\""),sep=" ")
  xx
}
#test <- createBubRfun(Rfun,"sra2fastq[1-8]","wgetDownload")
#system(test)

processBamFilesUseJobArray <- function(input.alignment.dir, output.dir)
{
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
                                           pattern = "accepted_hits.bam"))

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

#  cmd0 = "72:00 -n 8 -q general -u aimin.yan@med.miami.edu"

#  cmd1 = "bsub -P bbc -J \"bamProcess\" -o %J.bamProcess.log -e %J.bamProcess.err -W"

  cmd2 <- "cp"

#  y <- lapply(re, function(u)
#  {
    # /SRR2038198/Fs12/accepted_hits.bam

    u <- as.integer(index)

    dir.name.0 <- dirname(re[[u]])
    dir.name.1 <- dirname(dir.name.0)

    x <- basename(dir.name.0)
    y <- basename(dir.name.1)
    file.name <- basename(re[[u]])

    sample.name <- paste(y, x, file.name, sep = "-")

    cmd <- paste(cmd2, re[[u]], file.path(output.dir, sample.name),
                 sep = " ")

    system(cmd,intern = TRUE)
    cat(cmd,"\n\n")
  #})

  #print(y)

  # yyy <- lapply(y, function(u)
  # {
  #   yy <- system(u, intern = TRUE)
  #   yy
  # })
  #
  # yyyy <- unlist(yyy)
  #
  # return(yyyy)
}

convertBam2bedUsingJobArray <- function(input.bamfile.dir, output.bedfile.dir)
{
  res <- parserreadfiles(input.bamfile.dir, "bam")

  res <- res$input

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

  #m.id <- grep("login", system("hostname", intern = TRUE))

  # if (!dir.exists(output.bedfile.dir)) { dir.create(output.bedfile.dir,
  # recursive = TRUE) }

  #output.bedfile.dir <- file.path(output.bedfile.dir, "BedFileFromBam")

  if (!dir.exists(output.bedfile.dir))
  {
    dir.create(output.bedfile.dir, recursive = TRUE)
  }

 # cmd.l <- lapply(1:length(res), function(u, m.id, res, output.bedfile.dir)
#  {
    # cat(u,'\n') cmd9 <- 'grep' cmd10 <- '~/PathwaySplice/inst/extdata/' cmd11
    # <- '/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt' cmd12 <- '>' cmd13
    # <- paste0('/Counts.',n,'.genes.txt') xxx <- gsub(';','',xx)
    u <- as.integer(index)
    path_name = dirname(res[[u]])
    path_name2 <- basename(path_name)

    file_name = file_path_sans_ext(basename(res[[u]]))

    file_name <- paste0(path_name2, "-", file_name)
    # if (m.id == 1)
    # {
    #   job.name <- paste0("bam2bed.", u)
    #   cmd0 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
    #                                Memory = 25000, span.ptile = 16, job.name)

      cmd1 <- "bedtools bamtobed -i"
      cmd2 <- ">"

      cmd3 <- paste(cmd1, res[[u]], cmd2, file.path(output.bedfile.dir,
                                                          paste0(file_name, ".bed")), sep = " ")
    # } else
    # {
    #   cmd3 <- paste(cmd1, res[[u]], cmd2, file.path(output.bedfile.dir,
    #                                                 paste0(file_name, ".bed")), sep = " ")
    # }

    system(cmd3)

    cat(cmd3,"\n\n")
 # }, m.id, res, output.bedfile.dir)

#  re <- list(cmdl = cmd.l, output.bedfile.dir = output.bedfile.dir)

#  re

}

removeReadsOnExonIntronUsingJobArray <- function(input.bed.file.dir, annotation.bed.file.dir,
                                    output.bed.file.dir)
{
  re <- parserreadfiles(input.bed.file.dir, "bed")

  res <- re$input

  annotationBed <- parserreadfiles(annotation.bed.file.dir,"bed",sample.group=c("hg19_exons.bed","hg19_intron.bed"))

  #m.id <- grep("login", system("hostname", intern = TRUE))

  if (!dir.exists(output.bed.file.dir))
  {
    dir.create(output.bed.file.dir, recursive = TRUE)
  }

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")


#  cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
#                                          span.ptile, res, annotationBed, output.bed.file.dir)
#  {

    u <- as.integer(index)

    file_name = file_path_sans_ext(basename(res[[u]]))

    exon.intron <- paste(unlist(annotationBed$input),collapse=" ")

    # if (m.id == 1)
    # {
    #   job.name = paste0("bedRmExonIntron.", u)
    #   cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
    #                                Memory = 16000, span.ptile = 16, job.name)
    #   exon.intron <- paste(unlist(annotationBed$input),collapse=" ")
    #   cmd2 = paste("bedtools intersect -v -a", res[[u]], "-b", exon.intron,
    #                "\\>", file.path(output.bed.file.dir, paste0(file_name, "_rm_exon_intron.bed")),
    #                sep = " ")
    #   cmd3 = paste(cmd1, cmd2, sep = " ")
    # } else

      cmd3 = paste("bedtools intersect -v -a", res[[u]], "-b", exon.intron,
                   ">", file.path(output.bed.file.dir, paste0(file_name, "_rm_exon_intron.bed")),
                   sep = " ")
#}

    cmd <- cmd3

    cat(cmd, "\n\n")

    system(cmd)

    cmd
 # }, m.id, Wall.time, cores, Memory, span.ptile, res, annotationBed, output.bed.file.dir)

}

getCount4DownstreamUsingJobArray <- function(input.bed.file.dir, annotation.bed.file.dir,
                                output.count.file.dir)
{

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

  re <- parserreadfiles(input.bed.file.dir, "bed")

  res <- re$input

  annotationBed <- parserreadfiles(annotation.bed.file.dir,"bed",sample.group=c("hg19_gene.bed"))

  #m.id <- grep("login", system("hostname", intern = TRUE))
  exon.intron <- paste(unlist(annotationBed$input),collapse=" ")

  if (!dir.exists(output.count.file.dir))
  {
    dir.create(output.count.file.dir, recursive = TRUE)
  }

  #cmd.l <- lapply(1:length(res), function(u, m.id, Wall.time, cores, Memory,
  #                                        span.ptile, res, annotationBed, output.count.file.dir)
  #{
    u <- as.integer(index)

    file_name = file_path_sans_ext(basename(res[[u]]))

    # if (m.id == 1)
    # {
    #   job.name = paste0("bed2count.", u)
    #   cmd1 <- ChipSeq:::usePegasus("parallel", Wall.time = "72:00", cores = 32,
    #                                Memory = 16000, span.ptile = 16, job.name)
    #   exon.intron <- paste(unlist(annotationBed$input),collapse=" ")
    #   cmd2 = paste("bedtools window -a",exon.intron,"-b",res[[u]],"-l 0 -r 45000 -sw -c",
    #                "\\>", file.path(output.count.file.dir, paste0(file_name, "_downstream_count.txt")),
    #                sep = " ")
    #   cmd3 = paste(cmd1, cmd2, sep = " ")
    # } else
    #{
      cmd3 = paste("bedtools window -a",exon.intron,"-b",res[[u]],"-l 0 -r 45000 -sw -c",
                   ">", file.path(output.count.file.dir, paste0(file_name, "_downstream_count.txt")),
                   sep = " ")
    #}

    cmd <- cmd3

    cat(cmd, "\n\n")

    system(cmd)

  #  cmd
#  }, m.id, Wall.time, cores, Memory, span.ptile, res, annotationBed, output.count.file.dir)

}

CountAndDE <- function(output.count.dir,sample.info.file,output.res.dir) {

  if (!dir.exists(output.res.dir))
  {
    dir.create(output.res.dir, recursive = TRUE)
  }

  x <- read.table(file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"),header = TRUE)

  g.1 <- colnames(x)[2]
  g.2 <- unique(as.character(x[,2]))[2]
  g.3 <- unique(as.character(x[,2]))[1]

  res <- convertCountFile2Table(output.count.dir,"*.txt")

  res.new <- ThreeUTR:::matchAndDE(res,sample.info.file,group.comparision = c(g.1,g.2,g.3))

  sumResult(res.new,output.res.dir)

  write.table(res.new$re.FC.sorted,file = file.path(output.res.dir,"Results.csv"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)

  return(res.new)

}

#' res <- ThreeUTR:::convertCountFile2Table("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts","*.txt")

#' library(org.Hs.eg.db)
#' res <- ThreeUTR:::CountAndDE("~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Counts",file.path(system.file("extdata", package = "ThreeUTR"),"sample_infor.txt"),"~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Results")

#' output.res.dir <- "~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Results"
#'
#' res.sum <- ThreeUTR:::sumResult(res,"~/Dropbox (BBSR)/Aimin_project/Research/DoGs/Results")
#'
sumResult <-function(res,output.res.dir){

  if (!dir.exists(output.res.dir))
  {
    dir.create(output.res.dir, recursive = TRUE)
  }

  re <- res$re.DESeq
  table(re$padj<0.05)
  ## Order by adjusted p-value
  re <- re[order(re$padj), ]
  ## Merge with normalized count data
  redata <- merge(as.data.frame(re), as.data.frame(counts(res$dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(redata)[1] <- "Gene"
  head(redata)

  png(file.path(output.res.dir,"diffexpr-pvalue.png"))
  hist(re$pvalue, breaks=50, col="grey",main = "Histogram of p value from differential analysis")
  dev.off()

  attr(re,"metadata")$filterThreshold
  #attr(re,"metadata")$filterNumRej
  png(file.path(output.res.dir,"independent-filtering.png"))
  plot(attr(re,"metadata")$filterNumRej, type="b", xlab="quantiles of baseMean", ylab="number of rejections")
  dev.off()

  # Regularized log transformation for clustering/heatmaps, etc
  rld <- rlogTransformation(res$dds)
  head(assay(rld))
  hist(assay(rld))

  # Colors for plots below
  ## Ugly:
  ## (mycols <- 1:length(unique(condition)))
  ## Use RColorBrewer, better
  library(RColorBrewer)

  condition <- res$colData[,1]
  n <- length(res$colData[,1])

  (mycols <- brewer.pal(n, "Dark2")[1:length(unique(condition))])

  # Sample distance heatmap
  sampleDists <- as.matrix(dist(t(assay(rld))))
  library(gplots)
  png(file.path(output.res.dir,"qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(10, 10), main="Sample Distance Matrix")
  dev.off()

  png(file.path(output.res.dir,"qc-pca-1.png"), 1000, 1000, pointsize=20)
  rld_pca(rld,colors=mycols,intgroup=colnames(res$colData)[1],legendpos="bottomleft",xlim=c(-30, 30),ylim=c(-20,30))
  dev.off()

  png(file.path(output.res.dir,"qc-pca-2.png"), 1000, 1000, pointsize=20)
  DESeq2::plotPCA(rld, intgroup=colnames(res$colData)[1])
  dev.off()

  png(file.path(output.res.dir,"diffexpr-maplot.png"), 1500, 1000, pointsize=20)
  maplot(redata, main="MA Plot")
  dev.off()

  png(file.path(output.res.dir,"diffexpr-volcanoplot-1.png"), 1200, 1000, pointsize=20)
  volcanoplot(redata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-8, 8))
  dev.off()

  png(file.path(output.res.dir,"diffexpr-volcanoplot-2.png"), 1200, 1000, pointsize=20)
  deawVolcano(re)
  dev.off()
}

maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)

#  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))

  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=NULL,cex=textcx))

  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}

volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main=main, ...))
  abline(h=1.3, col = "blue")
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(padj), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(padj), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
 #   with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(padj), labs=Gene, cex=textcx, ...))

    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(padj), labs=NULL, cex=textcx, ...))

  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|Log2FoldChange|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))

}



processSpliceJunctionFilesUseJobArray <- function(input.alignment.dir, output.dir)
{
  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
                                           pattern = "junctions.bed"))

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")


  cmd2 <- "cp"

  u <- as.integer(index)

  dir.name.0 <- dirname(re[[u]])
  dir.name.1 <- dirname(dir.name.0)

  x <- basename(dir.name.0)
  y <- basename(dir.name.1)
  file.name <- basename(re[[u]])

  sample.name <- paste(y, x, file.name, sep = "-")

  cmd <- paste(cmd2, re[[u]], file.path(output.dir, sample.name),
               sep = " ")

  system(cmd,intern = TRUE)
  cat(cmd,"\n\n")
}

#' res.sum <- ThreeUTR:::generateSmallBedFile(res,strand.specific=TRUE,"~/Dropbox (BBSR)/Aimin_project/Research/DoGs/SmallExample")
#'
generateSmallBedFile <- function(res,strand.specific=FALSE,output.res.dir){

  if (!dir.exists(output.res.dir))
  {
    dir.create(output.res.dir, recursive = TRUE)
  }


  re <- res$re.FC.sorted
  table(re$padj<0.05)

  re1 <- re[which(re$padj<0.05),]

  gene.list <- unique(re1$gene)
  gene.list.choose <- sample(gene.list,100)

  re2 <- re1[which(re1$gene %in% gene.list.choose),]

  re3 <- re2[,which(colnames(re2) %in% c("seqnames","start","end","strand"))]

  if(strand.specific == TRUE){
  re3[which(re3$strand == "+"),]$end <- re3[which(re3$strand == "+"),]$end+45000
  re3[which(re3$strand == "-"),]$start <- re3[which(re3$strand == "-"),]$start+45000
}
  re4 <- re3[order(re3$seqnames),]

  write.table(re4[,1:3],file = file.path(output.res.dir,"sample_genes.bed"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

  re <- list(re2=re2,re3=re3,re4=re4)

  return(re)

}

generateSubSetBam <- function(input.alignment.dir,annotation.bed.file.dir,output.dir){

  if (!dir.exists(output.dir))
  {
    dir.create(output.dir, recursive = TRUE)
  }

  annotationBed <- parserreadfiles(annotation.bed.file.dir,"bed",sample.group=c("sample_genes.bed"))

  annotation.gene <- paste(unlist(annotationBed$input),collapse=" ")

  re <- file.path(input.alignment.dir, dir(input.alignment.dir, recursive = TRUE,
                                           pattern = "*.bam"))

  index <- system('echo $LSB_JOBINDEX',intern = TRUE)
  total <- system('echo $LSB_JOBINDEX_END',intern = TRUE)

  cat(index,"\n\n")
  cat(total,"\n\n")

  cmd2 <- "intersectBed -abam"

  u <- as.integer(index)

  dir.name.0 <- dirname(re[[u]])
  dir.name.1 <- dirname(dir.name.0)

  x <- basename(dir.name.0)
  y <- basename(dir.name.1)
  file.name <- basename(re[[u]])

  sample.name <- paste(y, x, file.name, sep = "-")

  sample.name <- file_path_sans_ext(sample.name)

  cmd <- paste(cmd2,re[[u]],"-b",annotation.gene,">",file.path(output.dir, paste0(sample.name,"_region.bam")),sep = " ")

  system(cmd,intern = TRUE)
  cat(cmd,"\n\n")

}

#results <- res$re.DESeq
#deawVolcano(results)

deawVolcano <- function(results) {
  library(dplyr)
  library(ggplot2)
  data <- data.frame(gene = row.names(results),
                     pvalue = -log10(results$padj),
                     lfc = results$log2FoldChange)
  data <- na.omit(data)

  data <- data %>% mutate(color = ifelse(data$lfc > 0 & data$pvalue > 1.3,
                          yes = "Treated",
                          no = ifelse(data$lfc < 0 & data$pvalue > 1.3,
                                      yes = "Untreated",
                                      no = "none")))


  # Color corresponds to fold change directionality
  colored <- ggplot(data, aes(x = lfc, y = pvalue)) +
    geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + # add gene points
    theme_bw(base_size = 16) + # clean up theme
    theme(legend.position = "none") + # remove legend
    ggtitle(label = "Volcano Plot", subtitle = "Colored by directionality") +  # add title
    xlab(expression(log[2]("Treated" / "Untreated"))) + # x-axis label
    ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
    geom_vline(xintercept = 0, colour = "black") + # add line at 0
    geom_hline(yintercept = 1.3, colour = "black") + # p(0.05) = 1.3
    annotate(geom = "text",
             label = "Untreated",
             x = -2, y = 85,
             size = 7, colour = "black") + # add Untreated text
    annotate(geom = "text",
             label = "Treated",
             x = 2, y = 85,
             size = 7, colour = "black") + # add Treated text
    scale_color_manual(values = c("Treated" = "#E64B35",
                                  "Untreated" = "#3182bd",
                                  "none" = "#636363")) # change colors

  # Plot figure
  #colored
  colored + scale_y_continuous(trans = "log1p")
}
