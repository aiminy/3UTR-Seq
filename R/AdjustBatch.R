#' AdjustBatch
#'
#' @param re.rMAT.com.total.test.2
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"
#'
#' re.3UTR<-AdjustBatch(re.rMAT.com.total.test.2,dir.name,"By_P","By_Gene","Interested_gene_only","Interested_gene_only_sig")
#'
#' re.3UTR.one.side<-AdjustBatch(re.rMAT.com.total.test.2,dir.name,"By_P_one_side","By_Gene_one_side",
#' "Interested_gene_only_one_side","Interested_gene_only_sig_one_side")
#'
#'
#' re.3UTR.one.side.norm<-AdjustBatch(re.rMAT.com.total.test.2,DE.norm.with.rpkm.norm.2,dir.name,"By_P_one_side_norm_2","By_Gene_one_side_norm_2",
#' "Interested_gene_only_one_side_norm_2","Interested_gene_only_sig_one_side_norm_2")
#'
AdjustBatch <- function(re.rMAT.com.total.test.2,DE.norm.with.rpkm.norm.2,dir.name,out_by_p,out_by_gene,out_gene_interested,out_gene_interested_sig) {

  #names(re.rMAT.com.total.test.2)

  df.NT<-counts(re.rMAT.com.total.test.2$dds)
  colData<-colData(re.rMAT.com.total.test.2$dds)

  cell<-factor(rep(c('emp','hela'),c(2,2)))
  cell=rep(cell,2)

  colData.1 <- data.frame(condition=factor(rep(c('Dox', 'WT'),c(4, 4))),cell=cell)

  #Use original design, not adjusted by any batch factor
  dds <- DESeqDataSetFromMatrix(df.NT,colData.1,formula(~condition))
  colData(dds)

  print(design(dds))
  dds <- DESeq(dds)
  res.dds<- results(dds)
  head(res.dds)
  print(sizeFactors(estimateSizeFactors(dds)))
  print(head(counts(dds,normalized = TRUE)))

  #Check hidden batch effects
  #dds <- estimateSizeFactors(dds)
  dds <- DESeqDataSetFromMatrix(df.NT, colData.1, formula(~condition))
  dat <- counts(dds)

  #idx <- rowMeans(dat) > 1
  #dat <- dat[idx,]

  mod <- model.matrix(~condition,colData(dds))
  mod0 <- model.matrix(~1,colData(dds))
  svseq.NT <- svaseq(dat, mod, mod0)

  par(mfrow=c(2,1),mar=c(3,5,3,1))
  stripchart(svseq.NT$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
  abline(h=0)
  stripchart(svseq.NT$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
  abline(h=0)

  #Adjust by surrogate variables 1 and 2
  ddssva <- dds
  ddssva$SV1 <- svseq.NT$sv[,1]
  ddssva$SV2 <- svseq.NT$sv[,2]
  design(ddssva) <- ~ SV1 + SV2 + condition
  ddssva <- DESeq(ddssva)
  print(sizeFactors(estimateSizeFactors(ddssva)))
  print(head(counts(ddssva,normalized = TRUE)))
  res.ddssva.2.sva.adjusted <- results(ddssva)
  head(res.ddssva.2.sva.adjusted)
  summary(res.ddssva.2.sva.adjusted)

  #Adjusted by cell type

  dds.cell <- DESeqDataSetFromMatrix(df.NT,colData.1,design=~ cell+condition)

  #design(dds.cell) <- ~ cell + condition
  print(design(dds.cell))

  dds.cell <- DESeq(dds.cell)
  res.dds.cell<- results(dds.cell)

  re.condition<-results(dds.cell, contrast=c("condition", "WT", "Dox"))
  head(re.condition)
  summary(re.condition)

  re.cell<-results(dds.cell, contrast=c("cell", "emp", "hela"))
  summary(re.cell)

  head(res.dds.cell)
  print(sizeFactors(estimateSizeFactors(dds.cell)))
  print(head(counts(dds.cell,normalized = TRUE)))

  # dds.cell.2 <- DESeq(dds.cell.2)
  # res.dds.cell.2<- results(dds.cell.2)
  # head(res.dds.cell.2)
  #
  # summary(res.dds.cell.2)
  # summary(res.dds.cell)
  #
  #
  # #Use limma based method
  # #dat.adjusted<-cbind(dat[,c(3,4,7,8)]*(2^re.cell$log2FoldChange),dat[,c(3,4,7,8)])
  # dds.ad <- DESeqDataSetFromMatrix(df.NT,colData.1,design=~ cell+condition)
  # dds.ad <- estimateSizeFactors(dds.ad)
  # norm.counts <- counts(dds.ad, normalized=TRUE)
  # log.norm.counts <- log2(norm.counts + 1)
  # y<-log.norm.counts
  # batch <- c("emp","emp","hela","hela","emp","emp","hela","hela")
  #
  # y2<-removeBatchEffect(y, batch)
  # #par(mfrow=c(1,2))
  # boxplot(as.data.frame(y),main="Original")
  # boxplot(as.data.frame(y2),main="Batch corrected")
  #
  # f.st134.frma <-factor(colData.1$condition)
  # design.st134.frma <- model.matrix(~0+f.st134.frma)
  # colnames(design.st134.frma) <- levels(f.st134.frma)
  #
  # fit.st134.all.probes.frma <- lmFit(y2, design.st134.frma)
  # cont.matrix.st134.frma <- makeContrasts(stWTDox="WT-Dox",levels=design.st134.frma)
  # fit2.st134.all.probes.frma  <- contrasts.fit(fit.st134.all.probes.frma, cont.matrix.st134.frma)
  # fit2.st134.all.probes.frma <- eBayes(fit2.st134.all.probes.frma)
  # TopTableSt34.all.probes.frma<-topTable(fit2.st134.all.probes.frma,coef=1,n=dim(y2)[1],sort.by="p")
  #
  # cutoff <- 0.3
  # wtResCont <- decideTests(fit2.st134.all.probes.frma, p.value = cutoff, method = "global")
  # summary(wtResCont)
  #
  # length(which(TopTableSt34.all.probes.frma$adj.P.Val<0.3&TopTableSt34.all.probes.frma$logFC<0))
  #
  # length(which(TopTableSt34.all.probes.frma$adj.P.Val<0.3&TopTableSt34.all.probes.frma$logFC>0))
  #
  # #matrix of SVs
  # svseq.NT$sv
  # #number of SVs
  # svseq.NT$n.sv
  # save(svseq.NT, file='/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/svseq.NT.RData')
  #
  # colData.sva <- cbind(colData,svseq.NT$sv)
  #
  # colData.sva <- data.frame(condition=factor(rep(c('Dox', 'WT'),
  #                                                c(4, 4)),svseq.NT$sv))
  #
  #
  # # preprocess data
  # colData.sva <- data.frame(condition=factor(rep(c('Dox', 'WT'),
  #                                                c(4, 4))),
  #                           svseq.NT$sv)
  #
  # colnames(colData.sva)[-1] = paste0("SV", 1:svseq.NT$sv)
  #
  # ncol <- dim(svseq.NT$sv)[2]
  #
  # SV1<-paste(paste0("SV", 1:svseq.NT$sv), collapse = "+")
  #
  # design.sva.1 <- formula(paste0("~ SV1 +","condition"))
  #
  # design.sva.2 <- formula(paste0("~ condition +","SV1"))
  #
  #
  # #design.sva <- formula(paste0("~ paste(paste0("SV", 1:svseq.NT$sv), collapse = "+") +", condition ))
  # colData.sva.1<-colData.sva[,c(2,1)]
  #
  # ddssva.NT <- DESeqDataSetFromMatrix(df.NT,colData.sva, design.sva.1)
  # rownames(colData(ddssva.NT)) <- colnames(df.NT)
  # ddssva.NT <- DESeq(ddssva.NT)
  # res.sva.NT <- results(ddssva.NT)
  # head(res.sva.NT)
  #
  # ddssva.NT.1 <- DESeqDataSetFromMatrix(df.NT,colData.sva.1, design.sva.1)
  # rownames(colData(ddssva.NT.1)) <- colnames(df.NT)
  # ddssva.NT.1 <- DESeq(ddssva.NT.1)
  # res.sva.NT.1 <- results(ddssva.NT.1)
  # head(res.sva.NT.1)
  #
  #
  # ddssva.NT.2 <- DESeqDataSetFromMatrix(df.NT,colData.sva, design.sva)
  # rownames(colData(ddssva.NT.2)) <- colnames(df.NT)
  #
  #
  # ddssva.NT.2 <- DESeq(ddssva.NT.2)
  # res.sva.NT.2 <- results(ddssva.NT.2)
  #
  # head(res.sva.NT)

  MatchGene <- function(res.sva.NT,re.rMAT.com.total.test.2) {
    temp<-re.rMAT.com.total.test.2$DE[,c(1,10:24)]
    rownames(temp)<-temp[,1]
    temp2<-merge(as.data.frame(res.sva.NT),temp,by=0,sort=FALSE)
    re.FC<-temp2
    re.FC.sorted<-re.FC[order(re.FC$pvalue),]
    return(re.FC.sorted)
  }


  re.FC.sorted<-MatchGene(res.ddssva.2.sva.adjusted,re.rMAT.com.total.test.2)

  DE.norm.with.rpkm.norm.3<-apply(DE.norm.with.rpkm.norm.2,2,as.numeric)

  rownames(DE.norm.with.rpkm.norm.3)<-rownames(DE.norm.with.rpkm.norm.2)

  DE.rpkm.norm<-cbind(DE.norm.with.rpkm.norm.3,DE.norm.with.rpkm.norm.3[,9:16]*(1000/4500))

  rownames(re.FC.sorted)=re.FC.sorted$Row.names

  re.FC.sorted<-merge(re.FC.sorted[,-1],DE.rpkm.norm,by=0,sort=FALSE)

  p.one.sided<-unlist(lapply(re.FC.sorted$stat,function(x){convert.z.score(x,two.side=F)}))

  p.one.sided.adjust<-p.adjust(p.one.sided,method = "BH")

  re.FC.sorted<-cbind(re.FC.sorted,p.one.sided,p.one.sided.adjust)

  re.FC.sorted.sorted.by.gene<- re.FC.sorted[order(re.FC.sorted$gene),]

  write.csv(re.FC.sorted,file=paste0(dir.name,"3UTR_DE_",out_by_p,".csv"),row.names = FALSE,quote=FALSE)

  write.csv(re.FC.sorted.sorted.by.gene,file=paste0(dir.name,"3UTR_DE_",out_by_gene,".csv"),row.names = FALSE,quote=FALSE)

  target<-re.FC.sorted[which(re.FC.sorted$pvalue<0.05&re.FC.sorted$log2FoldChange<0),]$gene

  #which(target %in% genes.interested[,1])

  re.genes.interested<-re.FC.sorted.interested<-re.FC.sorted[which(re.FC.sorted$gene %in% genes.interested[,1]),]

  re.genes.interested.sig<-re.FC.sorted.interested[which(re.FC.sorted.interested$pvalue<0.05&re.FC.sorted.interested$log2FoldChange<0),]

  write.csv(re.genes.interested,file=paste0(dir.name,"3UTR_DE_",out_gene_interested,".csv"),row.names = FALSE,quote=FALSE)

  write.csv(re.genes.interested.sig,file=paste0(dir.name,"3UTR_DE_",out_gene_interested_sig,".csv"),row.names = FALSE,quote=FALSE)

  #re.com.total.pos<-re.FC.sorted[which(re.FC.sorted$`FoldChange(WT/Dox)`<1),]
  #re.com.total.neg<-re.FC.sorted[which(re.FC.sorted$`FoldChange(WT/Dox)`>=1),]

  re<-list(re.FC.sorted=re.FC.sorted,
           re.FC.sorted.sorted.by.gene=re.FC.sorted.sorted.by.gene,
           re.genes.interested=re.genes.interested,
           re.genes.interested.sig=re.genes.interested.sig)

  return(re)
}
