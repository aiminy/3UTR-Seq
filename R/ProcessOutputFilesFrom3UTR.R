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
#' dir.name="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/"
#'
#' #new data
#' dir.name="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"
#' input.file.pattern="*downstream.count.hg19.strand.based.txt"
#'
#' com.input.file.pattern="*downstream.com.count.hg19.strand.based.txt"
#'
#' total.com.input.file.pattern="*.gene.downstream.count.hg19.strand.based.total2.txt"
#'
#' #use intergenic reads as normal factor
#' normal.factor="*.rm.exon.intron.hg19.bed"
#'
#' #use total reads as normal factor
#' normal.factor="*.bed"
#'
#' sink("testcode.txt")
#' re.rMAT<-ProcessOutputFilesFrom3UTR(dir.name,input.file.pattern,normal.factor,"intergenic")
#' re.rMAT.com<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,"com_intergenic")
#'
#' re.rMAT.com.total<-ProcessOutputFilesFrom3UTR(dir.name,total.com.input.file.pattern,normal.factor,"com_total")
#'
#' re.rMAT.com.total.interested<-re.rMAT.com.total$DE[which(re.rMAT.com.total$DE$gene %in% genes.interested[,1]),]
#'
#' re.com.total.pos<-re.rMAT.com.total.interested[which(re.rMAT.com.total.interested$`FoldChange(WT/Dox)`<1),]
#' re.com.total.neg<-re.rMAT.com.total.interested[which(re.rMAT.com.total.interested$`FoldChange(WT/Dox)`>=1),]
#'
#'
#' re.com.total.pos[which(re.com.total.pos$gene %in% c("TPCN2")),]
#'
#' re.rMAT.com.interested<-re.rMAT.com$DE[which(re.rMAT.com$DE$gene %in% genes.interested[,1]),]
#' re.pos<-re.rMAT.com.interested[which(re.rMAT.com.interested$`FoldChange(WT/Dox)`<1),]
#' re.neg<-re.rMAT.com.interested[which(re.rMAT.com.interested$`FoldChange(WT/Dox)`>=1),]
#'
#'
#' re.rMAT.total<-ProcessOutputFilesFrom3UTR(dir.name,input.file.pattern,normal.factor,"total")
#' re.rMAT.com.total<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,"com_total")
#'
#' sink()
#'
#' sink("test_com_total.txt")
#' re.rMAT.com.total.test<-ProcessOutputFilesFrom3UTR(dir.name,com.input.file.pattern,normal.factor,"com_total_test")
#' sink()
#'
#' re.rMAT.com.total.test.2<-ProcessOutputFilesFrom3UTR(dir.name,total.com.input.file.pattern,normal.factor,"com_total_test_3")
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
#' col.name<-c(paste0("D.",sapply(strsplit(colnames(DE.norm.with.rpkm.norm.2)[1:8],"\\."),"[[",2)),
#' paste0("R.",sapply(strsplit(colnames(DE.norm.with.rpkm.norm.2)[9:16],"\\."),"[[",3)))
#'
#' col.name.2<-gsub("2016-02-10-","",col.name)
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
#' abline(model, col = "red")
#'}
#'
#' plot(apply(DE.norm.with.rpkm.norm.2,2,as.numeric)[,1:3])
#'
#' M <- cor(DE.norm.with.rpkm.norm.3[,1:8])
#'
#' corrplot.mixed(M)
#'
#' corrplot(M,method="ellipse",order = "hclust", addrect = 2)
#'
#' corrplot(M,method="number",order = "hclust", addrect = 2)
#'
#' corrplot(M,method="number",order = "FPC", addrect = 2)
#'
#
#' corrplot(M,method="ellipse")
#'
#' corrplot(M,method="ellipse",type="upper",)
#' corrplot(M, method="number")
#'
#' corrplot(M,method="ellipse",order="hclust", addrect=4)

#'
#' size.factor.2<-estimateSizeFactors(re.rMAT.com.total.test.2$dds)
#' sizeFactors(size.factor.2)

#corrplot(M, order="hclust", addrect=2, col="wb", bg="gold2")


ProcessOutputFilesFrom3UTR<-function(dir.name,input.file.pattern,normal.factor,out){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",9)

  print(file.name.2)

  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=F)
    colnames(re)=c("Count","GeneName")
    re
  })

  temp.name<-strsplit(names(file.name.2),split="\\.")

  temp.name.2<-trimws(do.call("rbind",lapply(temp.name,"[[",1)))

  final.filtered.norm<-NormalizedCount(re.out,dir.name,normal.factor)

  final.filtered<-final.filtered.norm$final.filtered
  num.intergenic.reads.2<-final.filtered.norm$num.intergenic.reads.2

  n=length(num.intergenic.reads.2)-1
  #suppressPackageStartupMessages(library(DESeq2))

  print(head(final.filtered))

  countData <- apply(final.filtered[,c(2,6,10,12,4,8,14,16)], 2, as.numeric)

  rownames(countData)<-final.filtered[,1]

  colData <- data.frame(condition=factor(c(rep("Dox",4),rep("WT",4))))
  dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))

  #size.factor<-estimateSizeFactors(dds)

  #print(size.factor)

  re.DESeq<-results(DESeq(dds))

  print(re.DESeq)

  re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
  colnames(re.FC)[7]="FoldChange(WT/Dox)"

  txs.gene<-ReformatTxsGene()

  re.FC<-merge(re.FC,txs.gene$txs_genes_DF_2,by=0)

  #print(head(final.filtered))

  colnames(final.filtered)[1]<-"tx_name"

  final.filtered.2<-cbind(final.filtered[,grep(glob2rx("tx_name"),colnames(final.filtered))],
  final.filtered[,grep(glob2rx("Normalized*emp*Dox*"),colnames(final.filtered))],
  final.filtered[,grep(glob2rx("Normalized*hela*dox*"),colnames(final.filtered))],
  final.filtered[,grep(glob2rx("Normalized*emp*WT*"),colnames(final.filtered))],
  final.filtered[,grep(glob2rx("Normalized*hela*wt*"),colnames(final.filtered))])

  colnames(final.filtered.2)[1]<-"tx_name"

  re.FC<-merge(re.FC,final.filtered.2,by="tx_name")

  re.FC.sorted<-re.FC[order(re.FC$pvalue),]

  input.file.pattern.2<-sub("\\*","_",input.file.pattern)

  write.csv(re.FC.sorted,file=paste0(dir.name,"3UTR_DE_",out,input.file.pattern.2,".csv"))

  re2<-cbind(trimws(do.call("rbind",lapply(num.intergenic.reads.2[1:n],"[[",1))),do.call("rbind",lapply(num.intergenic.reads.2[1:n],"[[",9)))

  re.out.3<-list(ReadDownstream45kb=re.out,intergenic_reads=re2,MergedNorma=final.filtered,DE=re.FC.sorted,dds=dds)

  return(re.out.3)

}
