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
#' dir.name="/media/H_driver/2016/Ramin_azhang/"
#'
#' input.file.pattern="*downstream.count.plus.minus.and.minus.plus.3.FPKM.xls"
#'
#' input.file.pattern.4.interested.gene="final_list.csv"
#'
#' out_put="for_bioinfo_core/RNA_seq/Results4Check/"
#'
#' out_by_p=paste0(out_put,"DE_by_p.csv")
#' out_by_gene=paste0(out_put,"DE_by_gene.csv")
#' out_gene_interested=paste0(out_put,"DE_by_interested_gene.csv")
#' out_gene_interested_sig=paste0(out_put,"DE_by_interested_gene_sig.csv")
#'
#' re<-ProcessOutputFilesFrom3UTR_recount(dir.name,input.file.pattern,input.file.pattern.4.interested.gene,out_by_p,out_by_gene,out_gene_interested,out_gene_interested_sig)
#'

ProcessOutputFilesFrom3UTR_recount<-function(dir.name,input.file.pattern,input.file.pattern.4.interested.gene,out_by_p,out_by_gene,out_gene_interested,out_gene_interested_sig){

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",9)

  #print(file.name.2)

  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=F)
    colnames(re)=c("chrom","st","end","accession","mRNA_size","gene_strand","Frag_count","FPM","FPKM")
    re
  })

  temp.name<-strsplit(names(file.name.2),split="\\.")
  temp.name.2<-trimws(do.call("rbind",lapply(temp.name,"[[",1)))

  re.out.2<-do.call("cbind",lapply(re.out,function(x){x}))

  re.out.2.name<-gsub("gene.downstream.count.plus.minus.and.minus.plus.3.FPKM.xls.","",colnames(re.out.2))

  colnames(re.out.2)<-re.out.2.name

  chrom<-as.data.frame(apply(re.out.2[,grep("chrom",colnames(re.out.2))],1,unique))
  st<-as.data.frame(apply(re.out.2[,grep("\\.st",colnames(re.out.2))],1,unique))
  end<-as.data.frame(apply(re.out.2[,grep("end",colnames(re.out.2))],1,unique))
  accession<-as.data.frame(apply(re.out.2[,grep("accession",colnames(re.out.2))],1,unique))
  mRNA_size<-as.data.frame(apply(re.out.2[,grep("mRNA_size",colnames(re.out.2))],1,unique))
  gene_strand<-as.data.frame(apply(re.out.2[,grep("gene_strand",colnames(re.out.2))],1,unique))

  colnames(chrom)<-"chrom"
  colnames(st)<-"st"
  colnames(end)<-"end"
  colnames(accession)<-"accession"
  colnames(mRNA_size)<-"mRNA_size"
  colnames(gene_strand)<-"gene_strand"

  temp<-cbind.data.frame(chrom,st,end,accession,mRNA_size,gene_strand)

  Frag_count<-re.out.2[,grep("Frag_count",colnames(re.out.2))]
  FPM<-re.out.2[,grep("FPM",colnames(re.out.2))]
  FPKM<-re.out.2[,grep("FPKM",colnames(re.out.2))]

  re.out.3<-cbind(temp,Frag_count,FPM,FPKM)

  countData <- apply(re.out.3[,c(7,9,11,12,8,10,13,14)], 2, as.numeric)

  rownames(countData)<-re.out.3[,4]

  dat <- countData
  idx <- rowMeans(dat) > 1
  dat <- dat[idx,]


  cell<-factor(rep(c('emp','hela'),c(2,2)))
  cell=rep(cell,2)

  colData <- data.frame(condition=factor(rep(c('Dox', 'WT'),c(4, 4))),cell=cell)

  dds <- DESeqDataSetFromMatrix(dat, colData, formula(~ condition))

  #counts(dds)


  mod <- model.matrix(~condition,colData(dds))
  mod0 <- model.matrix(~1,colData(dds))
  svseq.NT <- svaseq(counts(dds), mod, mod0)

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

  #print(sizeFactors(estimateSizeFactors(ddssva)))
  #print(head(counts(ddssva,normalized = TRUE)))

  res.ddssva.2.sva.adjusted <- results(ddssva)
  #head(res.ddssva.2.sva.adjusted)
  #summary(res.ddssva.2.sva.adjusted)

  re.FC<-cbind(as.data.frame(res.ddssva.2.sva.adjusted),2^res.ddssva.2.sva.adjusted[,2])
  colnames(re.FC)[7]="FoldChange(WT/Dox)"

  MatchGene <- function(res.sva.NT,re.rMAT.com.total.test.2) {
    temp<-re.rMAT.com.total.test.2$DE[,c(1,23,24)]
    rownames(temp)<-temp[,1]
    #temp<-temp[,-1]
    temp2<-merge(as.data.frame(res.sva.NT),temp,by=0,sort=FALSE)
    re.FC<-temp2
    re.FC.sorted<-re.FC[order(re.FC$pvalue),]
    return(re.FC.sorted)
  }

  re.FC.sorted<-MatchGene(re.FC,re.rMAT.com.total.test.2)

  rownames(re.out.3)=re.out.3$accession
  re.FC.sorted<-re.FC.sorted[,-1]

  rownames(re.FC.sorted)=re.FC.sorted$tx_name

  re.out.3.temp<-re.out.3[,c(1:6,7,9,11,12,8,10,13,14,15,17,19,20,16,18,21,22,23,25,27,28,24,26,29,30)]

  re.FC.sorted<-merge(re.FC.sorted,re.out.3.temp,by=0,sort=FALSE)

  p.one.sided<-unlist(lapply(re.FC.sorted$stat,function(x){convert.z.score(x,two.side=F)}))

  p.one.sided.adjust<-p.adjust(p.one.sided,method = "BH")

  re.FC.sorted<-cbind(re.FC.sorted,p.one.sided,p.one.sided.adjust)

  re.FC.sorted.sorted.by.gene<- re.FC.sorted[order(re.FC.sorted$gene),]

  write.csv(re.FC.sorted,file=paste0(dir.name,out_by_p),row.names = FALSE,quote=FALSE)

  write.csv(re.FC.sorted.sorted.by.gene,file=paste0(dir.name,out_by_gene),row.names = FALSE,quote=FALSE)

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern.4.interested.gene))

  genes.interested<-read.csv(file.name,header=F)

  re.genes.interested<-re.FC.sorted[which(re.FC.sorted$gene %in% genes.interested[,1]),]

  re.genes.interested.sig<-re.genes.interested[which(re.genes.interested$pvalue<0.05&re.genes.interested$log2FoldChange<0),]

  write.csv(re.genes.interested,file=paste0(dir.name,out_gene_interested),row.names = FALSE,quote=FALSE)

  write.csv(re.genes.interested.sig,file=paste0(dir.name,out_gene_interested_sig),row.names = FALSE,quote=FALSE)

  re<-list(re.FC.sorted=re.FC.sorted,
           re.FC.sorted.sorted.by.gene=re.FC.sorted.sorted.by.gene,
           re.genes.interested=re.genes.interested,
           re.genes.interested.sig=re.genes.interested.sig)

  return(re)

}
