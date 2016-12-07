#' FilterTranscripts
#'
#' @return
#' @export
#'
#' @examples
#'
#' Re.txs<-FilterTranscripts("/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/","final_list.csv")
#'
#'
#'
FilterTranscripts <- function(input.dir,input.file.pattern) {

  file.name=paste0(input.dir,dir(input.dir,recursive = TRUE,pattern=input.file.pattern))

  #print(file.name)

  genes.interested<-read.csv(file.name,header=F)

  #print(genes.interested)

  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  library("Homo.sapiens")

  txs <- as.data.frame(transcripts(Homo.sapiens, columns=c("TXNAME","GENEID","SYMBOL")))

  txs.4.genes.interested<-txs[which(txs$SYMBOL %in% genes.interested[,1]),]

  cat(dim(txs.4.genes.interested))

  #print(head(txs.4.genes.interested))

  #
  # columns(TxDb)
  # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  # txs1<-transcripts(txdb)
  #
  # columns(txdb)
  # keys <- c("129787")
  # keys <- keys(hgu95av2.db, 'ENTREZID')
  # select(txdb, keys = keys, columns=c("TXCHROM","TXSTART","TXEND","TXSTRAND","TXNAME"), keytype="GENEID")
  # txs2<-select(txdb,keys = keys,columns=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND","TXNAME"),keytype="GENEID")
  # dim(txs2)
  # txs2
  #
  # txdb
   GRList <- transcriptsBy(txdb, by = "gene")
  # length(GRList)
  # which(names(GRList) %in% c("129787"))
  # GRList[[4237]]
  # as.data.frame(GRList[[4237]])
  #
   num.of.txs.4.gene<-unlist2(lapply(GRList,length))

   boxplot(num.of.txs.4.gene)

   hist(num.of.txs.4.gene)

  # which.max(num.of.txs.4.gene)
  #
  # index<-which(names(unlist(num.of.txs.4.gene))=="4237")
  #
  # GRList.data.frame<-lapply(GRList,as.data.frame)
  #
  # GRList.data.frame[[5922]]

  data2<-txs
  data3<-cbind(data2,unlist(data2$SYMBOL))
  colnames(data3)[9]="Gene"

  txs.num<-count(data3,"Gene")
  hist(txs.num)

  txs.min<-do.call(rbind, by(data3, data3$Gene, function(x) x[which.min(x$width),]))
  txs.max<-do.call(rbind, by(data3, data3$Gene, function(x) x[which.max(x$width),]))

  txs.min.start<-do.call(rbind, by(data3, data3$Gene, function(x) x[which.min(x$start),]))
  txs.max.end<-do.call(rbind, by(data3, data3$Gene, function(x) x[which.max(x$end),]))

  temp<-merge(txs.min.start,txs.max.end,by="Gene",suffixes = c(".min.start",".max.end"))
  temp2<-merge(temp,txs.min,by="Gene")
  temp3<-merge(temp2,txs.max,by="Gene",suffixes = c(".shortest",".longest"))
  temp4<-merge(temp3,txs.num,by="Gene")

  #temp3<-temp[order(temp3$Gene),]


  txs.num.4.genes.interested<-txs.num[which(txs.num$Gene %in% genes.interested[,1]),]
  txs.min.4.genes.interested<-txs.min[which(txs.min$Gene %in% genes.interested[,1]),]
  txs.max.4.genes.interested<-txs.max[which(txs.max$Gene %in% genes.interested[,1]),]

  re<-list(txs=txs,txs.num=txs.num,txs.min=txs.min,txs.max=txs.max,txs.min.start.max.end=temp4,
           txs.4.genes.interested=txs.4.genes.interested,
           txs.num.4.genes.interested=txs.num.4.genes.interested,
           txs.min.4.genes.interested=txs.min.4.genes.interested,
           txs.max.4.genes.interested=txs.max.4.genes.interested)

  return(re)

}

# index<-which(Re.txs$txs.min.start.max.end$seqnames.min.start!=Re.txs$txs.min.start.max.end$seqnames.max.end)
#
# dim(Re.txs$txs.min.start.max.end[index,])
# names(Re.txs$txs.min.start.max.end)
#
# head(Re.txs$txs.min.start.max.end[index,c(1,2,10)],275)

#data2<-head(Re.txs$txs.all,100)
#data2

#do.call(rbind, by(data2, unlist(data2$SYMBOL), function(x) x[which.min(x$width),]))

#do.call(rbind,tapply(1:nrow(data2), unlist(data2$SYMBOL), function(x) data2[x,][which.min(data2$width[x]),]))
#dara

