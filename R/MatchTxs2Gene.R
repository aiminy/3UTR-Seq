#' MatchTxs2Gene
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' dir.name="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"
#' output.file="InterestedGeneTxsSorted.csv"
#'
#'
#' Re.interested.gene<-MatchTxs2Gene.R(re.rMAT,Re.txs,dir.name,output.file)
#'
#'
#' Re.interested.gene[which(Re.interested.gene$SYMBOL %in% c("XYLB")),]
#'
#'
#'
MatchTxs2Gene.R <- function(re.rMAT,Re.txs,dir.name,out.file) {
  txs.DE<-re.rMAT$DE
  txs.all<-Re.txs$txs
  head(txs.DE)
  head(txs.all)
  rownames(txs.all)<-txs.all$TXNAME
  txs.DE.all<-merge(txs.DE,txs.all,by=0)

  txs.4.certain.gene<-txs.DE.all[which(txs.DE.all$SYMBOL %in% Re.txs$txs.4.genes.interested$SYMBOL),]

  gene<-unlist2(txs.4.certain.gene$SYMBOL)

  txs.4.certain.gene.2<-cbind(txs.4.certain.gene,gene)

  txs.4.certain.gene.3<-txs.4.certain.gene.2[order(txs.4.certain.gene.2$gene),]

  write.csv(txs.4.certain.gene.3,file=paste0(dir.name,out.file))

  return(txs.4.certain.gene.2)

}
