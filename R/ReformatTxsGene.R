#' ReformatTxsGene
#'
#' @return
#' @export
#'
#' @examples
#'
#' ReformatTxsGene()
#'
ReformatTxsGene <- function() {

  #Transcripts based
  txs.genes<-transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
  gene.symbol<-getSYMBOL(names(txs.genes), data='org.Hs.eg')

  txs.genes.2<-txs.genes
  names(txs.genes.2)<-gene.symbol

  txs.genes.3<-lapply(txs.genes.2,as.data.frame)

  ListReformat <- function(txs.genes.3) {
    txs.genes.3.name<-lapply(seq_along(txs.genes.3), function(i){
      name.gene<-names(txs.genes.3)[[i]]
      n.txs<-dim(txs.genes.3[[i]])[1]
      gene<-rep(name.gene,n.txs)
      y<-cbind(txs.genes.3[[i]],gene)
      y
    })
    names(txs.genes.3.name)<-names(txs.genes.3)
    return(txs.genes.3.name)
  }

  txs.genes.3.name.2<-ListReformat(txs.genes.3)

  txs.genes.3.name.2.DF<-do.call(rbind.data.frame,txs.genes.3.name.2)

  txs.genes.3.name.2.DF.2<-txs.genes.3.name.2.DF

  rownames(txs.genes.3.name.2.DF.2)<-txs.genes.3.name.2.DF$tx_name

  re<-list(txs_genes_list=txs.genes.3.name.2,txs_genes_DF=txs.genes.3.name.2.DF,txs_genes_DF_2=txs.genes.3.name.2.DF.2)

  return(re)

}
