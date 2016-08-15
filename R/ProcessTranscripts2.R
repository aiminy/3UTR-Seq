#' ProcessTranscript2
#'
#' @return
#' @export
#'
#' @examples
#'
#' ProcessTranscript2()
ProcessTranscript2 <- function() {

  #Gene based
  genes_hg19<-as.data.frame(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

  #Transcripts based
  txs.genes<-transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)

  gene.symbol<-getSYMBOL(names(txs.genes), data='org.Hs.eg')

  txs.genes.2<-txs.genes
  names(txs.genes.2)<-gene.symbol

  gene.symbol.2<-getSYMBOL(genes_hg19$gene_id,data = 'org.Hs.eg')
  gene.hg19.gene.symbol<-cbind(genes_hg19,gene.symbol.2)

  txs.genes.3<-lapply(txs.genes.2,as.data.frame)

  txs.genes.3.name<-sapply(txs.genes.3,function(x){

    gene<-names(x)
    #
    #
    #})

    #print(x)

    n.txs<-dim(x)[1]
    y<-cbind(x,rep(gene,n.txs))
    y
    },simplify=FALSE,USE.NAMES = TRUE)

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

  #txs.genes.3.name.22<-unlist2(txs.genes.3.name.2)

  txs.genes.33<-do.call(rbind.data.frame,txs.genes.3)
  txs.genes.33.sorted <- txs.genes.33[order(row.names(txs.genes.33)), ]

  #txs.genes.33.sorted.2<-cbind(txs.genes.33.sorted,strsplit(rownames(txs.genes.33.sorted),split="."))

  #txs.gene.333<-data.frame(Reduce(rbind,txs.genes.3))

  print(txs.genes.3[which(names(txs.genes.3)=="STRA6")])

  txs.genes.4<-lapply(txs.genes.3,function(x){
    y<-unique(x$seqnames)
    ny<-length(y)
    ny})
  txs.genes.5<-do.call("rbind", lapply(txs.genes.4, "[[", 1))

  head(txs.genes.5)

  gene.not.unique.chr<-unique(names(txs.genes.5[which(txs.genes.5[,1]!=1),]))

  gene.unique.chr<-unique(names(txs.genes.5[which(txs.genes.5[,1]==1),]))


  txs.genes.3.name.not.unique<-txs.genes.3.name.2[which(names(txs.genes.3.name.2) %in% gene.not.unique.chr)]

  txs.genes.3.name.unique<-txs.genes.3.name.2[which(names(txs.genes.3.name.2) %in% gene.unique.chr)]


  txs.genes.4.strand<-lapply(txs.genes.3.name.unique,function(x){
    y<-unique(x$strand)
    ny<-length(y)
    ny})
  txs.genes.5.strand<-do.call("rbind", lapply(txs.genes.4.strand, "[[", 1))

  temp<-names(txs.genes.5.strand[which(txs.genes.5.strand[,1]==2),])

  txs.genes.3.unique.name.not.strand<-txs.genes.3.name.unique[which(names(txs.genes.3.name.unique) %in% temp)]

  length(unique(names(txs.genes.3.unique.name.not.strand)))

  temp2<-names(txs.genes.5.strand[which(txs.genes.5.strand[,1]==1),])

  txs.genes.3.unique.name.unique.strand<-txs.genes.3.name.unique[which(names(txs.genes.3.name.unique) %in% temp2)]

  length(unique(names(txs.genes.3.unique.name.unique.strand)))


  input.dir="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/"
  input.file.pattern="final_list.csv"
  file.name=paste0(input.dir,dir(input.dir,recursive = TRUE,pattern=input.file.pattern))
  genes.interested<-read.csv(file.name,header=F)

  length(which(names(txs.genes.3.unique.name.unique.strand) %in% genes.interested[,1]))

  length(which(names(txs.genes.3.unique.name.not.strand) %in% genes.interested[,1]))

  length(which(names(txs.genes.3.name.not.unique) %in% genes.interested[,1]))

  txs.genes.3.name.not.unique[which(names(txs.genes.3.name.not.unique) %in% genes.interested[,1])]

  print(gene.hg19.gene.symbol[which(gene.hg19.gene.symbol$gene.symbol.2=="STRA6"),])

  re<-list(txs_genes_list=txs.genes.3.name.2,txs_genes_DF=txs.genes.3.name.2.DF)

  return(re)

}
