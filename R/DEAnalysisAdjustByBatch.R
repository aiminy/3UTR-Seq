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
#' rownames(df.NT)<-df.NT[,1]
#' df.NT.2<- df.NT[,-1]
#'
#' re.DoGs.adjusted.by.batch<-DEAnalysisAdjustByBatch(df.NT.2)
#'

DEAnalysisAdjustByBatch <- function(df.NT) {

  cell<-factor(rep(c('emp','hela'),c(2,2)))
  cell=rep(cell,2)

  colData <- data.frame(condition=factor(rep(c('Dox', 'WT'),c(4, 4))),cell=cell)

  dds <- DESeqDataSetFromMatrix(df.NT, colData, formula(~condition))
  dat <- counts(dds)

  mod <- model.matrix(~condition,colData(dds))
  mod0 <- model.matrix(~1,colData(dds))
  svseq.NT <- svaseq(dat, mod, mod0)

  #Adjust by surrogate variables 1 and 2
  ddssva <- dds
  ddssva$SV1 <- svseq.NT$sv[,1]
  ddssva$SV2 <- svseq.NT$sv[,2]
  design(ddssva) <- ~ SV1 + SV2 + condition
  ddssva <- DESeq(ddssva)
  re.DESeq <- results(ddssva)

  re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
  colnames(re.FC)[7]="FoldChange(WT/Dox)"

  txs.gene<-ReformatTxsGene()

  re.FC<-merge(re.FC,txs.gene$txs_genes_DF_2,by=0)

  re.FC.sorted<-re.FC[order(re.FC$pvalue),]

  return(re.FC.sorted)

}
