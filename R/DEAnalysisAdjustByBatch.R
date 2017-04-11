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

DEAnalysisAdjustByBatch <- function(df.NT,Nbatch) {


  print(pkg.env$sample)

  cell<-factor(rep(c('emp','hela'),Nbatch))
  cell=rep(cell,2)

  colData <- data.frame(condition=factor(rep(c('Dox', 'WT'),c(5, 5))),cell=cell)

  dds <- DESeqDataSetFromMatrix(df.NT, colData, formula(~condition))
  dat <- counts(dds)

  mod <- model.matrix(~condition,colData(dds))
  mod0 <- model.matrix(~1,colData(dds))
  svseq.NT <- svaseq(dat, mod, mod0)

  #How many significant surrogate variables is?

  #ncol <- dim(svseq.NT$sv)[2]

  ncol <- svseq.NT$n.sv

  cat(ncol,"\n")

  #ncol<-2
  #Adjust by surrogate variables 1 and 2
  ddssva <- dds

  if(ncol==1){
    colData(ddssva)<-cbind2(colData(ddssva),svseq.NT$sv)
    colnames(colData(ddssva))[dim(colData(ddssva))[2]] = paste0("SV",1)
  }else{
    for(i in 1:ncol){
    colData(ddssva)<-cbind2(colData(ddssva),svseq.NT$sv[,i])
    colnames(colData(ddssva))[dim(colData(ddssva))[2]] = paste0("SV",i)
    }
  }


  #ddssva$SV1 <- svseq.NT$sv[,1]
  #ddssva$SV2 <- svseq.NT$sv[,2]

  #design(ddssva) <- ~ SV1 + SV2 + condition

  cat("ddssva")
  print(head(ddssva))
  cat("ddsva_done")

  design(ddssva)<-formula(paste("~ ",paste(paste0("SV", 1:ncol), collapse = "+"),"+", "condition" ))

  #formula(paste("~ ",paste(paste0("SV", 1:2), collapse = "+"),"+", "condition" ))

  ddssva <- DESeq(ddssva)
  re.DESeq <- results(ddssva)

  re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
  colnames(re.FC)[7]="FoldChange(WT/Dox)"

  txs.gene<-ReformatTxsGene()

  re.FC<-merge(re.FC,txs.gene$txs_genes_DF_2,by=0)

  re.FC.sorted<-re.FC[order(re.FC$pvalue),]

  return(re.FC.sorted)

}
