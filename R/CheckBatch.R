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
#' CheckBatch(dir.name,"Btach_using_DoGs",df.NT.2)
#'
CheckBatch <- function(dir.name,output_bacth_check,df.NT) {

  cell<-factor(rep(c('emp','hela'),c(2,2)))
  cell=rep(cell,2)

  colData.1 <- data.frame(condition=factor(rep(c('Dox', 'WT'),c(4, 4))),cell=cell)

  #Check hidden batch effects
  #dds <- estimateSizeFactors(dds)
  dds <- DESeqDataSetFromMatrix(df.NT, colData.1, formula(~condition))
  dat <- counts(dds)

  #idx <- rowMeans(dat) > 1
  #dat <- dat[idx,]

  mod <- model.matrix(~condition,colData(dds))
  mod0 <- model.matrix(~1,colData(dds))
  svseq.NT <- svaseq(dat, mod, mod0)

  png(paste0(dir.name,output_bacth_check,".png"))
  par(mfrow=c(2,1),mar=c(3,5,3,1))
  stripchart(svseq.NT$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
  abline(h=0)
  stripchart(svseq.NT$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
  abline(h=0)
  dev.off()

}
