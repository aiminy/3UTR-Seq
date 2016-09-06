#' Title
#'
#' @param Re.unadjusted.adjusted
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' out.dir.name="/media/aiminyan/DATA/Ramin_azhang/Counts5CasesEachSample/"
#' out.file.pattern="DoGs_using_longest_txs_adjust_by_batch"
#'
#' re.longest.transcript<-TestUseLongestTxsOfGene(Re.unadjusted.adjusted,out.dir.name,out.file.pattern)
#'
#'
TestUseLongestTxsOfGene<-function(Re.unadjusted.adjusted,out.dir.name,out.file.pattern){

  tempData<-Re.unadjusted.adjusted$DE

  #which(tempData$gene==NA)

  tempdata.byGSym.2<-tempData

  rownames(tempdata.byGSym.2) = NULL

  data.byGSym = ddply(tempdata.byGSym.2, c("gene"),function(h){
    #summary = apply(h,2,max)
    y<-as.data.frame(h[which.max(h$width),])
    y
   }
  )

  df.NT<-data.byGSym[,c(1,9:16)]
  rownames(df.NT)<-df.NT[,1]
  df.NT.2<- df.NT[,-1]
  #'
  re.DoGs.adjusted.by.batch<-DEAnalysisAdjustByBatch(df.NT.2)


  re.DoGs.adjusted.by.batch[which(re.DoGs.adjusted.by.batch$Row.names=="uc010paz.2"),]

  #length(unique(as.character(re.DoGs.adjusted.by.batch$gene)))
  #unique(as.character(re.DoGs.adjusted.by.batch$gene))
  #[grep("uc002jah.2",unique(as.character(re.DoGs.adjusted.by.batch$Row.names)))]

  write.csv(re.DoGs.adjusted.by.batch,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern,".csv")
            ,row.names = FALSE,quote=FALSE)

  return(re.DoGs.adjusted.by.batch)

}
