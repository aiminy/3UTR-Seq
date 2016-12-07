#' ProcessOutputFilesFromDoGsOnly
#'
#' @param dir.name
#' @param input.file.pattern
#' @param out.dir.name
#' @param out.file.pattern.interested
#' @param out.file.pattern.positive.gene
#' @param out.file.pattern.negative.gene
#' @param out.file.pattern.all
#' @param dir.name.gene.list
#' @param pattern.4.gene.list
#' @param adjust_by_batch
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/"
#'
#' input.file.pattern="*count.2.txt"
#'
#' dir.name.gene.list="/media/H_driver/2016/Ramin_azhang/"
#' pattern.4.gene.list="final_list.csv"
#'
#'
#' out.dir.name="/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/"
#'
#' for (i in 1:10) {
#'
#' time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
#' out.dir.name=paste0("/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOnePermutationAt_",time.string,"/")
#' dir.create(out.dir.name)
#'
#' out.file.pattern.interested="DoGs_adjust_by_batch_interested_gene"
#' out.file.pattern.positive.gene="DoGs_adjust_by_batch_positive"
#' out.file.pattern.negative.gene="DoGs_adjust_by_batch_negative"
#' out.file.pattern.all= "DoGs_adjust_by_batch_all"
#'
#' Re.unadjusted.adjusted<-ProcessOutputFilesFromDoGsOnly(dir.name,input.file.pattern,out.dir.name,out.file.pattern.interested,
#' out.file.pattern.positive.gene,
#' out.file.pattern.negative.gene,
#' out.file.pattern.all,
#' dir.name.gene.list,
#' pattern.4.gene.list,
#' adjust_by_batch="YES")
#'
#' }
#'
#' save.image(file=paste0(out.dir.name,"re_save_2.RData"))
#' savehistory(file=paste0(out.dir.name,"re_save_2.Rhistory"))

ProcessOutputFilesFromDoGsOnly<-function(dir.name,
                                         input.file.pattern,
                                         out.dir.name,
                                         out.file.pattern.interested,
                                         out.file.pattern.positive.gene,
                                       out.file.pattern.negative.gene,
                                       out.file.pattern.all,
                                       dir.name.gene.list,
                                       pattern.4.gene.list,
                                       adjust_by_batch){
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",7)

  print(file.name.2)

  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=F)
    colnames(re)=c("Count","GeneName")
    re
  })

  temp.name<-strsplit(names(file.name.2),split="\\.")

  temp.name.2<-trimws(do.call("rbind",lapply(temp.name,"[[",1)))

  temp.name.3<-unique(temp.name.2[,1])

  temp.name.4<-as.list(temp.name.3)
  names(temp.name.4)<-temp.name.3

  txs.name<-lapply(re.out,"[[",2)
  txs.name<-unlist(txs.name)
  txs.name<-unique(as.character(txs.name))
  txs.name.count<-data.frame(Count=rep(0,length(txs.name)),GeneName=txs.name)

  #dim(txs.name.count)

    Process4EachSample <- function(one_sample,txs.name.count,re.out) {

      count.gene.matched<-re.out[(grep(one_sample,names(re.out)))]

      count.gene.matched.minus.gene<-count.gene.matched[(grep("minus.gene",names(count.gene.matched)))]
      count.gene.matched.plus.gene<-count.gene.matched[(grep("plus.gene",names(count.gene.matched)))]

      ReformatCount <- function(count.gene.matched.minus.gene,txs.name.count) {

        reformat.count.gene.matched<-lapply(count.gene.matched.minus.gene,function(x,txs.name.count){
          y<-txs.name.count
          if(dim(x[which(x$GeneName %in% y$GeneName),])[1]!=0)
          {
            y[match(x$GeneName,y$GeneName),]$Count<-x$Count
          }
          y
          },txs.name.count)


        reformat.count.gene.matched.2<-cbind(apply(do.call("cbind",lapply(reformat.count.gene.matched,"[",2)),1,unique),
        do.call("cbind",lapply(reformat.count.gene.matched,"[",1)))

        colnames(reformat.count.gene.matched.2)=c("GeneName",names(count.gene.matched.minus.gene))

        gene<-as.character(reformat.count.gene.matched.2$GeneName)

        count.DoGs.plus.read<-apply(as.data.frame(reformat.count.gene.matched.2[,c(grep("plus.read.DoGs.count.2.txt",colnames(reformat.count.gene.matched.2)))]),1,sum)
        count.DoGs.minus.read<-apply(as.data.frame(reformat.count.gene.matched.2[,c(grep("minus.read.DoGs.count.2.txt",colnames(reformat.count.gene.matched.2)))]),1,sum)


        reformat.count.gene.matched.3<-as.data.frame(cbind(gene,count.DoGs.plus.read,count.DoGs.minus.read))

        return(reformat.count.gene.matched.3)
      }

      re1<-ReformatCount(count.gene.matched.minus.gene,txs.name.count)
      re2<-ReformatCount(count.gene.matched.plus.gene,txs.name.count)
      re<-rbind(re1,re2)

      re<-re[order(re$gene),]

      return(re)
    }

   re.8.samples<-lapply(temp.name.4,function(u,txs.name.count,re.out){
     YY<-Process4EachSample(u,txs.name.count,re.out)
     YY},txs.name.count,re.out)

   DEAnalysis <- function(countData) {
     colData <- data.frame(condition=factor(c(rep("Dox",3),rep("WT",3))))

     dds <- DESeqDataSetFromMatrix(countData, colData, formula(~ condition))

     #size.factor<-estimateSizeFactors(dds)

     #print(size.factor)

     re.DESeq<-results(DESeq(dds))

     re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
     colnames(re.FC)[7]="FoldChange(WT/Dox)"

     txs.gene<-ReformatTxsGene()

     re.FC<-merge(re.FC,txs.gene$txs_genes_DF_2,by=0)

     re.FC.sorted<-re.FC[order(re.FC$pvalue),]

     return(re.FC.sorted)

   }

   ProcessEachCorner <- function(re.8.samples,i){

     gene<-apply(do.call("cbind",lapply(re.8.samples,"[",1)),1,unique)
     count.gene.plus.read.8.samples<-cbind(gene,do.call("cbind",lapply(re.8.samples,"[",i)))

     colnames(count.gene.plus.read.8.samples)=c("gene",names(re.8.samples))

     df<-count.gene.plus.read.8.samples
     df<-apply(df[,-1],2,as.numeric)
     index<-rowSums(df[, -1])>0
     dff<-count.gene.plus.read.8.samples[index,]
     rownames(dff)<-dff$gene
     dff<-dff[,-1]

     #print(head(dff))

    #countData <- apply(dff[,c(1,3,5,6,2,4,7,8)], 2, as.numeric)
    real.index<-c(1,3,4,2,5,6)

    permutation.index<-real.index
    #permutation.index=array(sample(real.index))

    countData <- apply(dff[,permutation.index], 2, as.numeric)

    rownames(countData)<-rownames(dff)

    re.FC<-DEAnalysis(countData)

    #re.DoGs.adjusted.by.batch<-DEAnalysisAdjustByBatch(countData)

    # re.FC<-cbind(as.data.frame(re.DESeq),2^re.DESeq[,2],counts(dds))
    # colnames(re.FC)[7]="FoldChange(WT/Dox)"
    #
    # txs.gene<-ReformatTxsGene()
    #
    #re.FC<-list(re.FC=re.FC,re.DoGs.adjusted.by.batch=re.DoGs.adjusted.by.batch)

    return(re.FC)

   }

  #re.BL<-ProcessEachCorner(re.8.samples,2) #BL
  #re.TL<-ProcessEachCorner(re.8.samples,3) #TL

  re.BR<-ProcessEachCorner(re.8.samples,2) #BR
  re.TR<-ProcessEachCorner(re.8.samples,3) #TR

  #head(re.BL)
  #head(re.BR)

  #head(re.TL)
  #head(re.TR)

  #re.TR[which(re.TR$Row.names =="uc003tzi.4"),]

  #Check for uc003tzi.4:gene(+)
  #re.BL[[1]][which(re.BL[[1]]$Row.names =="uc003tzi.4"),] #BL for gene(+)
  #re.TL[[1]][which(re.TL[[1]]$Row.names =="uc003tzi.4"),] #TL for gene(+) gene
  #re.BR[[1]][which(re.BR[[1]]$Row.names =="uc003tzi.4"),] #BR for gene(+)
  #re.TR[[1]][which(re.TR[[1]]$Row.names =="uc003tzi.4"),] #TR for gene(+) DoGs

  #Check for uc002szf.1:gene(-)
  #re.BL[[1]][which(re.BL[[1]]$Row.names =="uc002szf.1"),] #BR for gene(-)
  #re.TL[[1]][which(re.TL[[1]]$Row.names =="uc002szf.1"),] #
  #re.BR[[1]][which(re.BR[[1]]$Row.names =="uc002szf.1"),] #
  #re.TR[[1]][which(re.TR[[1]]$Row.names =="uc002szf.1"),] #

  #Check for uc001aac.4:gene(-)
  #re.BL[[1]][which(re.BL[[1]]$Row.names =="uc001aac.4"),]#BR for gene(-) Gene
  #re.TL[[1]][which(re.TL[[1]]$Row.names =="uc001aac.4"),]#TR for gene(-)
  #re.TR[[1]][which(re.BR[[1]]$Row.names =="uc001aac.4"),]#TL for gene(-)
  #re.BR[[1]][which(re.TR[[1]]$Row.names =="uc001aac.4"),]#BL for gene(-) DoGs

  #if(adjust_by_batch=="NO"){
  #re.BL.4.plus.gene.BR.4.minus.gene<-re.BL
  #re.TL.4.plus.gene.TR.4.minus.gene<-re.TL

  re.BR.4.plus.gene.BL.4.minus.gene<-re.BR
  re.TR.4.plus.gene.TL.4.minus.gene<-re.TR
#}else{
  #  re.BL.4.plus.gene.BR.4.minus.gene<-re.BL[[2]]
  #  re.TL.4.plus.gene.TR.4.minus.gene<-re.TL[[2]]
  #  re.BR.4.plus.gene.BL.4.minus.gene<-re.BR[[2]]
  #  re.TR.4.plus.gene.TL.4.minus.gene<-re.TR[[2]]
  #}

  #Get the counts for DoGs of plus and minus gene
  DoGs.4.plus.Gene<-re.TR.4.plus.gene.TL.4.minus.gene[which(re.TR.4.plus.gene.TL.4.minus.gene$strand=="+"),]
  DoGs.4.minus.Gene<-re.BR.4.plus.gene.BL.4.minus.gene[which(re.BR.4.plus.gene.BL.4.minus.gene$strand=="-"),]

  DoGs.4.plus.minus.Gene<-rbind(DoGs.4.plus.Gene, DoGs.4.minus.Gene)

  #DoGs.4.plus.Gene[which(DoGs.4.plus.Gene$Row.names =="uc003tzi.4"),] #TR for gene(+) DoGs
  #DoGs.4.minus.Gene[which( DoGs.4.minus.Gene$Row.names =="uc001aac.4"),]#BL for gene(-) DoGs

  GeneTypeBasedDE <- function(DoGs.4.plus.Gene) {
    Count.DoGs.4.plus.Gene<-DoGs.4.plus.Gene[,c(1,9:14)]
    rownames(Count.DoGs.4.plus.Gene)<-Count.DoGs.4.plus.Gene$Row.names
    Count.DoGs.4.plus.Gene.2<-Count.DoGs.4.plus.Gene[,-1]

    #head(Count.DoGs.4.plus.Gene.2)



    real.index<-c(1,2,3,4,5,6)
    permutation.index<-real.index
    permutation.index=array(sample(real.index))
    tmp<-Count.DoGs.4.plus.Gene.2

    Count.DoGs.4.plus.Gene.2<-tmp[,permutation.index]

    cat("get count\n")
    print(head(Count.DoGs.4.plus.Gene.2))
    cat("get count done \n")

    if(adjust_by_batch=="NO"){
    re.DESeq.DoGs.plus.gene<-DEAnalysis(Count.DoGs.4.plus.Gene.2)
    }else{
    re.DESeq.DoGs.plus.gene<-DEAnalysisAdjustByBatch(Count.DoGs.4.plus.Gene.2)
    }

    return(re.DESeq.DoGs.plus.gene)
  }

  print(head(DoGs.4.plus.Gene))
  re.DESeq.DoGs.plus.gene<-GeneTypeBasedDE(DoGs.4.plus.Gene)
  #head(re.DESeq.DoGs.plus.gene[which(re.DESeq.DoGs.plus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.gene$pvalue<0.05),])

  print(head(DoGs.4.minus.Gene))
  re.DESeq.DoGs.minus.gene<-GeneTypeBasedDE(DoGs.4.minus.Gene)
  #head(re.DESeq.DoGs.minus.gene[which(re.DESeq.DoGs.minus.gene$log2FoldChange<0&re.DESeq.DoGs.minus.gene$pvalue<0.05),])

  print(head(DoGs.4.plus.minus.Gene))
  re.DESeq.DoGs.plus.minus.gene<-GeneTypeBasedDE(DoGs.4.plus.minus.Gene)
  #head(re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.minus.gene$pvalue<0.05),])

  gene.interested<-ReadGeneList(dir.name.gene.list,pattern.4.gene.list)

  re.DESeq.DoGs.plus.minus.gene.interested<-re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$gene %in% gene.interested[,1]),]

  #re.DESeq.DoGs.plus.minus.gene.interested[which(re.DESeq.DoGs.plus.minus.gene.interested$gene %in% c("TPCN2")),]
  #re.FC.sorted<-re.FC[order(re.FC$pvalue),]

  ReformatGeneSymbol <- function(re.DoGs.adjusted.by.batch) {
    re2<-re.DoGs.adjusted.by.batch
    re2$gene<-as.character(re2$gene)
    re2$gene<-paste0("'",as.character(re2$gene),"'")
    return(re2)
  }

  re2<-ReformatGeneSymbol(re.DESeq.DoGs.plus.gene)
  write.csv(re2,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.positive.gene,".csv"),row.names = FALSE,quote=FALSE)

  re2<-ReformatGeneSymbol(re.DESeq.DoGs.minus.gene)
  write.csv(re2,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.negative.gene,".csv"),row.names = FALSE,quote=FALSE)

  re2<-ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene)
  write.csv(re2,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.all,".csv"),row.names = FALSE,quote=FALSE)
  #Re.unadjusted.adjusted$DE

  #write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.all,".csv"),row.names = FALSE,quote=FALSE)
  re2<-ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene.interested)
  write.csv(re2,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.interested,".csv"))

  #write.csv(re.DESeq.DoGs.plus.minus.gene,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.all,".csv"),row.names = FALSE,quote=FALSE)
  #Re.unadjusted.adjusted$DE

  #write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,"3UTR_DE_",out.file.pattern.all,".csv"),row.names = FALSE,quote=FALSE)


  re.out.3<-list(re.out=re.out,
                 re.8.samples=re.8.samples,
                 DE_positive_gene=re.DESeq.DoGs.plus.gene,
                 DE_negative_gene=re.DESeq.DoGs.minus.gene,
                 DE=re.DESeq.DoGs.plus.minus.gene,
                 DE_interested= re.DESeq.DoGs.plus.minus.gene.interested)

   return(re.out.3)

}
