DEAnalysis <- function(countData,wt.index,dox.index)
{

  a <- length(wt.index)
  b <- length(dox.index)

  colData <- data.frame(condition = factor(c(rep("Dox",
                                                 b), rep("WT", a))))

  dds <- DESeqDataSetFromMatrix(countData, colData, formula(~condition))

  # size.factor<-estimateSizeFactors(dds)

  # print(size.factor)

  re.DESeq <- results(DESeq(dds))

  re.FC <- cbind(as.data.frame(re.DESeq), 2^re.DESeq[,
                                                     2], counts(dds))
  colnames(re.FC)[7] = "FoldChange(WT/Dox)"

  txs.gene <- ReformatTxsGene()

  re.FC <- merge(re.FC, txs.gene$txs_genes_DF_2, by = 0)

  re.FC.sorted <- re.FC[order(re.FC$pvalue), ]

  return(re.FC.sorted)

}

#' getutrcount
#'
#' @param dir.name
#' @param input.file.pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name <- "/Volumes/Bioinformatics$/Aimin_project/UTR/NewCounts"
#' input.file.pattern <- "count.txt"
#'
#' res <- ThreeUTR:::getutrcount(dir.name, input.file.pattern,file.path(system.file("extdata",package = "ThreeUTR"),"sample_infor.txt"))
#'
getutrcount <- function(dir.name, input.file.pattern,sample.infor.file)
{
  file.name = file.path(dir.name, dir(dir.name, recursive = TRUE,
                                      pattern = input.file.pattern))
  file.name.2 <- as.list(file.name)

  names(file.name.2) = basename(file.name)


  file.name.3 <- lapply(file.name.2, function(u)
  {
    if (!file.size(u) == 0){
      re = u
    }else{
      re= NULL
    }
    re
  })

  filterByRmNull <- function(re.peaks.only.bed) {
    re.peaks.only.bed.2<-re.peaks.only.bed[lapply(re.peaks.only.bed,length) > 0]

    # names(re.peaks.only.bed.2)=unlist(lapply(1:length(re.peaks.only.bed.2),function(u,re.peaks.only.bed.2){
    #   tmp=re.peaks.only.bed.2
    #   x=tmp[[u]]
    #   path_name=dirname(x)
    #   file_name=basename(x)
    #   file_name
    # },re.peaks.only.bed.2))
    return(re.peaks.only.bed.2)
  }

  file.name.4 <- filterByRmNull(file.name.3)

  sample.infor <- read.table(sample.infor.file,header = TRUE)

  names(file.name.4)=unlist(lapply(1:length(file.name.4),function(u,file.name.4,sample.infor){

    tmp=file.name.4
    x=tmp[[u]]
    path_name=dirname(x)
    file_name=basename(x)
    pos = gregexpr('-', file_name)
    a = pos[[1]][1]+1
    b = pos[[1]][2]-1
    c = substr(file_name,a,b)
    con = sample.infor[match(c,sample.infor$Sample),]$Condition
    file_name_con = paste0(con,"-",file_name)

    file_name_con


  },file.name.4,sample.infor))

  re.out <- lapply(file.name.4, function(u)
  {
    if (!file.size(u) == 0){
      re = read.table(u, header = F)
      colnames(re) = c("Count", "GeneName")
      re
    }
  })

  names(re.out)=names(file.name.4)

  temp.name <- strsplit(names(file.name.4), split = "\\.")
  temp.name.2 <- trimws(do.call("rbind", lapply(temp.name,
                                                "[[", 1)))

  temp.name.3 <- unique(temp.name.2[, 1])

  temp.name.4 <- as.list(temp.name.3)
  names(temp.name.4) <- temp.name.3

  txs.name <- lapply(re.out, "[[", 2)
  txs.name <- unlist(txs.name)
  txs.name <- unique(as.character(txs.name))
  txs.name.count <- data.frame(Count = rep(0, length(txs.name)),
                               GeneName = txs.name)

  # dim(txs.name.count)

  Process4EachSample <- function(one_sample, txs.name.count,
                                 re.out)
  {

    count.gene.matched <- re.out[(grep(one_sample, names(re.out)))]

    count.gene.matched.minus.gene <- count.gene.matched[(grep("minus.gene",
                                                              names(count.gene.matched)))]
    count.gene.matched.plus.gene <- count.gene.matched[(grep("plus.gene",
                                                             names(count.gene.matched)))]

    ReformatCount <- function(count.gene.matched.minus.gene,
                              txs.name.count)
    {

      reformat.count.gene.matched <- lapply(count.gene.matched.minus.gene,
                                            function(x, txs.name.count)
                                            {
                                              y <- txs.name.count
                                              if (dim(x[which(x$GeneName %in% y$GeneName),
                                                        ])[1] != 0)
                                              {
                                                y[match(x$GeneName, y$GeneName), ]$Count <- x$Count
                                              }
                                              y
                                            }, txs.name.count)


      reformat.count.gene.matched.2 <- cbind(apply(do.call("cbind",
                                                           lapply(reformat.count.gene.matched, "[", 2)),
                                                   1, unique), do.call("cbind", lapply(reformat.count.gene.matched,
                                                                                       "[", 1)))

      colnames(reformat.count.gene.matched.2) = c("GeneName",
                                                  names(count.gene.matched.minus.gene))

      gene <- as.character(reformat.count.gene.matched.2$GeneName)

      count.DoGs.plus.read <- apply(as.data.frame(reformat.count.gene.matched.2[,
                                                                                c(grep("plus.read.DoGs.count.txt", colnames(reformat.count.gene.matched.2)))]),
                                    1, sum)
      count.DoGs.minus.read <- apply(as.data.frame(reformat.count.gene.matched.2[,
                                                                                 c(grep("minus.read.DoGs.count.txt", colnames(reformat.count.gene.matched.2)))]),
                                     1, sum)


      reformat.count.gene.matched.3 <- as.data.frame(cbind(gene,
                                                           count.DoGs.plus.read, count.DoGs.minus.read))

      return(reformat.count.gene.matched.3)
    }

    re1 <- ReformatCount(count.gene.matched.minus.gene, txs.name.count)
    re2 <- ReformatCount(count.gene.matched.plus.gene, txs.name.count)
    re <- rbind(re1, re2)

    re <- re[order(re$gene), ]

    return(re)
  }

  re.8.samples <- lapply(temp.name.4, function(u, txs.name.count,
                                               re.out)
  {
    YY <- Process4EachSample(u, txs.name.count, re.out)
    YY
  }, txs.name.count, re.out)

  ProcessEachCorner <- function(re.8.samples, i,sample.infor)
  {

    gene <- apply(do.call("cbind", lapply(re.8.samples, "[",
                                          1)), 1, unique)
    count.gene.plus.read.8.samples <- cbind(gene, do.call("cbind",
                                                          lapply(re.8.samples, "[", i)))

    colnames(count.gene.plus.read.8.samples) = c("gene",
                                                 names(re.8.samples))

    df <- count.gene.plus.read.8.samples
    # df <- apply(df[, -1], 2, as.numeric)
    # index <- rowSums(df[, -1]) > 0
    # dff <- count.gene.plus.read.8.samples[index, ]
    # rownames(dff) <- dff$gene
    # dff <- dff[, -1]
    #
    # X <- toupper(unique(as.character(sample.infor$Condition)))
    #
    # wt.index <- grep(X[1],toupper(colnames(dff)))
    #
    # dox.index <- grep(X[2],toupper(colnames(dff)))
    #
    # real.index <- c(dox.index,wt.index)
    #
    # n.sample <- length(real.index)
    #
    # permutation.index <- real.index
    #
    # countData <- apply(dff[, permutation.index], 2, as.numeric)
    #
    # rownames(countData) <- rownames(dff)
    #
    # # txs.gene <- ReformatTxsGene()
    # #
    # # xx <- txs.gene$txs_genes_DF_2
    # #
    # # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    # #
    # # results <- getBM(attributes = c("ucsc","ensembl_transcript_id"),
    # #                  filters = "ensembl_transcript_id", values = xx,
    # #                  mart = mart)
    # #
    # # xxx <- results[-which(results$ucsc==""),]
    # #
    # # rownames(xxx) <- xxx$ucsc
    # #
    # # re.FC <- merge(countData,xxx, by = 0)
    #
    # re.FC <- countData

    #return(re.FC)
    return( df )

  }

   re.BR <- ProcessEachCorner(re.8.samples, 2)  #BR
   re.TR <- ProcessEachCorner(re.8.samples, 3)  #TR


  #  re.BR.4.plus.gene.BL.4.minus.gene <- re.BR
  #  re.TR.4.plus.gene.TL.4.minus.gene <- re.TR
  #
  # # Get the counts for DoGs of plus and minus gene
  # DoGs.4.plus.Gene <- re.TR.4.plus.gene.TL.4.minus.gene[which(re.TR.4.plus.gene.TL.4.minus.gene$strand ==
  #                                                               "+"), ]
  # DoGs.4.minus.Gene <- re.BR.4.plus.gene.BL.4.minus.gene[which(re.BR.4.plus.gene.BL.4.minus.gene$strand ==
  #                                                                "-"), ]
  #
  # DoGs.4.plus.minus.Gene <- rbind(DoGs.4.plus.Gene, DoGs.4.minus.Gene)
  #
  # res <- list(DoGs.4.plus.Gene=DoGs.4.minus.Gene,
  #             DoGs.4.minus.Gene=DoGs.4.minus.Gene,
  #             DoGs.4.plus.minus.Gene=DoGs.4.plus.minus.Gene)

   res <- list(reBR=re.BR,reTR=re.TR,re.8.samples=re.8.samples)


  return(res)

}

#' Title
#'
#' @param res
#'
#' @return
#' @export
#'
#' @examples
#'
#' y <- ThreeUTR:::annotatetranscript(res$reBR)
#'
annotatetranscript <- function(res)
{

  x <- res
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "hsapiens_gene_ensembl",
                           host = 'ensembl.org')

  t2g <-  getBM(attributes = c("ucsc","ensembl_transcript_id","strand","ensembl_gene_id",
                                       "external_gene_name"),  filters = "ensembl_transcript_id",values = rownames(x),mart = mart)

  #xx <- t2g[-which(t2g$ensembl_transcript_id==""),]

  #rownames(xx) <- xx$ensembl_transcript_id

  #xxx <- merge(x,xx,by=0)

  xxx <- t2g

  xxx
}

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
#' dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/'
#'
#' input.file.pattern='*count.2.txt'
#'
#' dir.name.gene.list='/media/H_driver/2016/Ramin_azhang/'
#' pattern.4.gene.list='final_list.csv'
#'
#'
#' out.dir.name='/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOne/'
#'
#' for (i in 1:10) {
#'
#' time.string = gsub(':', '-', gsub(' ', '_', as.character(Sys.time())))
#' out.dir.name=paste0('/media/aiminyan/DATA/Ramin_azhang/Counts4DoGsOnlyRmOnePermutationAt_',time.string,'/')
#' dir.create(out.dir.name)
#'
#' out.file.pattern.interested='DoGs_adjust_by_batch_interested_gene'
#' out.file.pattern.positive.gene='DoGs_adjust_by_batch_positive'
#' out.file.pattern.negative.gene='DoGs_adjust_by_batch_negative'
#' out.file.pattern.all= 'DoGs_adjust_by_batch_all'
#'
#' Re.unadjusted.adjusted<-ProcessOutputFilesFromDoGsOnly(dir.name,input.file.pattern,out.dir.name,out.file.pattern.interested,
#' out.file.pattern.positive.gene,
#' out.file.pattern.negative.gene,
#' out.file.pattern.all,
#' dir.name.gene.list,
#' pattern.4.gene.list,
#' adjust_by_batch='YES')
#'
#' }
#'
#' save.image(file=paste0(out.dir.name,'re_save_2.RData'))
#' savehistory(file=paste0(out.dir.name,'re_save_2.Rhistory'))

ProcessOutputFilesFromDoGsOnly <- function(dir.name, input.file.pattern,
    out.dir.name, out.file.pattern.interested, out.file.pattern.positive.gene,
    out.file.pattern.negative.gene, out.file.pattern.all, dir.name.gene.list,
    pattern.4.gene.list, adjust_by_batch,permutation.set.up=NULL)
    {

    x <- getutrcount(dir.name, input.file.pattern)

    DoGs.4.plus.Gene <- x[[1]]
    DoGs.4.minus.Gene <- x[[2]]
    DoGs.4.plus.minus.Gene <- x[[3]]

    GeneTypeBasedDE <- function(DoGs.4.plus.Gene)
    {

        #a <- n+9-1

        Count.DoGs.4.plus.Gene <- DoGs.4.plus.Gene
        rownames(Count.DoGs.4.plus.Gene) <- Count.DoGs.4.plus.Gene$Row.names
        Count.DoGs.4.plus.Gene.2 <- Count.DoGs.4.plus.Gene[,
            -1]

        # head(Count.DoGs.4.plus.Gene.2)

        wt.index <- grep("WT",toupper(colnames(Count.DoGs.4.plus.Gene.2)))

        dox.index <- grep("DOX",toupper(colnames(Count.DoGs.4.plus.Gene.2)))

        real.index <- c(dox.index,wt.index)

        n <- length(real.index)

        if(!is.null(permutation.set.up)){
        real.index <- seq(1,n)
        permutation.index <- real.index
        permutation.index = array(sample(real.index))
        }else
        {
          permutation.index <- real.index
        }

        tmp <- Count.DoGs.4.plus.Gene.2

        Count.DoGs.4.plus.Gene.2 <- tmp[, permutation.index]

        cat("get count\n")
        print(head(Count.DoGs.4.plus.Gene.2))
        cat("get count done \n")

        if (adjust_by_batch == "NO")
        {
            re.DESeq.DoGs.plus.gene <- DEAnalysis(Count.DoGs.4.plus.Gene.2,wt.index,dox.index)
        } else
        {

          re.DESeq.DoGs.plus.gene <- DEAnalysisAdjustByBatch(Count.DoGs.4.plus.Gene.2,Nbatch=c(3,2))
        }

        return(re.DESeq.DoGs.plus.gene)
    }

    print(head(DoGs.4.plus.Gene))
    re.DESeq.DoGs.plus.gene <- GeneTypeBasedDE(DoGs.4.plus.Gene)
    # head(re.DESeq.DoGs.plus.gene[which(re.DESeq.DoGs.plus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.gene$pvalue<0.05),])

    print(head(DoGs.4.minus.Gene))
    re.DESeq.DoGs.minus.gene <- GeneTypeBasedDE(DoGs.4.minus.Gene)
    # head(re.DESeq.DoGs.minus.gene[which(re.DESeq.DoGs.minus.gene$log2FoldChange<0&re.DESeq.DoGs.minus.gene$pvalue<0.05),])

    print(head(DoGs.4.plus.minus.Gene))
    re.DESeq.DoGs.plus.minus.gene <- GeneTypeBasedDE(DoGs.4.plus.minus.Gene,n)
    # head(re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$log2FoldChange<0&re.DESeq.DoGs.plus.minus.gene$pvalue<0.05),])

    gene.interested <- ReadGeneList(dir.name.gene.list, pattern.4.gene.list)

    re.DESeq.DoGs.plus.minus.gene.interested <- re.DESeq.DoGs.plus.minus.gene[which(re.DESeq.DoGs.plus.minus.gene$gene %in%
        gene.interested[, 1]), ]

    # re.DESeq.DoGs.plus.minus.gene.interested[which(re.DESeq.DoGs.plus.minus.gene.interested$gene
    # %in% c('TPCN2')),]
    # re.FC.sorted<-re.FC[order(re.FC$pvalue),]

    ReformatGeneSymbol <- function(re.DoGs.adjusted.by.batch)
    {
        re2 <- re.DoGs.adjusted.by.batch
        re2$gene <- as.character(re2$gene)
        re2$gene <- paste0("'", as.character(re2$gene), "'")
        return(re2)
    }

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.positive.gene,
        ".csv")), row.names = FALSE, quote = FALSE)

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.minus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.negative.gene,
        ".csv")), row.names = FALSE, quote = FALSE)

    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.all,
        ".csv")), row.names = FALSE, quote = FALSE)
    # Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)
    re2 <- ReformatGeneSymbol(re.DESeq.DoGs.plus.minus.gene.interested)
    write.csv(re2, file = file.path(out.dir.name, paste0("3UTR_DE_", out.file.pattern.interested,
        ".csv")))

    # write.csv(re.DESeq.DoGs.plus.minus.gene,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE) Re.unadjusted.adjusted$DE

    # write.csv(Re.unadjusted.adjusted$DE,file=paste0(out.dir.name,'3UTR_DE_',out.file.pattern.all,'.csv'),row.names
    # = FALSE,quote=FALSE)


    re.out.3 <- list(dogcount=x,
        DE_positive_gene = re.DESeq.DoGs.plus.gene, DE_negative_gene = re.DESeq.DoGs.minus.gene,
        DE = re.DESeq.DoGs.plus.minus.gene, DE_interested = re.DESeq.DoGs.plus.minus.gene.interested)

    return(re.out.3)

}




idtransform <- function(values_to_be_transforemd) {

  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl",
                  host = 'ensembl.org')

  t2g <-  getBM(attributes = c("ucsc","ensembl_transcript_id","strand","ensembl_gene_id",
                               "external_gene_name"),  filters = "ucsc",values = values_to_be_transforemd,mart = mart)
  t2g


}

#' Title
#'
#' @param input.file
#' @param output.file
#'
#' @return
#' @export
#'
#' @examples
#' input.file <- "/Volumes/Bioinformatics$/Zhen-Gao/UTR/Gencode-3UTR-transcript.tsv"
#' output.file <- "/Volumes/Bioinformatics$/Zhen-Gao/UTR/Gencode-3UTR-transcript-transformed.tsv"
#' ThreeUTR:::gaoidtransform(input.file,output.file)
#'
gaoidtransform <- function(input.file, output.file) {

  res <- read.table(input.file,header=F)
  colnames(res) <- c("Chr","Start","End","utr","Symbol","Strand","ucsc")

  #head(res)

  y <- ThreeUTR:::idtransform(res$ucsc)

  yy <- merge(res,y,by="ucsc")
  write.table(yy, file = output.file, append = FALSE, quote = FALSE, sep = "\t",
                          eol = "\n", na = "NA", dec = ".", row.names = FALSE)
}
