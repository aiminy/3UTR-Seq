#' Compare2Pos
#'
#' @param in1
#' @param in2
#'
#' @return
#' @export Compare2Pos
#'
#' @examples
#'
#' #real data
#' in1="/Volumes/Bioinformatics$/2016/Ramin_azhang/ResultsDoGsOnlyRmOne/RealData/3UTR_DE_DoGs_adjust_by_batch_all.csv"
#'
#' #permutated data
#'
#' in1="/Volumes/Bioinformatics$/2016/Ramin_azhang/ResultsDoGsOnlyRmOne/Counts4DoGsOnlyRmOnePermutation/3UTR_DE_DoGs_adjust_by_batch_all.csv"
#' threshold=0.05
#' in2="/Volumes/Bioinformatics$/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/final_list.csv"
#'
#' Compare2Pos(in1,threshold,in2)
#'
Compare2Pos<-function(in1,threshold,in2){

  a<-read.table(in1,header=T,sep=",")
  b<-read.table(in2,header=F,sep=",")

  print(head(a))

  print(head(b))

   g1<-unique(as.character(a[which(a$padj<=threshold),22]))
   g2<-unique(as.character(b[,1]))

   g11<-g1[-which(is.na(g1))]

   g22<-g2

   print(g11)

   print(g22)

  re<-list(R=g11,pos=g22,a=a,b=b)

  names(re)[1]<-paste0("R","_",threshold)

  venn.plot <- venn.diagram(
    x = re[c(1,2)],
    filename = paste0(dirname(in1),"/",names(re)[1],names(re)[2],"_overlap_venn.tiff"),
    height = 3000,
    width = 3500,
    resolution = 1000,
    col = "black",
    lty = "dotted",
    lwd = 1,
    fill = c("red","blue"),
    alpha = 0.50,
    label.col = c(rep("black",3)),
    cex = 0.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("red","blue"),
    cat.cex = 0.5,
    cat.pos = 0.5,
    cat.dist = 0.05,
    cat.fontfamily = "serif"
  )

  return(re)

}
