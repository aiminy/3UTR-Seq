#' @ ProcessOutputFilesFrom3UTR
#'
#' @ read output from 3UTR
#'
#' @param dir_name
#' @param input_file_pattern
#'
#' @return
#' @export
#'
#' @examples
#'
#' dir.name="/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/"
#' input.file.pattern="*downstream.count.txt"
#'
#' re.rMAT<-ProcessOutputFilesFrom3UTR(dir.name,input.file.pattern)
#'
#'
#'
ProcessOutputFilesFrom3UTR<-function(dir.name,input.file.pattern){

  ProcessOutputFilesFrom_rMATS_read<-function(input_file){

    re=read.table(input_file,header=T)

    return(re)

  }

  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern))
  file.name.2<-as.list(file.name)

  names(file.name.2)=sapply(strsplit(file.name,split="\\/"),"[[",9)

  re.out<-lapply(file.name.2,function(u){
    re=read.table(u,header=F)
    colnames(re)=c("Count","GeneName")
    re
  })

#   #Define DE by FDR
#   if(De_defined_by_what=="FDR"){
#   re.out.2<-do.call(c,lapply(re.out, function(u){
#     x<-as.character(u[which(u$FDR<0.05),]$GeneID)
#     x
#     #x<-as.data.frame(t(x))
#     #colnames(x)<-colnames(Data4Goterm)
#     #x
#   }))
#   }else if(De_defined_by_what=="P_and_inclusion"){
#     re.out.2<-do.call(c,lapply(re.out, function(u){
#       x<-as.character(u[which(u$PValue<0.05&abs(u$IncLevelDifference)>0.20),]$GeneID)
#       x
#       #x<-as.data.frame(t(x))
#       #colnames(x)<-colnames(Data4Goterm)
#       #x
#     }))
#     }
#
# # print(names((re.out.2)))
#
#   n1=length(re.out.2[grep("A3SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
#   n2=length(re.out.2[grep("A5SS.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
#   n3=length(re.out.2[grep("MXE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
#   n4=length(re.out.2[grep("RI.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
#   n5=length(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))])
#
#   cat(n1,"\t",n2,"\t",n3,"\t",n4,"\t",n5,"\n")
#
#   re.out.3<-list(JunctionCountOnly=unique(re.out.2[grep("JunctionCountOnly",names(re.out.2))]),
#                  ReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
#                  SEMATSJunctionCountOnly=unique(re.out.2[grep("SE.MATS.JunctionCountOnly",names(re.out.2))]),
#                  SEReadsOnTargetAndJunctionCounts=unique(re.out.2[grep("SE.MATS.ReadsOnTargetAndJunctionCounts",names(re.out.2))]),
#                  SEs=unique(re.out.2[grep("SE.MATS",names(re.out.2))])
#                  )
#
   return(re.out)

}
