CompareToFeature.R <- function(Re.unadjusted.adjusted) {

  length(unique(as.character(re.DoGs.adjusted.by.batch$gene)))

  grep("-",unique(as.character(re.DoGs.adjusted.by.batch$gene)))

  emp1.WT<-read.table("/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4Feature/2016-02-10-emp1_WT.gene.downstream.count.rv.m.FC.txt",header=T)

  emp1.WT.feature<-emp1.WT[,c(1,7)]

  colnames(emp1.WT.feature)=c("Row.names","feature.Count")

  emp1.WT.inhouse<-Re.unadjusted.adjusted$DE[,c(1,13)]

  Re.unadjusted.adjusted$DE[which(Re.unadjusted.adjusted$DE[,1] %in% c("uc001aai.1")),]

  emp1.inhouse.feature<-merge(emp1.WT.feature,emp1.WT.inhouse)

  colnames(emp1.inhouse.feature)<-c("Transcript","Counting_based_featureCounts","Counting_based_bbc_count")

  plot(emp1.inhouse.feature[,3],emp1.inhouse.feature[,2],ylab="featureCounts_based",xlab="bbc_inhouse_developed_program_based")

  head(emp1.inhouse.feature)
}
