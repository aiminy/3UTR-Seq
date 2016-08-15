

cmd2="wc -l /media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results/*.rm.exon.intron.hg19.bed"

Re.out<-system(cmd2, intern = TRUE, ignore.stderr = TRUE)

re<-strsplit(Re.out,split="\\/")

re2<-cbind(trimws(do.call("rbind",lapply(re[1:4],"[[",1))),do.call("rbind",lapply(re[1:4],"[[",9)))
