#' Title
#'
#' @param re.out
#' @param dir.name
#' @param by_what
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' dir.name=/media/H_driver/2016/Ramin_azhang/for_bioinfo_core/RNA_seq/Results4NewData/"
#'
#' by_what="*.rm.exon.intron.hg19.bed"
#'
#' NormalizedCount(re.out,dir.name,by_what)

NormalizedCount <- function(re.out,dir.name,by_what) {


  pattern_check=grep(glob2rx("*.rm.exon.intron.hg19.bed"),by_what)

  if(length(pattern_check)!=0){
    cmd2=paste0("wc -l ",dir.name,by_what)
  }else{
    cmd2=paste0("wc -l $(ls ",dir.name,by_what," | grep -v hg19)")
  }

  num.intergenic.reads<-system(cmd2, intern = TRUE, ignore.stderr = TRUE)

  num.intergenic.reads.2<-strsplit(num.intergenic.reads,split="\\/")

  #cat("OK")
  #print(num.intergenic.reads.2)
  #cat("OK")

  n=length(num.intergenic.reads.2)-1
  num.intergenic.reads.3<-cbind(trimws(do.call("rbind",lapply(num.intergenic.reads.2[1:n],"[[",1))),do.call("rbind",lapply(num.intergenic.reads.2[1:n],"[[",9)))

  name.sample<-unlist(lapply(strsplit(num.intergenic.reads.3[,2],split="\\."),"[[",1))

  print(names(re.out))

  print(name.sample)

  normalized.count.gene.matched.final<-list()

  for(i in 1:8){
    count.gene.matched<-re.out[(grep(name.sample[i],names(re.out)))]
    num.intergenic.reads.matched<-as.numeric(num.intergenic.reads.3[grep(name.sample[i],num.intergenic.reads.3[,2]),1])

    #print(count.gene.matched)
    #print(num.intergenic.reads.matched)

    normalized.count.gene.matched<-(count.gene.matched[[1]][,1]/num.intergenic.reads.matched)*(10^6)
    count.gene.matched.2<-data.frame(cbind(as.character(count.gene.matched[[1]][,2]),count.gene.matched[[1]][,1],normalized.count.gene.matched))

    colnames(count.gene.matched.2)=c("Gene",paste0("Counts.",name.sample[i]),paste0("Normalized.Counts.",name.sample[i]))

    normalized.count.gene.matched.final[[i]]<-count.gene.matched.2
  }

  merge.all <- function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
  }

  output <- Reduce(merge.all,normalized.count.gene.matched.final)

  row.has.na <- apply(output, 1, function(x){any(is.na(x))})
  print(sum(row.has.na))

  final.filtered <- output[!row.has.na,]

  re<-list(final.filtered=final.filtered,num.intergenic.reads.2=num.intergenic.reads.2)

  return(re)

}
