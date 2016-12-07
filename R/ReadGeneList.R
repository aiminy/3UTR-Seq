#' Title
#'
#' @param dir.name
#' @param input.file.pattern.4.interested.gene
#'
#' @return
#' @export
#'
#' @examples
#'  dir.name="/media/H_driver/2016/Ramin_azhang/"
#'  input.file.pattern.4.interested.gene="final_list.csv"
#' gene.interested<-ReadGeneList(dir.name, input.file.pattern.4.interested.gene)
#'
ReadGeneList<- function(dir.name, input.file.pattern.4.interested.gene) {
  file.name=paste0(dir.name,dir(dir.name,recursive = TRUE,pattern=input.file.pattern.4.interested.gene))
  genes.interested<-read.csv(file.name,header=F)
  return(genes.interested)
}
