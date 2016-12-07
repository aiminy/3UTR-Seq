#' Title
#'
#' @param z
#' @param one.sided
#'
#' @return
#' @export
#'
#' @examples
#' convert.z.score(1.454648e-04,two.sided=T)
#' convert.z.score(1.454648e-04,two.sided=F)
#'
#' convert.z.score(-1.428500e-04,two.sided=T)
#' convert.z.score(-1.428500e-04,two.sided=F)
#' lapply(re.3UTR$re.FC.sorted$pvalue,function(x){
#'  y<-convert.z.score(x,two.side=F)
#'  y
#'  }
#'  )
#'
#' one.sided.p<-convert.z.score(re.3UTR$re.FC.sorted$pvalue,two.sided=F)
#'
#'
convert.z.score<-function(z, two.sided=T) {
  if((two.sided)) {
    pval = 2*pnorm(-abs(z));
  } else if(sign(z)==-1) {
    pval = pnorm(-abs(z));
  } else {
    pval = 1-pnorm(-abs(z));
  }
  return(pval);
}
