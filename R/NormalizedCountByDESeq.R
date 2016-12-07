#' Title
#'
#' @param mat
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' raw.count<-re.rMAT.com.total$DE[,10:17]
#'
#' raw.count.dds<-counts(re.rMAT.com.total.test.2$dds)
#'
#' raw.count.DE<-re.rMAT.com.total.test.2$DE[,c(1,10:17)]
#'
#' raw.count.DE.2<-raw.count.DE[match(rownames(raw.count.dds),raw.count.DE$tx_name),]
#'
#' head(raw.count.DE.2)
#'
#' raw.count.2<-apply(raw.count,2,as.numeric)
#'
#' head(raw.count)
#'
#'
#' nor.fac<-norm_factors(raw.count.2)
#'
#' estimateSizeFactorsForMatrix(raw.count.2)
#'
#'
norm_factors <- function(mat) {
  #nz <- apply(mat, 1, function(row) !any(round(row) == 0))
  #mat_nz <- mat[nz,]
  #p <- ncol(mat)
  #geo_means <- exp(apply(mat_nz, 1, function(row) mean(log(row)) ))

  loggeomeans <- rowMeans(log(mat))

  counts<-mat

  sf<-apply(counts, 2, function(cnts) {
    exp(median(log(cnts)-loggeomeans))})

  # <- sweep(mat_nz, 1, geo_means, `/`)

  #sf <- apply(s, 2, median)
  #scaling <- exp( (-1 / p) * sum(log(sf)))

  #sf * scaling
  sf
}
