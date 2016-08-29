sanitizeColData <- function(object) {
  if (is.null(mcols(colData(object)))) {
    mcols(colData(object)) <- DataFrame(type=rep("input",ncol(colData(object))),
                                        description=character(ncol(colData(object))))
  }
  class(mcols(colData(object))$type) <- "character"
  class(mcols(colData(object))$description) <- "character"
  mcols(colData(object))$type[ is.na(mcols(colData(object))$type) ] <- ""
  mcols(colData(object))$description[ is.na(mcols(colData(object))$description) ] <- ""
  object
}

estimateSizeFactors.DESeqDataSet <- function(object, type=c("ratio","iterate"),
                                             locfunc=stats::median, geoMeans, controlGenes, normMatrix) {
  type <- match.arg(type, c("ratio","iterate"))
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  object <- sanitizeColData(object)
  if (type == "iterate") {
    sizeFactors(object) <- estimateSizeFactorsIterate(object)
  } else {
    if ("avgTxLength" %in% assayNames(object)) {
      nm <- assays(object)[["avgTxLength"]]
      nm <- nm / exp(rowMeans(log(nm))) # divide out the geometric mean
      normalizationFactors(object) <- estimateNormFactors(counts(object),
                                                          normMatrix=nm,
                                                          locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
      message("using 'avgTxLength' from assays(dds), correcting for library size")
    } else if (missing(normMatrix)) {
      sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object), locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
    } else {
      normalizationFactors(object) <- estimateNormFactors(counts(object),
                                                          normMatrix=normMatrix,
                                                          locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
      message("using 'normMatrix', adding normalization factors which correct for library size")
    }
  }
  object
}
