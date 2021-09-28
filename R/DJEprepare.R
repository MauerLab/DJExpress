#' DJEprepare: Filtering junctions for differential expression analysis
#'
#' Filter junctions based on their expression. Returns output object for DJEanalize
#' @param annotate.out output object from DJEannotate()
#' @param Group1 vector or factor specifying sample names in the control group
#' @param minMean numeric: minimum of read count mean per junction
#' @param maxMean numeric: maximum of read count mean per junction
#' @param minVar numeric: minimum of read count variance per junction
#' @param maxVar numeric: maximum of read count variance per junction
#'
#' @return List object containing expression-based filtered junction counts, gene annotation and design matrix for differential expression analysis
#' @examples
#' DJEann <- system.file("extdata", "DJEann.rds", package = "DJExpress")
#' ann.out <- readRDS(DJEann)
#' Group1 <- colnames(ann.out$quant.annotated)[grep("GTEx", colnames(ann.out$quant.annotated))]
#' prep.out <- DJEprepare(ann.out, Group1)
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEprepare <-
  function(annotate.out,
           Group1,
           minMean = 10,
           maxMean = Inf,
           minVar = 0,
           maxVar = Inf) {
    junctExpr <- annotate.out$quant.annotated
    featureID <- annotate.out$featureID
    groupID <- annotate.out$groupID
    junctExprMean <- apply(junctExpr, 1, mean)
    junctExprVar <- apply(junctExpr, 1, stats::var)
    varMeanFilter <- junctExprMean >= minMean & junctExprMean <=
      maxMean & junctExprVar >= minVar & junctExprVar <= maxVar
    isFromGroup1 <- colnames(junctExpr) %in% Group1
    design       <- cbind(1, ifelse(isFromGroup1, 0, 1))
    return(
      list(
        JunctExprfilt = junctExpr[varMeanFilter, ],
        featureID = featureID[varMeanFilter],
        groupID = groupID[varMeanFilter],
        design = design
      )
    )
  }
