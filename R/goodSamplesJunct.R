#' goodSamplesJunct: Iterative filtering of samples and junctions with too many missing entries
#'
#' Returns a list of samples and junctions that pass criteria on maximum number of missing or low weight values. This is a wrapper function of WGCNA goodSamplesGenes.
#' @param datExpr expression data. A matrix or data frame in which columns are junctions and rows ar samples.
#' @param weights optional observation weights in the same format (and dimensions) as datExpr.
#' @param minFraction minimum fraction of non-missing samples for a gene to be considered good.
#' @param minNSamples minimum number of non-missing samples for a gene to be considered good.
#' @param minNGenes minimum number of good genes for the data set to be considered fit for analysis. If the actual number of good genes falls below this threshold, an error will be issued.
#' @param tol an optional 'small' number to compare the variance against. Defaults to the square of 1e-10 * max(abs(datExpr), na.rm = TRUE). The reason of comparing the variance to this number, rather than zero, is that the fast way of computing variance used by this function sometimes causes small numerical overflow errors which make variance of constant vectors slightly non-zero; comparing the variance to tol rather than zero prevents the retaining of such genes as 'good genes'.
#' @param minRelativeWeight observations whose relative weight is below this threshold will be considered missing. Here relative weight is weight divided by the maximum weight in the column (gene).
#' @param verbose integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no indentation, each unit adds two spaces.
#'
#' @seealso \code{\link[WGCNA]{goodSamplesGenes}}
#' @return A list with 1) A logical vector with one entry per sample that is TRUE if the sample is considered good and FALSE otherwise and B) A logical vector with one entry per gene that is TRUE if the gene is considered good and FALSE otherwise.
#' @examples
#' \dontrun{
#' data(DJEanlz)
#' datExpr0 = DJEanlz$v.norm
#' colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))] <- paste0("TCGA_",
#' seq(1,length(colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))]), 1))
#' Group1 <- colnames(DJEanlz$v.norm$E)[grep("SRR", colnames(DJEanlz$v.norm$E))]
#' datExpr0 = datExpr0[,-c(match(Group1, colnames(datExpr0)))]
#' datExpr0 = as.data.frame(t(datExpr0$E))
#' WGCNA::enableWGCNAThreads(nThreads = nThreads)
#' gsg = goodSamplesJunct(datExpr0, verbose = 3)
#' }
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
goodSamplesJunct <-
  function (datExpr,
            weights = NULL,
            minFraction = 1 / 2,
            minNSamples = 4,
            minNGenes = 4,
            tol = NULL,
            minRelativeWeight = 0.1,
            verbose = 1,
            indent = 0)
  {
    spaces = dynamicTreeCut::indentSpaces(indent)
    goodGenes = NULL
    goodSamples = NULL
    nBadGenes = 0
    nBadSamples = 0
    changed = TRUE
    iter = 1
    if (verbose > 0)
      dynamicTreeCut::printFlush(paste(
        spaces,
        "Flagging junctions and samples with too many missing values..."
      ))
    while (changed) {
      if (verbose > 0)
        dynamicTreeCut::printFlush(paste(spaces, " ..step", iter))
      goodGenes = WGCNA::goodGenes(
        datExpr,
        weights,
        goodSamples,
        goodGenes,
        minFraction = minFraction,
        minNSamples = minNSamples,
        minNGenes = minNGenes,
        minRelativeWeight = minRelativeWeight,
        tol = tol,
        verbose = verbose - 1,
        indent = indent +
          1
      )
      goodSamples = WGCNA::goodSamples(
        datExpr,
        weights,
        goodSamples,
        goodGenes,
        minFraction = minFraction,
        minNSamples = minNSamples,
        minNGenes = minNGenes,
        minRelativeWeight = minRelativeWeight,
        verbose = verbose - 1,
        indent = indent + 1
      )
      changed = ((sum(!goodGenes) > nBadGenes) | (sum(!goodSamples) >
                                                    nBadSamples))
      nBadGenes = sum(!goodGenes)
      nBadSamples = sum(!goodSamples)
      iter = iter + 1
    }
    allOK = (sum(c(nBadGenes, nBadSamples)) == 0)
    list(goodGenes = goodGenes,
         goodSamples = goodSamples,
         allOK = allOK)
  }
