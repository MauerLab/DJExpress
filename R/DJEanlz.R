#' output object of DJEanalyze function
#'
#' A list of objects containing the results of differential junction expression analysis. TCGA and GTEx example junction expression is used.
#'
#' @docType data
#'
#' @usage data(DJEanlz)
#'
#' @format An object of class \code{"list"}.
#'
#' @keywords datasets
#'
#' @references Kahles et al.Cancer Cell (2018) Aug 13;34(2):211-224.e6.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/30078747/}{PubMed})
#'
#' @examples
#' data(DJEanlz)
#' iPlot.out <- DJEplotSplice(DJEanlz, geneID="RPSA")
"DJEanlz"
