#' topDJE: Top table of differentially used junctions
#'
#' Top table ranking the most differentially spliced junctions. This is a wrapper function of limma topSplice.
#' @param fit MArrayLM fit object produced by dje
#' @param coef the coefficient (column) of fit for which differentially splicing is assessed.
#' @param test character string, possible values are "simes", "F" or "t". "F" gives F-tests for each gene. "t" gives t-tests for each exon. "simes" gives genewise p-values derived from the t-tests after Simes adjustment for each gene.
#' @param number integer, maximum number of rows to output.
#' @param FDR numeric, only show junctions or genes with false discovery rate less than this cutoff.
#' @param sort.by character string specifying which column to sort results by. Possible values for "p", "logFC", "Njunctions" or "none". "logFC" is only available if test="t" and "Njunctions" is only available if test="simes" or test="F".
#'
#' @seealso \code{\link[limma]{topSplice}}
#' @return A data.frame with any annotation columns found in fit plus the following columns as in limma::topSplice
#' @examples
#' \dontrun{
#' fit.norm0 <- system.file("extdata", "dje.out.rds", package = "DJExpress")
#' fit.norm <- readRDS(fit.norm0)
#' topdje.out <- topDJE(fit = fit.norm)
#' }
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
topDJE <-
  function (fit,
            coef = ncol(fit),
            test = "simes",
            number = 10L,
            FDR = 1,
            sort.by = "p")
  {
    if (is.null(fit$gene.genes$Njunctions))
      stop("fit should be fit object produced by dje")
    coef <- coef[1]
    test <- match.arg(test, choices = c("simes", "F", "f", "t"))
    if (test == "f")
      test <- "F"
    sort.by <- match.arg(sort.by, choices = c("p", "none", "logFC",
                                              "Njunctions"))
    if (sort.by == "logFC" & test != "t")
      stop("Sorting by logFC only available with Simes test")
    if (sort.by == "Njunctions" & test == "t")
      stop("Sorting by Njunctions only available with gene-level tests")
    switch(test,
           t = {
             out <- fit$genes
             out$logFC <- as.matrix(fit$coefficients)[, coef]
             out$t <- as.matrix(fit$t)[, coef]
             out$P.Value <- as.matrix(fit$p.value)[, coef]
           },
           F = {
             out <- fit$gene.genes
             out$F <- as.matrix(fit$gene.F)[, coef]
             out$P.Value <- as.matrix(fit$gene.F.p.value)[, coef]
           },
           simes = {
             out <- fit$gene.genes
             out$P.Value <- as.matrix(fit$gene.simes.p.value)[, coef]
           })
    out$FDR <- stats::p.adjust(out$P.Value, method = "BH")
    if (FDR < 1)
      out <- out[out$FDR <= FDR,]
    number <- min(number, nrow(out))
    if (number <= 1L)
      return(out)
    o <- switch(
      sort.by,
      p = order(out$P.Value, decreasing = FALSE),
      logFC = order(abs(out$logFC), decreasing = TRUE),
      Njunctions = order(out$Njunctions,-out$P.Value, decreasing = TRUE),
      none = 1:nrow(out)
    )
    o <- o[1:number]
    out[o,]
  }
