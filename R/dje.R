#' dje: Test for Differential Junction Expression
#'
#' Given a linear model fit at the junction level, test for differences in junction expression between experimental conditions. This is a wrapper function of limma diffSplice.
#' @param fit an MArrayLM fitted model object produced by lmFit or contrasts.fit. Rows should correspond to junctions.
#' @param geneid gene identifiers. Either a vector of length nrow(fit) or the name of the column of fit$genes containing the gene identifiers. Rows with the same ID are assumed to belong to the same gene.
#' @param junctionID junction identifiers. Either a vector of length nrow(fit) or the name of the column of fit$genes containing the junction identifiers.
#' @param robust logical, should the estimation of the empirical Bayes prior parameters be robustified against outlier sample variances?
#' @param verbose logical, if TRUE some diagnostic information about the number of genes and junctions is output.
#' @section Details:
#' An object of class MArrayLM containing both junction level and gene level tests. Results are sorted by geneid and by junctionid within gene.
#' See limma::diffSplice documentation for details.

#' @section Value:
#' This function tests for differential junction usage for each gene and for each column of fit following limma diffSplice implementation.
#' @examples
#' \dontrun{
#' fit.norm0 <- system.file("extdata", "fit.norm.rds", package = "DJExpress")
#' fit.norm <- readRDS(fit.norm0)
#' groupID0 <- system.file("extdata", "groupID.rds", package = "DJExpress")
#' groupID <- readRDS(groupID0)
#' junctionID0 <- system.file("extdata", "junctionID.rds", package = "DJExpress")
#' junctionID <- readRDS(junctionID0)
#' dje.out <- DJExpress::dje(fit = fit.norm, geneid=groupID, junctionID=junctionID)
#' }
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
dje <-
  function (fit,
            geneid,
            junctionID = NULL,
            robust = FALSE,
            verbose = TRUE)
  {
    junction.genes <- fit$genes
    if (is.null(junction.genes))
      junction.genes <- data.frame(junctionID = 1:nrow(fit))
    if (length(geneid) == 1) {
      genecolname <- as.character(geneid)
      geneid <- junction.genes[[genecolname]]
    }
    else {
      junction.genes$GeneID <- geneid
      genecolname <- "GeneID"
    }
    if (is.null(junctionID)) {
      junctioncolname <- NULL
    }
    else {
      if (length(junctionID) == 1) {
        junctioncolname <- as.character(junctionID)
        junctionID <- junction.genes[[junctioncolname]]
      }
      else {
        junction.genes$junctionID <- junctionID
        junctioncolname <- "junctionID"
      }
    }
    if (anyNA(geneid)) {
      isna <- which(is.na(geneid))
      geneid[isna] <- paste0("NA", 1:length(isna))
    }
    if (is.null(junctionID))
      o <- order(geneid)
    else
      o <- order(geneid, junctionID)
    geneid <- geneid[o]
    junction.genes <- junction.genes[o, , drop = FALSE]
    junction.coefficients <- fit$coefficients[o, , drop = FALSE]
    junction.stdev.unscaled <- fit$stdev.unscaled[o, , drop = FALSE]
    junction.df.residual <- fit$df.residual[o]
    junction.s2 <- fit$sigma[o] ^ 2
    junction.stat <- cbind(1, junction.df.residual, junction.s2)
    gene.sum <- rowsum(junction.stat, geneid, reorder = FALSE)
    gene.njunctions <- gene.sum[, 1]
    gene.df.residual <- gene.sum[, 2]
    gene.s2 <- gene.sum[, 3] / gene.sum[, 1]
    if (verbose) {
      cat("Total number of junctions: ", length(geneid), "\n")
      cat("Total number of genes: ", length(gene.njunctions), "\n")
      cat("Number of genes with 1 junction: ",
          sum(gene.njunctions ==
                1),
          "\n")
      cat("Mean number of junctions in a gene: ", round(mean(gene.njunctions),
                                                        0), "\n")
      cat("Max number of junctions in a gene: ",
          max(gene.njunctions),
          "\n")
    }
    squeeze <- limma::squeezeVar(var = gene.s2, df = gene.df.residual,
                                 robust = robust)
    gene.keep <- gene.njunctions > 1
    ngenes <- sum(gene.keep)
    if (ngenes == 0)
      stop("No genes with more than one junction")
    junction.keep <- rep(gene.keep, gene.njunctions)
    geneid <- geneid[junction.keep]
    junction.genes <- junction.genes[junction.keep, , drop = FALSE]
    junction.coefficients <-
      junction.coefficients[junction.keep, , drop = FALSE]
    junction.stdev.unscaled <-
      junction.stdev.unscaled[junction.keep, , drop = FALSE]
    junction.df.residual <- junction.df.residual[junction.keep]
    gene.njunctions <- gene.njunctions[gene.keep]
    gene.df.test <- gene.njunctions - 1
    gene.df.residual <- gene.df.residual[gene.keep]
    if (robust)
      squeeze$df.prior <- squeeze$df.prior[gene.keep]
    gene.df.total <- gene.df.residual + squeeze$df.prior
    gene.df.total <- pmin(gene.df.total, sum(gene.df.residual))
    gene.s2.post <- squeeze$var.post[gene.keep]
    u2 <- 1 / junction.stdev.unscaled ^ 2
    u2.rowsum <- rowsum(u2, geneid, reorder = FALSE)
    gene.betabar <-
      rowsum(junction.coefficients * u2, geneid, reorder = FALSE) / u2.rowsum
    g <- rep(1:ngenes, times = gene.njunctions)
    junction.coefficients <- junction.coefficients - gene.betabar[g,
                                                                  , drop = FALSE]
    junction.t <-
      junction.coefficients / junction.stdev.unscaled / sqrt(gene.s2.post[g])
    gene.F <-
      rowsum(junction.t ^ 2, geneid, reorder = FALSE) / gene.df.test
    junction.1mleverage <- 1 - (u2 / u2.rowsum[g, , drop = FALSE])
    junction.coefficients <- junction.coefficients / junction.1mleverage
    junction.t <- junction.t / sqrt(junction.1mleverage)
    junction.p.value <- 2 * stats::pt(abs(junction.t), df = gene.df.total[g],
                               lower.tail = FALSE)
    gene.F.p.value <-
      stats::pf(gene.F,
         df1 = gene.df.test,
         df2 = gene.df.total,
         lower.tail = FALSE)
    out <- methods::new("MArrayLM", list())
    out$genes <- junction.genes
    out$genecolname <- genecolname
    out$junctioncolname <- junctioncolname
    out$coefficients <- junction.coefficients
    out$t <- junction.t
    out$p.value <- junction.p.value
    out$gene.df.prior <- squeeze$df.prior
    out$gene.df.residual <- gene.df.residual
    out$gene.df.total <- gene.df.total
    out$gene.s2 <- gene.s2[gene.keep]
    out$gene.s2.post <- gene.s2.post
    out$gene.F <- gene.F
    out$gene.F.p.value <- gene.F.p.value
    gene.lastjunction <- cumsum(gene.njunctions)
    gene.firstjunction <- gene.lastjunction - gene.njunctions + 1
    no <- logical(nrow(junction.genes))
    isdup <-
      vapply(junction.genes, duplicated, no)[-gene.firstjunction,
                                             , drop = FALSE]
    isgenelevel <- apply(isdup, 2, all)
    out$gene.genes <- junction.genes[gene.lastjunction, isgenelevel,
                                     drop = FALSE]
    out$gene.genes$Njunctions <- gene.njunctions
    out$gene.firstjunction <- gene.firstjunction
    out$gene.lastjunction <- gene.lastjunction
    penalty <- rep_len(1L, length(g))
    penalty[gene.lastjunction] <- 1L - gene.njunctions
    penalty <- cumsum(penalty)[-gene.lastjunction]
    penalty <- penalty / rep(gene.njunctions - 1L, gene.njunctions - 1L)
    g2 <- g[-gene.lastjunction]
    out$gene.simes.p.value <- gene.F.p.value
    for (j in 1:ncol(fit)) {
      o <- order(g, junction.p.value[, j])
      p.adj <-
        pmin(junction.p.value[o, j][-gene.lastjunction] / penalty,
             1)
      o <- order(g2, p.adj)
      out$gene.simes.p.value[, j] <- p.adj[o][gene.firstjunction -
                                                0L:(ngenes - 1L)]
    }
    out
  }
