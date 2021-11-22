#' DJEanalyze: Test for Differential Junction Usage
#'
#' Normalize junction expression and performs differential junction usage analysis using limma::diffSplice methods. Returns output object for DJEplotSplice(), DJEvsTrait() and JCNAprepare() functions.
#' @param prepare.out output object from DJEprepare()
#' @param Group1 vector or factor specifying basic control sample names
#' @param junct.annot gene-annotated junction IDs
#' @param normalize.method limma voom's normalization method to be applied to the logCPM values (default "none")
#' @param calcNormFactors logical, should edgeR's library Size Normalization be calculated?
#' @param calMethod Normalization method to be used. Options are: "TMM","RLE","upperquartile" or "none" (default "none")
#' @param plot logical, should median/variance trend plot be displayed?
#' @param FDR numeric: adjusted p-value cutoff in the linear model fit to define significance
#' @param logFC numeric: log2-fold-change value cutoff corresponding to the effect or contrast in the linear model fit to define significance
#' @param level numeric: confidence interval for Linear Regression of absolute vs gene-wise logFC in junction expression
#' @param topJunct numeric: number of top significant junctions that should be labeled in volcano plot
#' @param targetGene vector or factor specifying the genes whose significant junctions should be labeled in volcano plot
#' @param legend.position character indicating the position of labels in volcano plot. Allowed values are: “left”, “top”, “right”, “bottom”. Default "right".
#' @param volc.axes.face character indicating the font face of axes tick labels in volcano plot. Allowed values are: "plain", "italic", "bold", "bold.italic"
#' @param volc.axes.col character indicating the color of axes tick labels in volcano plot
#' @param volc.axes.size numeric indicating the size of axes tick labels in volcano plot
#' @param volc.axes.angle numeric indicating the angle of axes tick labels in volcano plot
#' @param scale_color_manual vector of aesthetic values (colors) to indicate significance of data values in volcano plot (default are #005BA2 for downregulated, #D44F4F for upregulated and #B3B4B5 for non-significant junctions)
#' @section Details:
#' This wrapper function receives DJEprepare() output object and performs junction expression normalization and differential junction usage analysis. The differential analysis is an implementation of limma::diffSplice() method.
#' See limma::diffSplice documentation for details.
#' @return A list of objects containing:
#' v.norm: log-cpm values output by limma:voom
#' ex.norm: Differential Junction Expression analysis output
#' dje.out: Annotated ex.norm data set with additional information, including basic statistics (e.g. median, zero counts, etc) and DJE group for each junction.
#' dje.sig: Significant hits in dje.out based on FDR and logFC cutoffs.
#' volcano.plot: Volcano plot of differential junction expression
#' logFC.plot.junctions: Table with junctions shown in logFC.plot. They are defined as junctions passing FDR cutoff for differential usage as well as differential expression.
#' logFC.plot: Regression plot of Absolute logFC ~ Relative logFC for differentially used (compared to average junction expression in the gene) and differentially expressed junctions (basal vs tested sample group)
#  model.fit: Confidence and prediction intervals for Linear Regression (Absolute logFC ~ Relative logFC)
#  group.par: Labeling of junctions in regression plot based on fitted model, FDR and logFC cutoffs
#' @importFrom ggplot2 element_text
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @importFrom methods new
#' @seealso \code{\link{limma}}
#' @examples
#' DJEprep <- system.file("extdata", "DJEprep.rds", package = "DJExpress")
#' prep.out <- readRDS(DJEprep)
#' Group1 <- colnames(prep.out$JunctExprfilt)[grep("GTEx", colnames(prep.out$JunctExprfilt))]
#' anlz.out <- DJExpress::DJEanalyze(prep.out, Group1)
#' @export
DJEanalyze <-
  function(prepare.out,
           Group1,
           junct.annot = NULL,
           normalize.method = "none",
           calcNormFactors = FALSE,
           calMethod = "none",
           plot = FALSE,
           FDR = 0.05,
           logFC = 2,
           level = 0.9,
           topJunct = NULL,
           targetGene = NULL,
           legend.position = "right",
           volc.axes.face = "plain",
           volc.axes.col = "black",
           volc.axes.size = 11,
           volc.axes.angle = 0,
           scale_color_manual = c("#005BA2", "#D44F4F", "#B3B4B5")) {
    logFC.ebayes <- Significant <- GeneID <- NULL
    junctExpr <- prepare.out$JunctExprfilt
    featureID <- prepare.out$featureID
    groupID <- prepare.out$groupID
    design <- prepare.out$design
    dge <- edgeR::DGEList(junctExpr)
    if (calcNormFactors) {
      dge = edgeR::calcNormFactors(dge, method = calMethod)
    }
    if (plot) {
      v.norm <-
        limma::voom(dge, design, plot = TRUE, normalize.method = normalize.method)
    } else{
      v.norm <-
        limma::voom(dge, design, plot = FALSE, normalize.method = normalize.method)
    }
    fit.norm <- limma::lmFit(v.norm, design)
    ex.norm <-
      DJExpress::dje(fit.norm, geneid = groupID, junctionID = rownames(v.norm))
    topSplice.out <-
      DJExpress::topDJE(
        ex.norm,
        coef = 2,
        test = "t",
        number = nrow(fit.norm)
      )
    medianExp.group1 <-
      apply(junctExpr[, match(Group1, colnames(junctExpr))],
            1, stats::median, na.rm = TRUE)
    meanExp.group1 <-
      apply(junctExpr[, match(Group1, colnames(junctExpr))],
            1, mean, na.rm = TRUE)
    medianExp.group2 <-
      apply(junctExpr[, -c(match(Group1, colnames(junctExpr)))],
            1, stats::median, na.rm = TRUE)
    meanExp.group2 <-
      apply(junctExpr[, -c(match(Group1, colnames(junctExpr)))],
            1, mean, na.rm = TRUE)
    medianNormExp.group1 <-
      apply(v.norm$E[, match(Group1, colnames(junctExpr))],
            1, stats::median, na.rm = TRUE)
    meanNormExp.group1 <-
      apply(v.norm$E[, match(Group1, colnames(junctExpr))],
            1, mean, na.rm = TRUE)
    medianNormExp.group2 <-
      apply(v.norm$E[, -c(match(Group1, colnames(junctExpr)))],
            1, stats::median, na.rm = TRUE)
    meanNormExp.group2 <-
      apply(v.norm$E[, -c(match(Group1, colnames(junctExpr)))],
            1, mean, na.rm = TRUE)
    group1WithCounts <-
      rowSums(junctExpr[, match(Group1, colnames(junctExpr))] != 0)
    group2WithCounts <-
      rowSums(junctExpr[, -c(match(Group1, colnames(junctExpr)))] != 0)
    zeroCounts.group1 <-
      rowSums(junctExpr[, match(Group1, colnames(junctExpr))] == 0)
    zeroCounts.group2 <-
      rowSums(junctExpr[, -c(match(Group1, colnames(junctExpr)))] == 0)
    junct.addinfo <- as.data.frame(
      cbind(
        rownames(junctExpr),
        medianExp.group1,
        meanExp.group1,
        medianExp.group2,
        meanExp.group2,
        medianNormExp.group1,
        meanNormExp.group1,
        medianNormExp.group2,
        meanNormExp.group2,
        group1WithCounts,
        group2WithCounts,
        zeroCounts.group1,
        zeroCounts.group2
      )
    )
    junct.addinfo$neojunction <- NA
    junct.addinfo$neojunction[which(group1WithCounts == 0 &
                                      group2WithCounts != 0)] <- "neojunction"
    if (!is.null(junct.annot)) {
      can.j <-
        data.table::chmatch(rownames(junct.addinfo), as.character(junct.annot$junction), nomatch =
                              NA_integer_)
      junct.addinfo <-
        as.data.frame(cbind(junct.addinfo, annotation = junct.annot[can.j, 2]))
      unn <- which(is.na(junct.addinfo$annotation))
      junct.addinfo$annotation[unn] <- "unannotated"
      junct.addinfo$annotation[-c(unn)] <- NA
    }
    matched.ids <-
      data.table::chmatch(topSplice.out$junctionID, as.character(junct.addinfo$V1), nomatch =
                            NA_integer_)
    topSplice.out <-
      as.data.frame(cbind(topSplice.out, junct.addinfo[matched.ids, -c(1)]))
    rownames(topSplice.out) = NULL
    ebfit = limma::eBayes(fit.norm)
    Top.ebayes = limma::topTable(ebfit, coef = 2, number = nrow(fit.norm))
    Top.ebayes <-
      as.data.frame(cbind(rownames(Top.ebayes), Top.ebayes))
    colnames(Top.ebayes) <-
      c(
        "confirmed_ID",
        "logFC.ebayes",
        "AveExpr.ebayes",
        "t.ebayes",
        "P.Value.ebayes",
        "adj.P.Val.ebayes",
        "B.ebayes"
      )
    matched.ids.2 <-
      data.table::chmatch(topSplice.out$junctionID,
                          as.character(Top.ebayes$confirmed_ID),
                          nomatch = NA_integer_)
    topSplice.out <-
      as.data.frame(cbind(topSplice.out, Top.ebayes[matched.ids.2, -c(1)]))
    rownames(topSplice.out) = NULL
    SigTopTable = topSplice.out[which(as.numeric(as.character(topSplice.out$FDR))
                                      <= FDR), ]

    SigTopTable = SigTopTable[which(as.numeric(as.character(SigTopTable$adj.P.Val.ebayes))
                                    <= FDR), ]

    SigTopTable2 = topSplice.out[which(as.numeric(as.character(topSplice.out$FDR))
                                      <= FDR & abs(as.numeric(as.character(topSplice.out$logFC))) >= logFC), ]


    if(nrow(SigTopTable)<3){
      warning("Less than 3 differentially expressed junctions were found")
      return(
        list(
          v.norm = v.norm,
          ex.norm = ex.norm,
          dje.out = topSplice.out,
          dje.sig = NULL,
          logFC.plot.junctions = NULL,
          logFC.plot = NULL,
          volcano.plot = NULL,
          model.fit = NULL,
          group.par = NULL
        )
      )
    }else{
      model1 <- stats::lm(logFC.ebayes ~ logFC, data = SigTopTable)
      xRange = data.frame(logFC = seq(
        min(SigTopTable$logFC, na.rm = T),
        max(SigTopTable$logFC, na.rm = T),
        0.0001
      ))
      pred_interval <-
        stats::predict(
          model1,
          newdata = data.frame(logFC = xRange),
          interval = "prediction",
          level = level
        )
      modelConfInt <- stats::predict(model1,
                                     SigTopTable,
                                     level = level,
                                     interval = "prediction")
      insideInterval <-
        modelConfInt[, 'lwr'] <= SigTopTable[['logFC.ebayes']] &
        SigTopTable[['logFC.ebayes']] <= modelConfInt[, 'upr']
      #Group 0
      insideInterval.NOlogFC <-
        insideInterval & abs(SigTopTable$logFC) <= logFC & abs(SigTopTable$logFC.ebayes) <= logFC
      #Group 1
      insideInterval.logFC.neg <-
        insideInterval & SigTopTable$logFC <= -logFC &
        SigTopTable$logFC.ebayes <= -logFC
      insideInterval.logFC.pos <-
        insideInterval & SigTopTable$logFC >= logFC &
        SigTopTable$logFC.ebayes >= logFC
      #Group2
      insideInterval.NOeb.neg <-
        insideInterval & SigTopTable$logFC <= -logFC &
        abs(SigTopTable$logFC.ebayes) <= logFC
      insideInterval.NOeb.pos <-
        insideInterval & SigTopTable$logFC >= logFC &
        abs(SigTopTable$logFC.ebayes) <= logFC
      insideInterval.NOdif.neg <-
        insideInterval & SigTopTable$logFC.ebayes <= -logFC &
        abs(SigTopTable$logFC) <= logFC
      insideInterval.NOdif.pos <-
        insideInterval & SigTopTable$logFC.ebayes >= logFC &
        abs(SigTopTable$logFC) <= logFC
      outsideInterval.NOdif <-
        !insideInterval & abs(SigTopTable$logFC) <= logFC
      #Group3
      insideInterval.logFC.neg.pos <-
        insideInterval & SigTopTable$logFC <= -logFC &
        SigTopTable$logFC.ebayes >= logFC
      insideInterval.logFC.pos.neg <-
        insideInterval & SigTopTable$logFC >= logFC &
        SigTopTable$logFC.ebayes < -logFC
      outsideInterval.logFC.neg <-
        !insideInterval & SigTopTable$logFC < -logFC &
        SigTopTable$logFC.ebayes <= -logFC
      outsideInterval.logFC.pos <-
        !insideInterval & SigTopTable$logFC >= logFC &
        SigTopTable$logFC.ebayes >= logFC
      outsideInterval.logFC.neg.pos <-
        !insideInterval & SigTopTable$logFC <= -logFC &
        SigTopTable$logFC.ebayes >= logFC
      outsideInterval.logFC.pos.neg <-
        !insideInterval & SigTopTable$logFC >= logFC &
        SigTopTable$logFC.ebayes <= -logFC

      group.par <- list(
        Group0 = insideInterval.NOlogFC,
        Group1.neg = insideInterval.logFC.neg,
        Group1.pos = insideInterval.logFC.pos,
        Group2.RelativeNeg = insideInterval.NOeb.neg,
        Group2.RelativePos = insideInterval.NOeb.pos,
        Group2.AbsoluteNeg = insideInterval.NOdif.neg,
        Group2.AbsolutePos = insideInterval.NOdif.pos,
        Group2.NoRelative = outsideInterval.NOdif,
        Group3.RelativeNeg.AbsolutePos = insideInterval.logFC.neg.pos,
        Group3.RelativePos.AbsoluteNeg = insideInterval.logFC.pos.neg,
        Group3.RelativeNeg.AbsolutePos.2 = outsideInterval.logFC.neg.pos,
        Group3.RelativePos.AbsoluteNeg.2 = outsideInterval.logFC.pos.neg,
        Group3.RelativeNeg.AbsoluteNeg = outsideInterval.logFC.neg,
        Group3.RelativePos.AbsolutePos = outsideInterval.logFC.pos
      )
      plot(logFC.ebayes~logFC, data=SigTopTable[ insideInterval,],
           pch = 19, col = "white",
           ylab = "Absolute logFC (Group 2 vs Group 1)",
           xlab = "Relative logFC (this junction vs others)",
           ylim=c(min(SigTopTable$logFC.ebayes, na.rm = T)-0.5,max(SigTopTable$logFC.ebayes, na.rm = T)+0.5))
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'grey80',
        data = SigTopTable[ insideInterval.NOlogFC,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = '#D44F4F',
        data = SigTopTable[ insideInterval.logFC.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = '#005BA2',
        data = SigTopTable[ insideInterval.logFC.neg,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'lightblue2',
        data = SigTopTable[ insideInterval.NOeb.neg,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'rosybrown2',
        data = SigTopTable[ insideInterval.NOeb.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'grey50',
        data = SigTopTable[ insideInterval.NOdif.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'grey50',
        data = SigTopTable[ insideInterval.NOdif.neg,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'grey50',
        data = SigTopTable[ outsideInterval.NOdif,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = '#D44F4F',
        data = SigTopTable[ insideInterval.logFC.pos.neg,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = '#005BA2',
        data = SigTopTable[ insideInterval.logFC.neg.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'darkgreen',
        data = SigTopTable[ outsideInterval.logFC.neg,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'darkgreen',
        data = SigTopTable[ outsideInterval.logFC.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'darkgreen',
        data = SigTopTable[ outsideInterval.logFC.neg.pos,]
      )
      points(
        logFC.ebayes~logFC,
        pch = 19,
        col = 'darkgreen',
        data = SigTopTable[ outsideInterval.logFC.pos.neg,]
      )
      abline(h = 0, col = "grey80")
      abline(v = 0, col = "grey80")
      abline(h = 2, col = "grey80")
      abline(h = -2, col = "grey80")
      abline(v = 2, col = "grey80")
      abline(v = -2, col = "grey80")
      matplot(
        xRange,
        pred_interval,
        lty=c(1,2,2),
        lwd=c(1,1,1),#vector of line types and widths
        type="l",       #type of plot for each column of y
        add = TRUE,
        col = c("black", "black", "black")
      )
      p <- grDevices::recordPlot()
      grDevices::dev.off()
      model.fit <- as.data.frame(cbind(xRange, pred_interval))
      SigTopTable$group <- NA
      SigTopTable$group[which(insideInterval.NOlogFC == "TRUE")] <-
        "Group0 - No significant DTU and DJE"
      SigTopTable$group[which(insideInterval.logFC.pos == "TRUE")] <-
        "Group 1 (positive logFC)"
      SigTopTable$group[which(insideInterval.logFC.neg == "TRUE")] <-
        "Group 1 (negative logFC)"
      SigTopTable$group[which(insideInterval.NOeb.neg == "TRUE")] <-
        "Group 2 (negative DTU logFC / No significant DJE)"
      SigTopTable$group[which(insideInterval.NOeb.pos == "TRUE")] <-
        "Group 2 (positive DTU logFC / No significant DJE)"
      SigTopTable$group[which(insideInterval.NOdif.neg == "TRUE")] <-
        "Group 2 (No significant DTU / negative DJE logFC)"
      SigTopTable$group[which(insideInterval.NOdif.pos == "TRUE")] <-
        "Group 2 (No significant DTU / positive DJE logFC)"
      SigTopTable$group[which(outsideInterval.NOdif == "TRUE")] <-
        "Group 2 (No significant DTU / significant DJE logFC outside conf.int.)"
      SigTopTable$group[which(insideInterval.logFC.neg.pos == "TRUE")] <-
        "Group 3 (negative DTU logFC / positive DJE logFC)"
      SigTopTable$group[which(insideInterval.logFC.pos.neg == "TRUE")] <-
        "Group 3 (positive DTU logFC / negative DJE logFC)"
      SigTopTable$group[which(outsideInterval.logFC.neg.pos == "TRUE")] <-
        "Group 3 (negative DTU logFC / positive DJE logFC outside conf.int.)"
      SigTopTable$group[which(outsideInterval.logFC.pos.neg == "TRUE")] <-
        "Group 3 (positive DTU logFC / negative DJE logFC outside conf.int.)"
      SigTopTable$group[which(outsideInterval.logFC.neg == "TRUE")] <-
        "Group 3 (negative DTU logFC / negative DJE logFC outside conf.int.)"
      SigTopTable$group[which(outsideInterval.logFC.pos == "TRUE")] <-
        "Group 3 (positive DTU logFC / positive DJE logFC outside conf.int.)"
      SigTopTable <- SigTopTable[order(abs(SigTopTable$logFC), decreasing = TRUE),]
      topSplice.out$FDR[which(topSplice.out$FDR == 0)] = 1.0e-322
      topSplice.out$Significant <-
        ifelse(
          topSplice.out$FDR <= FDR &
            topSplice.out$logFC >= logFC,
          paste0("FDR","<", FDR, "& logFC >", logFC),
          "Not Sig"
        )
      topSplice.out$Significant[which(topSplice.out$FDR <= FDR &
                                        topSplice.out$logFC <= -logFC)] = paste0("FDR","<", FDR, "& logFC <", -logFC)
      if (!is.null(topJunct) | !is.null(targetGene)){
        if (!is.null(topJunct)) {
          topJunct0 <-
            c(
              which(topSplice.out$FDR <= FDR &
                      topSplice.out$logFC > logFC)[c(1:topJunct)],
              which(topSplice.out$FDR <= FDR &
                      topSplice.out$logFC < -logFC)[c(1:topJunct)]
            )
          volcano.p <- ggplot2::ggplot(topSplice.out, ggplot2::aes(x = logFC, y = -log10(FDR))) +
            ggplot2::geom_point(ggplot2::aes(color = Significant)) +
            ggplot2::ylim(0, max(-log10(topSplice.out$FDR), na.rm = T) + max(-log10(topSplice.out$FDR), na.rm = T) *
                            30 / 100) +
            ggplot2::scale_color_manual(values = scale_color_manual) +
            ggrepel::geom_text_repel(
              data = topSplice.out[topJunct0,],
              ggplot2::aes(label = GeneID),
              box.padding = 0.4,
              size = 5
            ) +
            ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = legend.position,
                                                               axis.text.x = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                          size=volc.axes.size, angle=volc.axes.angle),
                                                               axis.text.y = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                          size=volc.axes.size, angle=volc.axes.angle))
          return(
            list(
              v.norm = v.norm,
              ex.norm = ex.norm,
              dje.out = topSplice.out,
              dje.sig = SigTopTable2,
              logFC.plot.junctions = SigTopTable,
              logFC.plot = p,
              volcano.plot = volcano.p,
              model.fit = model.fit,
              group.par = group.par
            )
          )
        }
        if (!is.null(targetGene)) {
          targetGene0 <-
            which(topSplice.out$FDR <= FDR & topSplice.out$GeneID == targetGene)
          if (length(targetGene0) == 0)
            warning("target Gene has no significant DE junctions")
          if (length(targetGene0) != 0) {
            topSplice.out$Significant[targetGene0] = "targetGene"
            scale_color_manual = c(scale_color_manual, "#030303")
            volcano.p <- ggplot2::ggplot(topSplice.out, ggplot2::aes(x = logFC, y = -log10(FDR))) +
              ggplot2::geom_point(ggplot2::aes(color = Significant)) +
              ggplot2::ylim(0, max(-log10(topSplice.out$FDR), na.rm = T) + max(-log10(topSplice.out$FDR), na.rm = T) *
                              30 / 100) +
              ggplot2::scale_color_manual(values = scale_color_manual) +
              ggrepel::geom_text_repel(
                data = topSplice.out[targetGene0,],
                ggplot2::aes(label = GeneID),
                box.padding = 0.4,
                size = 5
              ) +
              ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = legend.position,
                                                                 axis.text.x = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                                                     size=volc.axes.size, angle=volc.axes.angle),
                                                                 axis.text.y = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                            size=volc.axes.size, angle=volc.axes.angle))
            return(
              list(
                v.norm = v.norm,
                ex.norm = ex.norm,
                dje.out = topSplice.out,
                dje.sig = SigTopTable2,
                logFC.plot.junctions = SigTopTable,
                logFC.plot = p,
                volcano.plot = volcano.p,
                model.fit = model.fit,
                group.par = group.par
              )
            )
          }
        }
      }else{
        volcano.p <- ggplot2::ggplot(topSplice.out, ggplot2::aes(x = logFC, y = -log10(FDR))) +
          ggplot2::geom_point(ggplot2::aes(color = Significant)) +
          ggplot2::ylim(0, max(-log10(topSplice.out$FDR), na.rm = T) + max(-log10(topSplice.out$FDR), na.rm = T) *
                          30 / 100) +
          ggplot2::scale_color_manual(values = scale_color_manual) +
          ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = legend.position,
                                                             axis.text.x = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                        size=volc.axes.size, angle=volc.axes.angle),
                                                             axis.text.y = element_text(face=volc.axes.face, color=volc.axes.col,
                                                                                        size=volc.axes.size, angle=volc.axes.angle))
        return(
          list(
            v.norm = v.norm,
            ex.norm = ex.norm,
            dje.out = topSplice.out,
            dje.sig = SigTopTable2,
            logFC.plot.junctions = SigTopTable,
            logFC.plot = p,
            volcano.plot = volcano.p,
            model.fit = model.fit,
            group.par = group.par
          )
        )
      }
    }
  }
