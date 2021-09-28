#' JCNAModTrait: Module membership vs. Junction significance plot
#'
#' Plots module membership vs junction significance to identify junctions with high significance for selected trait as well as high module membership in target modules
#' @param pass.out output object from JCNA1pass() or JCNA2pass()
#' @param nThreads numeric: number of threads to allow
#' @param trait selected trait for junction module association
#' @param module selected junction module for trait association
#' @param cor.method  a character string specifying the method to be used for correlation as in WGCNA. Options are "bicor", "pearson", "kendall" or "spearman".
#'
#' @return Interactive Highcharts plot with module membership vs junction significance for trait
#' @examples
#' J1p <- system.file("extdata", "J1pass.rds", package = "DJExpress")
#' J1pass <- readRDS(J1p)
#' jMT <- JCNAModTrait(J1pass, trait ="DDX1",module = "brown", cor.method = "bicor")
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom WGCNA bicor
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
JCNAModTrait <-
  function (pass.out,
            nThreads = 2,
            trait,
            module,
            cor.method = c("bicor", "pearson", "kendall", "spearman")) {
    net <- pass.out$net
    allTraits <- pass.out$datTraits
    datExpr <- pass.out$datExpr
    Junctrait <- pass.out$Junctrait
    moduleColors <- pass.out$moduleColors
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    print("Calculating module membership values and trait correlations")
    allTraits <- pass.out$datTraits
    MEs0 = WGCNA::moduleEigengenes(datExpr, pass.out$moduleColors)$eigengenes
    MEs = WGCNA::orderMEs(MEs0)
    trait0 = as.data.frame(allTraits[, match(trait, colnames(allTraits))])
    names(trait0) = trait
    modNames = substring(names(MEs), 3)
    if (cor.method == "bicor"){
      JunctMM = as.data.frame(WGCNA::bicor(datExpr, MEs, use = "p"))
    }else{
      JunctMM = as.data.frame(WGCNA::cor(datExpr, MEs, method = cor.method,  use = "p"))
    }
    MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(JunctMM), nrow(datExpr)))
    names(JunctMM) = paste("MM", modNames, sep = "")
    names(MMPvalue) = paste("p.MM", modNames, sep = "")
    if (cor.method == "bicor"){
      JunctTS = as.data.frame(WGCNA::bicor(datExpr, trait0, use = "p"))
    }else{
      JunctTS = as.data.frame(WGCNA::cor(datExpr, trait0, method = cor.method,  use = "p"))
    }
    JSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(JunctTS), nrow(datExpr)))
    names(JunctTS) = paste("JS.", names(trait0), sep = "")
    names(JSPvalue) = paste("p.JS.", names(trait0), sep = "")
    module = module
    column = match(module, modNames)
    moduleGenes = moduleColors == module
    highplot <-
      as.data.frame(cbind(
        JunctMM = JunctMM[moduleGenes, column],
        JunctTS = JunctTS[moduleGenes, 1],
        junctionID = colnames(datExpr)[moduleGenes]
      ))
    highplot$JunctMM <- abs(as.numeric(as.character(highplot$JunctMM)))
    highplot$JunctTS <- abs(as.numeric(as.character(highplot$JunctTS)))
    cor.coef <- stats::cor.test(highplot$JunctMM, highplot$JunctTS)
    MMplot <- highcharter::highchart() %>%
      highcharter::hc_add_series(
        highplot,
        type = "scatter",
        showInLegend = F,
        highcharter::hcaes(x = JunctMM, y = JunctTS),
        color = module
      ) %>%
      highcharter::hc_title(
        text = paste("Module membership vs. Junction significance for trait\n"),
        align = "center",
        style = list(color = "black", useHTML = TRUE)
      ) %>%
      highcharter::hc_subtitle(
        text = paste(
          "rho:",
          round(cor.coef$estimate, digits = 2),
          "(p-value:",
          round(cor.coef$p.value, digits = 4),
          ")"
        ),
        align = "center",
        style = list(color = "black", useHTML = TRUE)
      ) %>%
      highcharter::hc_xAxis(title = list(text = paste("Module Membership in", module, "module"))) %>%
      highcharter::hc_yAxis(title = list(text = paste("Junction significance for trait", trait))) %>%
      highcharter::hc_tooltip(pointFormat = "junction: {point.junctionID} ",
                              crosshairs = TRUE)
    return(list(
      MMplot = MMplot,
      JunctMM = JunctMM,
      MMPvalue = MMPvalue,
      JSPvalue = JSPvalue
    ))
  }
