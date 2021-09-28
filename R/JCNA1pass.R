#' JCNA1pass: 1-pass junction co-expression network analysis.
#'
#' As a wrapper of the blockwiseModules function in WGCNA package, this function performs automatic network construction and module on junction expression datasets in a block-wise manner.
#' Detected modules are not corrected for gene-based trait association effect.
#' @param prepare.out output object from JCNAprepare()
#' @param nThreads numeric: number of threads to allow.
#' @param networkType network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See adjacency.
#' @param cor.method  a character string specifying the method to be used for correlation as in WGCNA. Options are "bicor", "pearson", "kendall" or "spearman".
#' @param replaceMissingAdjacencies logical: should missing values in the calculation of adjacency be replaced by 0?
#' @param deepSplit integer value between 0 and 4. Provides a simplified control over how sensitive module detection should be to module splitting, with 0 least and 4 most sensitive. See cutreeDynamic for more details.
#' @param weights optional observation weights in the same format (and dimensions) as datExpr. These weights are used in correlation calculation.
#' @param detectCutHeight dendrogram cut height for module detection. See cutreeDynamic for more details.
#' @param minModuleSize minimum module size for module detection. See cutreeDynamic for more details.
#' @param blocks if input blocks was given, its copy; otherwise a vector of length equal number of genes giving the block label for each gene. Note that block labels are not necessarily sorted in the order in which the blocks were processed (since we do not require this for the input blocks). See blockOrder below.
#' @param maxBlockSize integer giving maximum block size for module detection. Ignored if blocks above is non-NULL. Otherwise, if the number of genes in datExpr exceeds maxBlockSize, genes will be pre-clustered into blocks whose size should not exceed maxBlockSize.
#' @param blockSizePenaltyPower number specifying how strongly blocks should be penalized for exceeding the maximum size. Set to a lrge number or Inf if not exceeding maximum block size is very important.
#' @param nPreclusteringCenters number of centers for pre-clustering. Larger numbers typically results in better but slower pre-clustering.
#' @param randomSeed integer to be used as seed for the random number generator before the function starts. If a current seed exists, it is saved and restored upon exit. If NULL is given, the function will not save and restore the seed.
#' @param checkMissingData logical: should data be checked for excessive numbers of missing entries in genes and samples, and for genes with zero variance? See details.
#' @param TOMType one of "none", "unsigned", "signed", "signed Nowick", "unsigned 2", "signed 2" and "signed Nowick 2". If "none", adjacency will be used for clustering. See TOMsimilarityFromExpr for details.
#' @param TOMDenom a character string specifying the TOM variant to be used. Recognized values are "min" giving the standard TOM described in Zhang and Horvath (2005), and "mean" in which the min function in the denominator is replaced by mean. The "mean" may produce better results but at this time should be considered experimental.
#' @param suppressTOMForZeroAdjacencies logical: should TOM be set to zero for zero adjacencies?
#' @param suppressNegativeTOM logical: should the result be set to zero when negative? Negative TOM values can occur when TOMType is "signed Nowick".
#' @param maxCoreScatter maximum scatter of the core for a branch to be a cluster, given as the fraction of cutHeight relative to the 5th percentile of joining heights. See cutreeDynamic for more details.
#' @param minGap minimum cluster gap given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. See cutreeDynamic for more details.
#' @param maxAbsCoreScatter maximum scatter of the core for a branch to be a cluster given as absolute heights. If given, overrides maxCoreScatter. See cutreeDynamic for more details.
#' @param minAbsGap minimum cluster gap given as absolute height difference. If given, overrides minGap. See cutreeDynamic for more details.
#' @param minSplitHeight minimum split height given as the fraction of the difference between cutHeight and the 5th percentile of joining heights. Branches merging below this height will automatically be merged. Defaults to zero but is used only if minAbsSplitHeight below is NULL.
#' @param minAbsSplitHeight minimum split height given as an absolute height. Branches merging below this height will automatically be merged. If not given (default), will be determined from minSplitHeight above.
#' @param stabilityLabels optional matrix of cluster labels that are to be used for calculating branch dissimilarity based on split stability. The number of rows must equal the number of genes in multiExpr; the number of columns (clusterings) is arbitrary. See branchSplitFromStabilityLabels for details.
#' @param useBranchEigennodeDissim logical: should branch eigennode (eigenjunction) dissimilarity be considered when merging branches in Dynamic Tree Cut?
#' @param minBranchEigennodeDissim minimum consensus branch eigennode (eigenjunction) dissimilarity for branches to be considerd separate. The branch eigennode dissimilarity in individual sets is simly 1-correlation of the eigennodes; the consensus is defined as quantile with probability consensusQuantile.
#' @param stabilityCriterion one of c("Individual fraction", "Common fraction"), indicating which method for assessing stability similarity of two branches should be used. We recommend "Individual fraction" which appears to perform better; the "Common fraction" method is provided for backward compatibility since it was the (only) method available prior to WGCNA version 1.60.
#' @param minStabilityDissim minimum stability dissimilarity criterion for two branches to be considered separate. Should be a number between 0 (essentially no dissimilarity required) and 1 (perfect dissimilarity or distinguishability based on stabilityLabels). See branchSplitFromStabilityLabels for details.
#' @param pamStage logical. If TRUE, the second (PAM-like) stage of module detection will be performed. See cutreeDynamic for more details.
#' @param pamRespectsDendro logical, only used when pamStage is TRUE. If TRUE, the PAM stage will respect the dendrogram in the sense an object can be PAM-assigned only to clusters that lie below it on the branch that the object is merged into. See cutreeDynamic for more details.
#' @param reassignThreshold p-value ratio threshold for reassigning junctions between modules. See Details.
#' @param minCoreKME a number between 0 and 1. If a detected module does not have at least minModuleKMESize genes with eigenjunction connectivity at least minCoreKME, the module is disbanded (its genes are unlabeled and returned to the pool of genes waiting for mofule detection).
#' @param minCoreKMESize see minCoreKME above.
#' @param minKMEtoStay junctions whose eigenjunction connectivity to their module eigenjunction is lower than minKMEtoStay are removed from the module.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param impute logical: should imputation be used for module eigengene calculation? See moduleEigengenes for more details.
#' @param trapErrors logical: should errors in calculations be trapped?
#' @param numericLabels logical: should the returned modules be labeled by colors (FALSE), or by numbers (TRUE)?
#' @param corType character string specifying the correlation to be used. Allowed values are (unique abbreviations of) "pearson" and "bicor", corresponding to Pearson and bidweight midcorrelation, respectively. Missing values are handled using the pairwise.complete.obs option.
#' @param maxPOutliers only used for corType=="bicor". Specifies the maximum percentile of data that can be considered outliers on either side of the median separately. For each side of the median, if higher percentile than maxPOutliers is considered an outlier by the weight function based on 9*mad(x), the width of the weight function is increased such that the percentile of outliers on that side of the median equals maxPOutliers. Using maxPOutliers=1 will effectively disable all weight function broadening; using maxPOutliers=0 will give results that are quite similar (but not equal to) Pearson correlation.
#' @param quickCor real number between 0 and 1 that controls the handling of missing data in the calculation of correlations. See details.
#' @param pearsonFallback Specifies whether the bicor calculation, if used, should revert to Pearson when median absolute deviation (mad) is zero. Recongnized values are (abbreviations of) "none", "individual", "all". If set to "none", zero mad will result in NA for the corresponding correlation. If set to "individual", Pearson calculation will be used only for columns that have zero mad. If set to "all", the presence of a single zero mad will cause the whole variable to be treated in Pearson correlation manner (as if the corresponding robust option was set to FALSE). Has no effect for Pearson correlation. See bicor.
#' @param cosineCorrelation logical: should the cosine version of the correlation calculation be used? The cosine calculation differs from the standard one in that it does not subtract the mean.
#' @param loadTOM logical: should Topological Overlap Matrices be loaded from previously saved files (TRUE) or calculated (FALSE)? It may be useful to load previously saved TOM matrices if these have been calculated previously, since TOM calculation is often the most computationally expensive part of network construction and module identification. See saveTOMs and saveTOMFileBase below for when and how TOM files are saved, and what the file names are. If loadTOM is TRUE but the files cannot be found, or do not contain the correct TOM data, TOM will be recalculated.
#' @param saveTOMs logical: should the consensus topological overlap matrices for each block be saved and returned?
#' @param saveTOMFileBase character string containing the file name base for files containing the consensus topological overlaps. The full file names have "block.1.RData", "block.2.RData" etc. appended. These files are standard R data files and can be loaded using the load function.
#' @param useInternalMatrixAlgebra logical: should WGCNA's own, slow, matrix multiplication be used instead of R-wide BLAS? Only useful for debugging.
#' @param useCorOptionsThroughout logical: should correlation options passed to network analysis also be used in calculation of kME? Set to FALSE to reproduce results obtained with WGCNA 1.62 and older.
#' @param verbose integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no indentation, each unit adds two spaces.
#' @seealso \code{\link[WGCNA]{blockwiseModules}}
#'
#' @return A list object containing:
#'
#' net: Constructed junction network and detected modules using blockwiseModules from WGCNA
#'
#' datExpr: Junction expression dataset
#'
#' datTraits: External traits dataset
#'
#' moduleColors: Assigned module name for each junction in datExpr
#'
#' module.den: Junction dendrogram divided by blockGenes in net
#'
#' moduleTraitCor: Matrix of correlation coefficients (Module eigengenes vs traits)
#'
#' moduleTraitPvalue: Matrix of correlation P-values (Module eigengenes vs traits)
#'
#' juncModuleMembership: Matrix of correlation coefficients (Junction expression vs Module eigengenes)
#'
#' juncMMPvalue: Matrix of correlation P-values (Junction expression vs Module eigengenes)
#'
#' ModuleTrait: Heatmap with Junction Modules vs Trait associations
#'
#' Junctrait: Data frame with correlation coefficients and associated P-values (Junction expression vs traits)
#'
#' sig.cors.trait: List with significant module-trait associations (coefficient >0.2 & P-value < 0.05)
#' @examples
#' \dontrun{
#' JCNAprep <- system.file("extdata", "Jprep.rds", package = "DJExpress")
#' Jprep <- readRDS(JCNAprep)
#' J1pass <- JCNA1pass(Jprep, cor.method = "bicor", nThreads = 2)
#' }
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom WGCNA bicor
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
JCNA1pass <-
  function (prepare.out,
            nThreads = 2,
            networkType = "signed",
            cor.method = c("bicor", "pearson", "kendall", "spearman"),
            replaceMissingAdjacencies = FALSE,
            deepSplit = 2,
            weights = NULL,
            detectCutHeight = 0.995,
            minModuleSize = min(20, ncol(datExpr) / 2),
            blocks = NULL,
            maxBlockSize = 5000,
            blockSizePenaltyPower = 5,
            nPreclusteringCenters = as.integer(min(ncol(datExpr) /
                                                     20,
                                                   100 * ncol(datExpr) /
                                                     maxBlockSize)),
            randomSeed = 54321,
            checkMissingData = TRUE,
            TOMType = "signed",
            TOMDenom = "min",
            suppressTOMForZeroAdjacencies = FALSE,
            suppressNegativeTOM = FALSE,
            maxCoreScatter = NULL,
            minGap = NULL,
            maxAbsCoreScatter = NULL,
            minAbsGap = NULL,
            minSplitHeight = NULL,
            minAbsSplitHeight = NULL,
            stabilityLabels = NULL,
            useBranchEigennodeDissim = FALSE,
            minBranchEigennodeDissim = mergeCutHeight,
            stabilityCriterion = c("Individual fraction", "Common fraction"),
            minStabilityDissim = NULL,
            pamStage = TRUE,
            pamRespectsDendro = TRUE,
            reassignThreshold = 0,
            minCoreKME = 0.5,
            minCoreKMESize = minModuleSize / 3,
            minKMEtoStay = 0.3,
            mergeCutHeight = 0.25,
            impute = TRUE,
            trapErrors = FALSE,
            numericLabels = TRUE,
            corType = "bicor",
            maxPOutliers = 1,
            quickCor = 0,
            pearsonFallback = "individual",
            cosineCorrelation = FALSE,
            loadTOM = FALSE,
            saveTOMs = TRUE,
            saveTOMFileBase = "JCNA_blockwiseTOM",
            useInternalMatrixAlgebra = FALSE,
            useCorOptionsThroughout = TRUE,
            verbose = 3,
            indent = 0)
  {
    sft <- NULL
    datExpr <- prepare.out$datExpr
    print("Constructing the junction network and modules")
    selected.p <-
      max(sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq < 0.90)])
    if (selected.p == -Inf) {
      selected.p = 6
    }
    if (selected.p == 30) {
      warning("Scale-free topology fit index failed to reach power values above 0.8")
      warning("Choosing power based on the number of samples")
      if (nrow(datExpr) <= 20 & networkType == "signed") {
        selected.p = 18
      }
      if (nrow(datExpr) <= 20 &
          networkType == "unsigned" | networkType == "signed hybrid") {
        selected.p = 9
      }
      if (nrow(datExpr) > 20 &
          nrow(datExpr) <= 30 & networkType == "signed") {
        selected.p = 16
      }
      if (nrow(datExpr) > 20 &
          nrow(datExpr) <= 30 &
          networkType == "unsigned" | networkType == "signed hybrid") {
        selected.p = 8
      }
      if (nrow(datExpr) > 30 &
          nrow(datExpr) <= 40 & networkType == "signed") {
        selected.p = 14
      }
      if (nrow(datExpr) > 30 &
          nrow(datExpr) <= 40 &
          networkType == "unsigned" | networkType == "signed hybrid") {
        selected.p = 7
      }
      if (nrow(datExpr) > 40 & networkType == "signed") {
        selected.p = 12
      }
      if (nrow(datExpr) > 40 &
          networkType == "unsigned" | networkType == "signed hybrid") {
        selected.p = 6
      }
    }
    WGCNA::enableWGCNAThreads(nThreads = nThreads)
    cor <- WGCNA::cor
    if(cor.method=="bicor" | cor.method=="pearson"){
      net <- WGCNA::blockwiseModules(datExpr, corType= cor.method, power = selected.p)
    }else{
      net <- WGCNA::blockwiseModules(datExpr, corType= "bicor", power = selected.p)
    }
    cor<-stats::cor
    print(paste0(" modules identified:", as.character(length(
      table(net$colors)
    ))))
    print(table(net$colors))
    print("Label 0 is reserved for genes outside of all modules")
    moduleColors = WGCNA::labels2colors(net$colors)
    module.den <- list()
    for (i in 1:length(net$blockGenes)) {
        WGCNA::plotDendroAndColors(
          net$dendrograms[[i]],
          moduleColors[net$blockGenes[[i]]],
          "Module colors",
          dendroLabels = FALSE,
          hang = 0.03,
          addGuide = TRUE,
          guideHang = 0.05
        )
      module.den[[i]] <- grDevices::recordPlot()
      grDevices::dev.off()
    }
    moduleLabels = net$colors
    MEs = net$MEs
    geneTree = net$dendrograms[[1]]
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    print("Calculating module membership values and trait correlations")
    allTraits <- prepare.out$datTraits
    MEs0 = WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = WGCNA::orderMEs(MEs0)
    if (cor.method == "bicor"){
      moduleTraitCor = WGCNA::bicor(MEs, allTraits, use = "p")
    }else{
      moduleTraitCor = WGCNA::cor(MEs, allTraits, method = cor.method,  use = "p")
    }
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
    modNames = substring(names(MEs), 3)
    if (cor.method == "bicor"){
      juncModuleMembership = as.data.frame(WGCNA::bicor(datExpr, MEs, use = "p"))
    }else{
      juncModuleMembership = as.data.frame(WGCNA::cor(datExpr, MEs, method = cor.method,  use = "p"))
    }
    MMPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(juncModuleMembership), nSamples))
    names(juncModuleMembership) = paste("MM", modNames, sep="")
    names(MMPvalue) = paste("p.MM", modNames, sep="")
    sig.cors.trait <- list()
    for (i in 1:ncol(moduleTraitCor)) {
      sig.cor <-
        which(abs(round(moduleTraitCor[, i], 1)) >= 0.2 &
                abs(moduleTraitPvalue[, i]) <= 0.05)
      sig.cors.trait[[i]] <- sig.cor
    }
    names(sig.cors.trait) <- colnames(moduleTraitCor)
    textMatrix = paste(signif(moduleTraitCor, 2),
                       "\n(",
                       signif(moduleTraitPvalue, 1),
                       ")",
                       sep = "")
      WGCNA::labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = names(allTraits),
      yLabels = names(MEs),
      ySymbols = names(MEs),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(70),
      textMatrix = textMatrix,
      setStdMargins = FALSE,
      cex.text = 0.5,
      zlim = c(-1, 1),
      main = paste("Module-trait relationships")
    )
    ModuleTrait <- grDevices::recordPlot()
    grDevices::dev.off()

    print("Defining junction significance for each trait")
    i <- 1
    trait <- as.data.frame(allTraits[, i])
    names(trait) = colnames(allTraits)[i]
    if (cor.method == "bicor"){
      JunctraitSig = as.data.frame(WGCNA::bicor(datExpr, trait, use = "p"))
    }else{
      JunctraitSig = as.data.frame(WGCNA::cor(datExpr, trait, method = cor.method,  use = "p"))
    }
    GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(JunctraitSig), nSamples))
    names(JunctraitSig) = paste("GS.", names(trait), sep = "")
    names(GSPvalue) = paste("p.GS.", names(trait), sep = "")
    Junctrait <- as.data.frame(cbind(JunctraitSig, GSPvalue))
    for (i in 2:ncol(allTraits)) {
      trait <-  as.data.frame(allTraits[, i])
      names(trait) = colnames(allTraits)[i]
      if (cor.method == "bicor"){
        JunctraitSig = as.data.frame(WGCNA::bicor(datExpr, trait, use = "p"))
      }else{
        JunctraitSig = as.data.frame(WGCNA::cor(datExpr, trait, method = cor.method, use = "p"))
      }
      GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(JunctraitSig), nSamples))
      names(JunctraitSig) = paste("GS.", names(trait), sep = "")
      names(GSPvalue) = paste("p.GS.", names(trait), sep = "")
      trait.sig <- as.data.frame(cbind(JunctraitSig, GSPvalue))
      Junctrait <- cbind(Junctrait, trait.sig)
    }
    print("done")
    return(
      list(
        net = net,
        datExpr = datExpr,
        datTraits = allTraits,
        moduleColors = moduleColors,
        module.den = module.den,
        moduleTraitCor = moduleTraitCor,
        moduleTraitPvalue = moduleTraitPvalue,
        juncModuleMembership = juncModuleMembership,
        juncMMPvalue = MMPvalue,
        ModuleTrait = ModuleTrait,
        Junctrait = Junctrait,
        sig.cors.trait = sig.cors.trait
      )
    )
  }
