#' JCNAGenePrepare: WGCNA annotation for junction expression network construction
#'
#' Performs weighted junction co-expression network analysis using gene expression and returns output object for JCNA2pass.
#' @param pass1.out output object from JCNA1pass()
#' @param genExpr Gene expression data. A matrix (preferred) or data frame in which columns are genes and rows ar samples. NAs are allowed, but not too many. See WGCNA::checkMissingData and details.
#' @param Group1 vector or factor specifying basic control sample names
#' @param gtf Reference transcriptome in a gtf file. Used to define gene identity of junctions
#' @param nThreads numeric: number of threads to allow.
#' @param abline.threshold numeric: height cut value in sample dendrogram used to remove offending sample(s)
#' @param networkType network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See adjacency.
#' @param level numeric: confidence interval for Linear Regression of Junction vs Gene trait association value
#' @param cor.method  a character string specifying the method to be used for correlation as in WGCNA. Options are "bicor", "pearson", "kendall" or "spearman".
#'
#' @return A list object containing:
#' Gene.sampletree: Sample dendrogram based on gene expression
#' Gene.hc: Sample dendrogram with trait heatmap
#' Gene.NetTop: Analysis of network topology for various soft-thresholding powers produced by pickSoftThreshold in WGCNA.
#' Left panel shows scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis).
#' The right panel displays mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis).
#' Gene.net: Constructed gene network and detected modules using blockwiseModules from WGCNA
#' Gene.ModuleTrait: Heatmap with Gene Modules vs Trait associations
#' Genetrait: Data frame with correlation coefficients and associated P-values (Gene expression vs traits)
#' GeneToJunct: Data frame with junction coordinates and correspondent gene
#' JunctGeneTrait: List of data frames with integrated information about junction and gene associations per trait
#' @examples
#' \dontrun{
#' gtf0 <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
#' gexp <- system.file("extdata", "genExpr.rds", package = "DJExpress")
#' filgenExpr <- readRDS(gexp)
#' J1pass0 <- system.file("extdata", "J1pass.rds", package = "DJExpress")
#' J1pass <- readRDS(J1pass0)
#' JgPrep <- JCNAgenePrepare(pass1.out = J1pass, genExpr = filgenExpr,
#'                           gtf=gtf0, networkType = "unsigned", cor.method = "bicor")
#' }
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom WGCNA bicor
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
JCNAgenePrepare <-
  function (pass1.out,
            genExpr,
            Group1 = NULL,
            gtf,
            nThreads = 2,
            abline.threshold = NULL,
            networkType = "unsigned",
            level = 0.95,
            cor.method = c("bicor", "pearson", "kendall", "spearman")) {
    net <- pass1.out$net
    MEs <- pass1.out$net$MEs
    Junctrait <- pass1.out$Junctrait
    traitData <- pass1.out$datTraits
    print("Calculating co-expression networks from gene expression")
    options(stringsAsFactors = FALSE)
    genExpr0 = genExpr
    if (!is.null(Group1)) {
      genExpr0 = genExpr0[, -c(match(Group1, colnames(genExpr0)))]
      genExpr0 = as.data.frame(t(genExpr0$E))
    }
    genExpr0 <- as.data.frame(t(genExpr0))
    match.samples <-
      match(rownames(traitData), rownames(genExpr0))
    genExpr <- genExpr0[match.samples[which(!is.na(match.samples))], ]
    traitData <- traitData[which(!is.na(match.samples)), ]
    WGCNA::enableWGCNAThreads(nThreads = nThreads)
    gsg = goodSamplesJunct(genExpr, verbose = 3)
    gsg$allOK
    if (!gsg$allOK)
    {
      if (sum(!gsg$goodGenes) > 0)
        dynamicTreeCut::printFlush(paste("Removing genes:", paste(names(genExpr)[!gsg$goodGenes], collapse = ", ")))
      if (sum(!gsg$goodSamples) > 0)
        dynamicTreeCut::printFlush(paste("Removing samples:", paste(rownames(genExpr)[!gsg$goodSamples], collapse = ", ")))
        genExpr = genExpr[gsg$goodSamples, gsg$goodGenes]
    }
    sampleTree = fastcluster::hclust(stats::dist(genExpr), method = "average")
    WGCNA::sizeGrWindow(12, 9)
    graphics::par(cex = 0.6)
    graphics::par(mar = c(0, 4, 2, 0))
    plot(
      sampleTree,
      main = "Sample clustering to detect outliers",
      sub = "",
      xlab = "",
      cex.lab = 1.5,
      cex.axis = 1.5,
      cex.main = 2
    )
    warning("consider using a abline.threshold value to remove outlier samples")
    if (!is.null(abline.threshold)) {
      graphics::abline(h = abline.threshold, col = "red")
      gene.sampletree <- grDevices::recordPlot()
      grDevices::dev.off()
      clust = WGCNA::cutreeStatic(sampleTree, cutHeight = abline.threshold, minSize = 10)
      table(clust)
      keepSamples = (clust == 1)
      datExpr = genExpr[keepSamples,]
      traitData = traitData[keepSamples,]
    }else{
      WGCNA::sizeGrWindow(12, 9)
      graphics::par(cex = 0.6)
      graphics::par(mar = c(0, 4, 2, 0))
      plot(
        sampleTree,
        main = "Sample clustering to detect outliers",
        sub = "",
        xlab = "",
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 2
      )
      gene.sampletree <- grDevices::recordPlot()
      grDevices::dev.off()
      datExpr = genExpr
      traitData = traitData
    }
    sampleTree2 = fastcluster::hclust(stats::dist(datExpr), method = "average")
    traitColors = WGCNA::numbers2colors(traitData, signed = FALSE)

    WGCNA::plotDendroAndColors(
      sampleTree2,
      traitColors,
      Group1Labels = names(traitData),
      dendroLabels = FALSE,
      main = "Sample dendrogram and trait heatmap"
    )
    hc <- grDevices::recordPlot()
    grDevices::dev.off()
    WGCNA::disableWGCNAThreads()
    options(warn = -1)
    options(stringsAsFactors = FALSE)
    powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
    sft = WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    WGCNA::sizeGrWindow(9, 5)
    graphics::par(mfrow = c(1, 2))
    cex1 = 0.9
    plot(
      sft$fitIndices[, 1],
      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      xlab = "Soft Threshold (power)",
      ylab = "Scale Free Topology Model Fit,signed R^2",
      type = "n",
      main = paste("Scale independence")
    )
    graphics::text(
      sft$fitIndices[, 1],
      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers,
      cex = cex1,
      col = "red"
    )
    graphics::abline(h = 0.90, col = "red")
    plot(
      sft$fitIndices[, 1],
      sft$fitIndices[, 5],
      xlab = "Soft Threshold (power)",
      ylab = "Mean Connectivity",
      type = "n",
      main = paste("Mean connectivity")
    )
    graphics::text(
      sft$fitIndices[, 1],
      sft$fitIndices[, 5],
      labels = powers,
      cex = cex1,
      col = "red"
    )
    NetTop <- grDevices::recordPlot()
    grDevices::dev.off()
    print("Constructing the gene network and modules")
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
    moduleColors = WGCNA::labels2colors(net$colors)
    moduleLabels = net$colors
    MEs = net$MEs
    geneTree = net$dendrograms
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    MEs0 = WGCNA::moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = WGCNA::orderMEs(MEs0)
    if (cor.method == "bicor"){
      moduleTraitCor = WGCNA::bicor(MEs, traitData, use = "p")
    }else{
      moduleTraitCor = WGCNA::cor(MEs, traitData, method = cor.method, use = "p")
    }
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
    textMatrix = paste(signif(moduleTraitCor, 2),
                       "\n(",
                       signif(moduleTraitPvalue, 1),
                       ")",
                       sep = "")
    WGCNA::labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = names(traitData),
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
    print("Defining gene significance for each trait")
    i <- 1
    trait <- as.data.frame(traitData[, i])
    names(trait) = colnames(traitData)[i]
    if (cor.method == "bicor"){
      GenetraitSig = as.data.frame(WGCNA::bicor(datExpr, trait, use = "p"))
    }else{
      GenetraitSig = as.data.frame(WGCNA::cor(datExpr, trait, method = cor.method, use = "p"))
    }
    GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(GenetraitSig), nSamples))
    names(GenetraitSig) = paste("GS.", names(trait), sep = "")
    names(GSPvalue) = paste("p.GS.", names(trait), sep = "")
    Genetrait <- as.data.frame(cbind(GenetraitSig, GSPvalue))
    for (i in 2:ncol(traitData)) {
      trait <-  as.data.frame(traitData[, i])
      names(trait) = colnames(traitData)[i]
      if (cor.method == "bicor"){
        GenetraitSig = as.data.frame(WGCNA::bicor(datExpr, trait, use = "p"))
      }else{
        GenetraitSig = as.data.frame(WGCNA::cor(datExpr, trait, method = cor.method, use = "p"))
      }
      GSPvalue = as.data.frame(WGCNA::corPvalueStudent(as.matrix(GenetraitSig), nSamples))
      names(GenetraitSig) = paste("GS.", names(trait), sep = "")
      names(GSPvalue) = paste("p.GS.", names(trait), sep = "")
      trait.sig <- as.data.frame(cbind(GenetraitSig, GSPvalue))
      Genetrait <- cbind(Genetrait, trait.sig)
    }
    JdatExpr <- pass1.out$datExpr
    df <- data.frame(x = colnames(JdatExpr))
    junct_coord <-
      df %>% tidyr::separate(x, c("chr", "start", "end", "strand"), ":")
    junct_coord$strand[which(junct_coord$strand == 0)] = ""
    junct_coord$strand[which(junct_coord$strand == 1)] = "+"
    junct_coord$strand[which(junct_coord$strand == 2)] = "-"
    junct_coord$junction_id = colnames(JdatExpr)
    if (unlist(base::strsplit(gtf, "[.]"))[length(unlist(base::strsplit(gtf, "[.]")))] ==
        "gz") {
      ggtf <-
        utils::read.delim(gzfile(gtf), header = FALSE, comment.char = "#")
    } else{
      ggtf <- utils::read.delim(gtf, header = FALSE, comment.char = "#")
    }
    if (length(grep("gene", ggtf$V3)) == 0)
      stop(
        "gtf doesn't seem to have the right format. Please provide a gtf with <gene> as available feature in third column (e.g. genecode)."
      )
    ggtf_genes <- ggtf[grep("gene", ggtf$V3), ]
    ggtf_genes.id <-
      stringr::str_match(ggtf_genes$V9, "gene_name (.*?);")
    gtf_coord = as.data.frame(ggtf_genes[, c(1, 4, 5, 7)])
    gtf_coord$gene_id = ggtf_genes.id[, 2]
    if (grepl("_", gtf_coord$gene_id[1])) {
      gtf_coord$gene_id <- sub("\\_.*", "", gtf_coord$gene_id)
    }
    colnames(gtf_coord) = c("chr", "start", "end", "strand", "gene_id")
    gtf_coord$start = as.numeric(as.vector(gtf_coord$start))
    gtf_coord$end = as.numeric(as.vector(gtf_coord$end))
    junct_coord$start = as.numeric(as.vector(junct_coord$start))
    junct_coord$end = as.numeric(as.vector(junct_coord$end))
    gtf_grange <- with(gtf_coord,
                       GenomicRanges::GRanges(chr, IRanges::IRanges(start, end, names = gene_id), strand))
    junct_grange <- with(junct_coord,
                         GenomicRanges::GRanges(chr, IRanges::IRanges(start, end, names = junction_id), strand))
    overlaps1 <-
      as.data.frame(
        IRanges::findOverlaps(
          junct_grange,
          gtf_grange,
          type = "within",
          ignore.strand = FALSE
        )
      )
    junct.pos = as.integer(overlaps1$queryHits)
    gtf.pos = as.integer(overlaps1$subjectHits)
    junct.genes = cbind(junct_coord[junct.pos, ], gtf_coord[gtf.pos, ])
    rownames(junct.genes) = NULL
    junct.genes <- junct.genes[!(
      duplicated(junct.genes$junction_id) |
        duplicated(junct.genes$junction_id, fromLast =
                     TRUE)
    ), ]
    matched.ids <- data.table::chmatch(colnames(JdatExpr),
                                       junct.genes$junction_id, nomatch =
                                         NA_integer_)
    junct.genes <- as.data.frame(cbind(junctionID = colnames(JdatExpr),
                                       junct.genes[matched.ids, ]))
    junctionID.2 <-
      junct.genes$junctionID[which(is.na(junct.genes$gene_id))]
    junct_coord.2 <- junct_coord[which(is.na(junct.genes$gene_id)), ]
    junct_coord.2$start = as.numeric(as.vector(junct_coord.2$start))
    junct_coord.2$end = as.numeric(as.vector(junct_coord.2$end))

    junct_grange.2 <- with(junct_coord.2,
                           GenomicRanges::GRanges(chr, IRanges::IRanges(start, end, names = junction_id), strand))
    overlaps2 <-
      as.data.frame(
        IRanges::findOverlaps(
          junct_grange.2,
          gtf_grange,
          type = "within",
          ignore.strand = FALSE
        )
      )
    junct.pos2 = as.integer(overlaps2$queryHits)
    gtf.pos2 = as.integer(overlaps2$subjectHits)
    junct.genes2 = cbind(junct_coord.2[junct.pos2, ], gtf_coord[gtf.pos2, ])
    rownames(junct.genes2) = NULL
    # Remove ambiguous annotation
    junct.genes3 <-
      junct.genes2[!(
        duplicated(junct.genes2$junction_id) |
          duplicated(junct.genes2$junction_id, fromLast =
                       TRUE)
      ), ]
    rownames(junct.genes3) = NULL
    matched.ids.2 <- data.table::chmatch(colnames(JdatExpr),
                                         junct.genes3$junction_id, nomatch =
                                           NA_integer_)
    junct.genes3 <- as.data.frame(cbind(junctionID = colnames(JdatExpr),
                                        junct.genes3[matched.ids.2, ]))
    rownames(junct.genes3) = NULL
    colnames(junct.genes) <-
      c(
        "featureID",
        "chr",
        "start",
        "end",
        "strand",
        "junction_id",
        "chr_val",
        "start_val",
        "end_val",
        "strand_val",
        "gene_id"
      )
    colnames(junct.genes3) <-
      c(
        "featureID",
        "chr",
        "start",
        "end",
        "strand",
        "junction_id",
        "chr_val",
        "start_val",
        "end_val",
        "strand_val",
        "gene_id"
      )
    missing.juct <- which(!is.na(junct.genes3$gene_id))
    junct.genes[missing.juct, ] <- junct.genes3[missing.juct, ]
    rownames(junct.genes) <- NULL
    JunctMM <- as.data.frame(cbind(junct.genes, pass1.out$Junctrait))
    print("Defining gene-Trait vs junction-Trait correlation")
    sig.cors.trait <- J1pass$sig.cors.trait
    JunctGeneTrait <- list()
    for (i in 1:length(sig.cors.trait)) {
      trait <- names(sig.cors.trait)[i]
      sig.junct0 <- JunctMM[, c(11, grep(trait, colnames(JunctMM)))]
      rownames(sig.junct0) <- JunctMM$featureID
      match.gene <-
        match(sig.junct0$gene_id, rownames(Genetrait))
      sig.gene0 <-
        Genetrait[match.gene, grep(trait, colnames(Genetrait))]
      sig.junct0[, 2] <- as.numeric(as.character(sig.junct0[, 2]))
      sig.gene0[, 1] <- as.numeric(as.character(sig.gene0[, 1]))
      colnames(sig.junct0)[c(2, 3)] <- c("sig.junct", "p.sig.junct")
      colnames(sig.gene0) <- c("sig.gene", "p.sig.gene")
      joint.sig <- as.data.frame(cbind(sig.junct0, sig.gene0))
      model1 <- stats::lm(sig.gene ~ sig.junct, data = joint.sig)
      xRange = data.frame(sig.junct = seq(
        min(joint.sig$sig.junct, na.rm = T),
        max(joint.sig$sig.junct, na.rm = T),
        0.00001
      ))
      pred_interval <-
        stats::predict(
          model1,
          newdata = data.frame(sig.junct = xRange),
          interval = "prediction",
          level = level
        )
      modelConfInt <- stats::predict(model1,
                                     joint.sig,
                                     level = level,
                                     interval = "prediction")
      insideInterval <-
        modelConfInt[, 'lwr'] < joint.sig[['sig.gene']] &
        joint.sig[['sig.gene']] < modelConfInt[, 'upr']
      joint.sig$gene.association <- NA
      joint.sig$gene.association[which(insideInterval == "TRUE")] <-
        "gene_dependent"
      joint.sig$gene.association[which(insideInterval == "FALSE")] <-
        "gene_independent"
      joint.sig$gene.association[which(is.na(joint.sig$sig.gene))] <-
        "no_geneInfo"
      JunctGeneTrait[[i]] <- joint.sig
    }
    names(JunctGeneTrait) <- colnames(traitData)
    return(
      list(
        Gene.sampletree = gene.sampletree,
        Gene.hc = hc,
        Gene.NetTop = NetTop,
        Gene.net = net,
        Gene.ModuleTrait = ModuleTrait,
        Genetrait = Genetrait,
        GeneToJunct = junct.genes,
        JunctGeneTrait = JunctGeneTrait
      )
    )
  }
