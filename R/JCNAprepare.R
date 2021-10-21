#' JCNAprepare: WGCNA annotation for junction expression-based network construction
#'
#' Filters samples and junctions in DJEanalize expression data object for junction co-expression network analysis.
#' @param input.type input object type. Options are DJEanalize() output object or junction read counts table.
#' @param workDir path to folder where individual junction quantification files (e.g. from STAR alignment) are located. This folder should only contain junction quantification files.
#' @param data.type One of c("sample", "matrix"), indicating whether the input files are matrices per individual samples or joint into a single expression matrix. Only used when input.type = "junction.counts".
#' @param aligner One of c("STAR", "other"), indicating the alignment tool used to produce junction quantification. Only used if data.type = "matrix".
#' #' If "other" is indicated, files in workDir path should contain junction IDs (with the format chr:start:end:strand) in the first column read counts in the second column.
#' @param analize.out output object from DJEanalize(). Only used if input.type = "DJEanalize.out".
#' @param Group1 vector or factor specifying basic control sample names. Only used if input.type = "DJEanalize.out".
#' @param gtf Reference transcriptome in genecode gtf file format. Used to define gene ID for junctions if input.type = "junction.counts".
#' @param minMean numeric: minimum of read count mean per junction. Only used if input.type = "junction.counts".
#' @param maxMean numeric: maximum of read count mean per junction. Only used if input.type = "junction.counts".
#' @param minVar numeric: minimum of read count variance per junction. Only used if input.type = "junction.counts".
#' @param maxVar numeric: maximum of read count variance per junction. Only used if input.type = "junction.counts".
#' @param calcNormFactors logical, should edgeR's library Size Normalization be calculated?. Only used if input.type = "junction.counts".
#' @param calMethod Normalization method to be used. Options are: "TMM","RLE","upperquartile" or "none" (default "none"). Only used if input.type = "junction.counts".
#' @param MVtrend.plot logical, should median/variance trend plot be displayed?. Only used if input.type = "junction.counts".
#' @param traitData a numeric vector or a matrix of external sample traits. Samples should be as rows and traits as columns.
#' @param TSdendrogram logical, should sample dendrogram plot with trait annotation be displayed?.
#' @param nThreads numeric: number of threads to allow.
#' @param abline.threshold numeric: height cut value in sample dendrogram used to remove offending sample(s).
#'
#' @return List object containing sample dendrogram, expression data matrix, trait data matrix and topology network plots
#' @examples
#' \dontrun{
#' data(DJEanlz)
#' SF <- system.file("extdata", "SF.expr.rds", package = "DJExpress")
#' SF.exp <- readRDS(SF)
#' colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))] <- paste0("TCGA_",
#' seq(1,length(colnames(DJEanlz$v.norm)[grep("TCGA", colnames(DJEanlz$v.norm))]), 1))
#' Group1 <- colnames(DJEanlz$v.norm$E)[grep("SRR", colnames(DJEanlz$v.norm$E))]
#' Jprep <- JCNAprepare(analize.out=DJEanlz, Group1 = Group1,
#' traitData = SF.exp, abline.threshold=60, input.type = "DJEanalize.out")
#' }
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
JCNAprepare <-
  function (input.type = c("DJEanalize.out", "junction.counts"),
            workDir = NULL,
            data.type = c("sample", "matrix"),
            aligner = c("STAR", "other"),
            analize.out = NULL,
            Group1=NULL,
            gtf = NULL,
            minMean = 10,
            maxMean = Inf,
            minVar = 0,
            maxVar = Inf,
            calcNormFactors = FALSE,
            calMethod = "none",
            MVtrend.plot = FALSE,
            traitData = NULL,
            TSdendrogram = FALSE,
            nThreads = 2,
            abline.threshold = NULL)
  {
    if (input.type == "DJEanalize.out") {
      if (is.null(Group1)) {
        stop("Group1 not provided")
      }
      options(stringsAsFactors = FALSE)
      datExpr0 = analize.out$v.norm
      datExpr0 = datExpr0[,-c(match(Group1, colnames(datExpr0)))]
      datExpr0 = as.data.frame(t(datExpr0$E))
      WGCNA::enableWGCNAThreads(nThreads = nThreads)
      gsg = goodSamplesJunct(datExpr0, verbose = 3)
      gsg$allOK
      if (!gsg$allOK)
      {
        if (sum(!gsg$goodGenes) > 0)
          dynamicTreeCut::printFlush(paste("Removing junctions:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))

        if (sum(!gsg$goodSamples) > 0)
          dynamicTreeCut::printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))

        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
      }
      sampleTree = fastcluster::hclust(stats::dist(datExpr0), method = "average")
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
      if (is.null(abline.threshold)) {
        sample.den <- grDevices::recordPlot()
      } else{
        graphics::abline(h = abline.threshold, col = "red")
        sample.den <- grDevices::recordPlot()
        clust = WGCNA::cutreeStatic(sampleTree, cutHeight = abline.threshold, minSize = 10)
        table(clust)
        keepSamples = (clust == 1)
        datExpr = datExpr0[keepSamples,]
        if (!is.null(traitData)) {
          traitData = traitData[keepSamples,]
        }
      }
      if (TSdendrogram == TRUE) {
        if (is.null(traitData)) {
          stop("TSdendrogram set TRUE but no traitData provided")
        }
        traitColors = WGCNA::numbers2colors(traitData, signed = FALSE)
        sampleTree2 = fastcluster::hclust(stats::dist(datExpr), method = "average")
        WGCNA::plotDendroAndColors(
          sampleTree2,
          traitColors,
          groupLabels = names(traitData),
          dendroLabels = FALSE,
          main = "Sample dendrogram and trait heatmap"
        )
        TSamp.den <- grDevices::recordPlot()

      }
      WGCNA::disableWGCNAThreads()
      options(warn = -1)
      options(stringsAsFactors = FALSE)
      powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
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
      if (TSdendrogram == TRUE) {
        return(
          list(
            sample.den = sample.den,
            TSamp.den = TSamp.den,
            datExpr = datExpr,
            datTraits = traitData,
            NetTop = NetTop,
            sft =sft
          )
        )
      } else{
        return(
          list(
            sample.den = sample.den,
            datExpr = datExpr,
            datTraits = traitData,
            NetTop = NetTop,
            sft =sft
          )
        )

      }
    } else{
      if (is.null(workDir)) {
        stop("path to junction quantification files not provided")
      }
      if (is.null(gtf)) {
        stop("path to transcriptome annotation file not provided")
      }
      if (!is.null(workDir)) {
        if (!dir.exists(workDir))
          stop("folder with junction quantification not found")
        x <- NULL
        files = list.files(workDir)
        if (aligner == "STAR") {
          i <- 1
          junctions <-
            readr::read_tsv(paste0(workDir, "/", files[i]),
                            comment = "",
                            col_names = FALSE)
          junctions <-
            as.data.frame(cbind(
              junction_id = paste(
                paste0(junctions$X1),
                junctions$X2,
                junctions$X3,
                junctions$X4,
                sep = ":"
              ),
              junctions[, 7, drop = FALSE]
            ))
          junctions$junction_id <-
            as.character(junctions$junction_id)
          colnames(junctions)[2] <- files[i]
          for (i in 2:length(files)) {
            added <-
              readr::read_tsv(paste0(workDir, "/", files[i]),
                              comment = "",
                              col_names = FALSE)
            added <-
              as.data.frame(cbind(
                junction_id = paste(paste0(added$X1), added$X2, added$X3,
                                    added$X4, sep = ":"),
                added[, 7, drop = FALSE]
              ))
            added$junction_id <- as.character(added$junction_id)
            junctions <-
              dplyr::full_join(junctions, added, by = c("junction_id"))
            colnames(junctions)[i + 1] <-
              files[i]
          }
          rownames(junctions) <- junctions$junction_id
          junctions <- junctions[, -c(1)]
          for (i in 1:ncol(junctions)) {
            junctions[, i][which(junctions[, i] < min.expressed |
                                   is.na(junctions[, i]))] = 0
          }
        } else{
          i <- 1
          junctions <-
            readr::read_tsv(paste0(workDir, "/", files[i]),
                            comment = "",
                            col_names = FALSE)
          junctions <-
            as.data.frame(junctions)
          colnames(junctions)[1] <- "junction_id"
          junctions$junction_id <-
            as.character(junctions$junction_id)
          colnames(junctions)[2] <- files[i]
          for (i in 2:length(files)) {
            added <-
              readr::read_tsv(paste0(workDir, "/", files[i]),
                              comment = "",
                              col_names = FALSE)
            added <-
              as.data.frame(added)
            colnames(added)[1] <- "junction_id"
            added$junction_id <- as.character(added$junction_id)
            junctions <-
              dplyr::full_join(junctions, added, by = c("junction_id"))
            colnames(junctions)[i + 1] <-
              files[i]
          }
          rownames(junctions) <- junctions$junction_id
          junctions <- junctions[, -c(1)]
          for (i in 1:ncol(junctions)) {
            junctions[, i][which(junctions[, i] < min.expressed |
                                   is.na(junctions[, i]))] = 0
          }
        }
      } else{
        junctions <-
          readr::read_tsv(paste0(workDir, "/", files[i]),
                          comment = "",
                          col_names = FALSE)
        junctions <-
          as.data.frame(junctions)
        colnames(junctions)[1] <- "junction_id"
        if (class(junctions$junction_id) == "numeric") {
          stop(
            "First column in matrix is numeric. Please ensure that the first column in your expression matrix correspond to junction IDs"
          )
        }
        junctions$junction_id <- as.character(junctions$junction_id)
        rownames(junctions) <- junctions$junction_id
        junctions <- junctions[, -c(1)]
        for (i in 1:ncol(junctions)) {
          junctions[, i][which(junctions[, i] < min.expressed |
                                 is.na(junctions[, i]))] = 0
        }
      }
      df <- data.frame(x = rownames(junctions))
      junct_coord <- df %>%
        tidyr::separate(x, c("chr", "start", "end", "strand"), ":")
      junct_coord$strand[which(junct_coord$strand == 0)] = ""
      junct_coord$strand[which(junct_coord$strand == 1)] = "+"
      junct_coord$strand[which(junct_coord$strand == 2)] = "-"
      junct_coord$junction_id = rownames(junctions)
      junctions <-
        junctions[which(junct_coord$strand == "+" |
                          junct_coord$strand == "-"), ]
      junct_coord <-
        junct_coord[which(junct_coord$strand == "+" |
                            junct_coord$strand == "-"), ]
      options(warn = -1)
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
      matched.ids <- data.table::chmatch(rownames(junctions),
                                         junct.genes$junction_id, nomatch =
                                           NA_integer_)
      junct.genes <-
        as.data.frame(cbind(junctionID = rownames(junctions),
                            junct.genes[matched.ids, ]))
      junctionID.2 <-
        junct.genes$junctionID[which(is.na(junct.genes$gene_id))]
      junct_coord.2 <-
        junct_coord[which(is.na(junct.genes$gene_id)), ]
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
      junct.genes3 <-
        junct.genes2[!(
          duplicated(junct.genes2$junction_id) |
            duplicated(junct.genes2$junction_id, fromLast =
                         TRUE)
        ), ]
      rownames(junct.genes3) = NULL
      matched.ids.2 <- data.table::chmatch(rownames(junctions),
                                           junct.genes3$junction_id, nomatch =
                                             NA_integer_)
      junct.genes3 <-
        as.data.frame(cbind(junctionID = rownames(junctions),
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
      junctions.2 <-
        junctions[c(which(!is.na(junct.genes$gene_id))), ]
      featureID <- rownames(junctions.2)
      groupID <-
        junct.genes$gene_id[c(which(!is.na(junct.genes$gene_id)))]
      junctions.2Mean <- apply(junctions.2, 1, mean)
      junctions.2Var <- apply(junctions.2, 1, stats::var)
      varMeanFilter <-
        junctions.2Mean >= minMean & junctions.2Mean <=
        maxMean &
        junctions.2Var >= minVar & junctions.2Var <= maxVar
      isFromGroup1 <- colnames(junctions.2) %in% Group1
      design       <- cbind(1, ifelse(isFromGroup1, 0, 1))
      junctExpr = junctions.2[varMeanFilter,]
      featureID = featureID[varMeanFilter]
      groupID = groupID[varMeanFilter]
      if (calcNormFactors) {
        dge = edgeR::calcNormFactors(dge, method = calMethod)
      }
      if (MVtrend.plot) {
        datExpr0 <-
          limma::voom(dge, design, plot = TRUE, normalize.method = normalize.method)
      } else{
        datExpr0 <-
          limma::voom(dge, design, plot = FALSE, normalize.method = normalize.method)
      }
      options(stringsAsFactors = FALSE)
      WGCNA::enableWGCNAThreads(nThreads = nThreads)
      gsg = goodSamplesJunct(datExpr0, verbose = 3)
      gsg$allOK
      if (!gsg$allOK)
      {
        if (sum(!gsg$goodGenes) > 0)
          dynamicTreeCut::printFlush(paste("Removing junctions:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))

        if (sum(!gsg$goodSamples) > 0)
          dynamicTreeCut::printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))

        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
      }
      sampleTree = fastcluster::hclust(stats::dist(datExpr0), method = "average")
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
      if (is.null(abline.threshold)) {
        sample.den <- grDevices::recordPlot()
      } else{
        graphics::abline(h = abline.threshold, col = "red")
        sample.den <- grDevices::recordPlot()
        clust = WGCNA::cutreeStatic(sampleTree, cutHeight = abline.threshold, minSize = 10)
        table(clust)
        keepSamples = (clust == 1)
        datExpr = datExpr0[keepSamples,]
        if (!is.null(traitData)) {
          traitData = traitData[keepSamples,]
        }
      }
      if (TSdendrogram == TRUE) {
        if (is.null(traitData)) {
          stop("TSdendrogram set TRUE but no traitData provided")
        }
        traitColors = WGCNA::numbers2colors(traitData, signed = FALSE)
        sampleTree2 = fastcluster::hclust(stats::dist(datExpr), method = "average")
        WGCNA::plotDendroAndColors(
          sampleTree2,
          traitColors,
          groupLabels = names(traitData),
          dendroLabels = FALSE,
          main = "Sample dendrogram and trait heatmap"
        )
        TSamp.den <- grDevices::recordPlot()
      }
      WGCNA::disableWGCNAThreads()
      options(warn = -1)
      options(stringsAsFactors = FALSE)
      powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
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
      if (TSdendrogram == TRUE) {
        return(
          list(
            sample.den = sample.den,
            TSamp.den = TSamp.den,
            datExpr = datExpr,
            datTraits = traitData,
            NetTop = NetTop,
            sft = sft
          )
        )

      } else{
        return(
          list(
            sample.den = sample.den,
            datExpr = datExpr,
            datTraits = traitData,
            NetTop = NetTop,
            sft = sft
          )
        )

      }
    }
  }
