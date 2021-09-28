#' DJEannotate: Gene annotation for junctions used in differential expression analysis
#'
#' Matches junction ID in expression matrix with respective gene name and returns output object for DJEprepare.
#' @param import.out output object from DJEimport()
#' @param gtf Reference transcriptome in genecode gtf file format. Used to define gene ID for junctions
#'
#' @return List object containing junction expression and assigned gene IDs
#' @examples
#' DJEimp <- system.file("extdata", "DJEimp.rds", package = "DJExpress")
#' imp.out <- readRDS(DJEimp)
#' gtf0 <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
#' ann.out <- DJEannotate(imp.out, gtf0)
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEannotate <- function(import.out, gtf) {
  options(warn = -1)
  junctions <- import.out$quant
  junct_coord <- import.out$coord
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
  junct.genes <- as.data.frame(cbind(junctionID = rownames(junctions),
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
  junctions.2 <- junctions[c(which(!is.na(junct.genes$gene_id))), ]
  featureID <- rownames(junctions.2)
  groupID <-
    junct.genes$gene_id[c(which(!is.na(junct.genes$gene_id)))]
  return(list(
    quant.annotated = junctions.2,
    featureID = featureID,
    groupID = groupID
  ))
}
