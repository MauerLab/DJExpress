#' JAn: Gene correspondence for splice junctions based on transcriptome annotation
#'
#' Identifies splice junctions based on the provided transcriptome annotation (gtf file)
#' @param gtf Reference transcriptome in a gtf file. Used to define gene identity of junctions
#'
#' @return Data frame with junction IDs and Associated gene ID
#' @examples
#' \dontrun{
#' gtf <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
#' Jan.out <- JAn(gtf)
#' }
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
JAn <- function (gtf) {
  if (unlist(base::strsplit(gtf, "[.]"))[length(unlist(base::strsplit(gtf, "[.]")))] ==
      "gz") {
    ggtf <-
      utils::read.delim(gzfile(gtf), header = FALSE, comment.char = "#")
  } else{
    ggtf <- utils::read.delim(gtf, header = FALSE, comment.char = "#")
  }
  if (length(grep("gene", ggtf$V3))==0)
    stop(
      "gtf doesn't seem to have the right format. Please provide a gtf with <gene> as available feature in third column (e.g. genecode)."
    )
  gene.ids <-
    stringr::str_match(ggtf$V9, "gene_id \\s*(.*?)\\s*;")[, 2]
  if (grepl("_", gene.ids[1])) {
    gene.ids <- sub("\\_.*", "", gene.ids)
  }
  gene.ids.u <- unique(gene.ids)
  junct.annot <- c()
  for (i in 1:length(gene.ids.u)) {
    this.gene <- which(gene.ids == gene.ids.u[i])
    ggtf.gene <- ggtf[this.gene, ]
    if (is.na(match("transcript", ggtf.gene$V3))) {
      if (ggtf.gene$V7[1] == "+") {
        start <- ggtf.gene$V5 + 1
        start <- start[-c(length(start))]
        end <- ggtf.gene$V4 - 1
        end <- end[-c(1)]
        junctions <-
          as.data.frame(cbind(
            gene = gene.ids.u[i],
            junction = paste(ggtf.gene$V1[1], start, end, "1", sep = ":")
          ))
        junct.annot <- rbind(junct.annot, junctions)
      }
      if (ggtf.gene$V7[1] == "-") {
        ggtf.gene <- ggtf.gene[order(ggtf.gene$V4), ]
        start <- ggtf.gene$V5 + 1
        start <- start[-c(length(start))]
        end <- ggtf.gene$V4 - 1
        end <- end[-c(1)]
        junctions <-
          as.data.frame(cbind(
            gene = gene.ids.u[i],
            junction = paste(ggtf.gene$V1[1], start, end, "2", sep = ":")
          ))
        junct.annot <- rbind(junct.annot, junctions)
      }
    } else{
      if (ggtf.gene$V7[1] == "+") {
        transc.pos <- which(ggtf.gene$V3 == "transcript")
        for (j in 1:length(transc.pos)) {
          if (j != length(transc.pos)) {
            transc.exons <-
              ggtf.gene[c(c(transc.pos[j] + 1):c(transc.pos[j + 1] - 1)), ]
            transc.exons <-
              transc.exons[which(transc.exons$V3 == "exon"), ]
            start <- transc.exons$V5 + 1
            start <- start[-c(length(start))]
            end <- transc.exons$V4 - 1
            end <- end[-c(1)]
            junctions <-
              as.data.frame(cbind(
                gene = gene.ids.u[i],
                junction = paste(transc.exons$V1[1], start, end, "1", sep = ":")
              ))
            junct.annot <- rbind(junct.annot, junctions)
          } else{
            transc.exons <- ggtf.gene[c(c(transc.pos[j] + 1):c(nrow(ggtf.gene))), ]
            transc.exons <-
              transc.exons[which(transc.exons$V3 == "exon"), ]
            start <- transc.exons$V5 + 1
            start <- start[-c(length(start))]
            end <- transc.exons$V4 - 1
            end <- end[-c(1)]
            junctions <-
              as.data.frame(cbind(
                gene = gene.ids.u[i],
                junction = paste(transc.exons$V1[1], start, end, "1", sep = ":")
              ))
            junct.annot <- rbind(junct.annot, junctions)

          }
        }
      }
      if (ggtf.gene$V7[1] == "-") {
        transc.pos <- which(ggtf.gene$V3 == "transcript")
        for (j in 1:length(transc.pos)) {
          if (j != length(transc.pos)) {
            transc.exons <-
              ggtf.gene[c(c(transc.pos[j] + 1):c(transc.pos[j + 1] - 1)), ]
            transc.exons <-
              transc.exons[which(transc.exons$V3 == "exon"), ]
            transc.exons <- transc.exons[order(transc.exons$V4), ]
            start <- transc.exons$V5 + 1
            start <- start[-c(length(start))]
            end <- transc.exons$V4 - 1
            end <- end[-c(1)]
            junctions <-
              as.data.frame(cbind(
                gene = gene.ids.u[i],
                junction = paste(transc.exons$V1[1], start, end, "2", sep = ":")
              ))
            junct.annot <- rbind(junct.annot, junctions)
          } else{
            transc.exons <- ggtf.gene[c(c(transc.pos[j] + 1):c(nrow(ggtf.gene))), ]
            transc.exons <-
              transc.exons[which(transc.exons$V3 == "exon"), ]
            transc.exons <- transc.exons[order(transc.exons$V4), ]
            start <- transc.exons$V5 + 1
            start <- start[-c(length(start))]
            end <- transc.exons$V4 - 1
            end <- end[-c(1)]
            junctions <-
              as.data.frame(cbind(
                gene = gene.ids.u[i],
                junction = paste(transc.exons$V1[1], start, end, "2", sep = ":")
              ))
            junct.annot <- rbind(junct.annot, junctions)
          }
        }
      }
    }
  }
  junct.annot <-
    junct.annot[which(duplicated(junct.annot$junction) == "FALSE"), ]
  junct.annot <- junct.annot[-c(grep(":::", junct.annot$junction)), ]
  return(junct.annot)
}
