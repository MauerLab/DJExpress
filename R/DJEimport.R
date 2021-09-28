#' DJEimport: Construction of junction expression matrix
#'
#' Uses sample-level junction quantification files to construct expression matrix. Returns output object for DJEannotate.
#' @param workDir path to folder where individual junction quantification files (e.g. from STAR alignment) are located. This folder should only contain junction quantification files.
#' @param data.type One of c("sample", "matrix"), indicating whether the input files are matrices per individual samples or joint into a single expression matrix.
#' @param aligner One of c("STAR", "other"), indicating the alignment tool used to produce junction quantification. Only used if data.type = "matrix".
#' If "other" is indicated, files in workDir path should contain junction IDs (with the format chr:start:end:strand) in the first column read counts in the second column.
#' @param min.expressed Numeric. Minimal number of junction reads to be considered expressed in sample. Default is 3 reads.
#'
#' @return Matrix of samples as columns and junction quantifications as rows and associated matrix with genomic coordinates.
#' @examples
#' in.file <- system.file("extdata", "junct.quant", package = "DJExpress")
#' out.file <- DJEimport(in.file, aligner="STAR")
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEimport <- function(workDir = NULL,
                      data.type = c("sample", "matrix"),
                      aligner = c("STAR", "other"),
                      min.expressed = 3) {
  if (!is.null(workDir)){
    if (is.null(workDir)){
      stop("path to junction quantification files not provided")
    }
    if (!dir.exists(workDir)) stop("folder not found")
    x <- NULL
    files = list.files(workDir)
    if(aligner=="STAR"){
      i <- 1
      junctions <-
        readr::read_tsv(paste0(workDir,"/", files[i]), comment = "", col_names = FALSE)
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
      junctions$junction_id <- as.character(junctions$junction_id)
      colnames(junctions)[2] <- files[i]
      for (i in 2:length(files)) {
        added <- readr::read_tsv(paste0(workDir,"/", files[i]), comment = "", col_names = FALSE)
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
        junctions[, i][which(junctions[, i] < min.expressed | is.na(junctions[, i]))] = 0
      }
    }else{
      i <- 1
      junctions <-
        readr::read_tsv(paste0(workDir,"/", files[i]), comment = "", col_names = FALSE)
      junctions <-
        as.data.frame(junctions)
      colnames(junctions)[1] <- "junction_id"
      junctions$junction_id <- as.character(junctions$junction_id)
      colnames(junctions)[2] <- files[i]
      for (i in 2:length(files)) {
        added <- readr::read_tsv(paste0(workDir,"/", files[i]), comment = "", col_names = FALSE)
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
        junctions[, i][which(junctions[, i] < min.expressed | is.na(junctions[, i]))] = 0
      }
    }
  }else{
    junctions <- readr::read_tsv(paste0(workDir,"/", files[i]), comment = "", col_names = FALSE)
    junctions <-
      as.data.frame(junctions)
    colnames(junctions)[1] <- "junction_id"
    if(class(junctions$junction_id)=="numeric"){
      stop("First column in matrix is numeric. Please ensure that the first column in your expression matrix correspond to junction IDs")
    }
    junctions$junction_id <- as.character(junctions$junction_id)
    rownames(junctions) <- junctions$junction_id
    junctions <- junctions[, -c(1)]
    for (i in 1:ncol(junctions)) {
      junctions[, i][which(junctions[, i] < min.expressed | is.na(junctions[, i]))] = 0
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
  return(list(quant = junctions,
              coord = junct_coord))
}
