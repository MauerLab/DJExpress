#' DJEvsTrait: Junction relationship to external trait
#'
#' Correlation between junction expression and sample traits
#' @param analize.out output object from DJEanalize()
#' @param Group1 vector or factor specifying basic control sample names
#' @param traitData a numeric vector or a matrix of external sample traits
#' @param nThreads numeric: number of threads to allow.
#' @param select.junctions a vector of junction IDs used to generate SpliceRadar plots when test.type = "Correlation". Junction ID vector length should be < 5
#' @param coeff numeric: correlation coefficient cutoff for significance
#' @param pvalue numeric: correlation p-value cutoff for significance
#' @param plot.heatmap logical, should junction vs trait correlation heatmap plot be saved as output?. Default is FALSE
#' @param test.type  a character string specifying the type of junction-trait association test to use. Options are "Correlation", "Regression" or "both".
#' @param cor.method  a character string specifying the method to be used for correlation as in WGCNA. Options are "bicor", "pearson", "kendall" or "spearman".
#' @param covariates a numeric vector or a matrix of external sample covariates to be used in regression models. Valid when test.type = "Regression" or "both".
#' @param useModel  integer specifying the model type to be used for junction-trait association as in MatrixEQTL. Options are modelLINEAR, modelANOVA or modelLINEAR_CROSS. Valid when test.type = "Regression" or "both".
#'
#' @return When "Correlation" is selected as test.type, DJEvsTrait returns the matrices with correlation coefficients and associated P-values for each junction-trait pair.
#' sig.cor contains top significant associations between junctions in select.junctions and traits for SpliceRadar plotting.
#' When selected "Regression" as test.type, DJEvsTrait returns a data frame with significant junction-trait associations based on selected model in useModel.
#' "both" test.type returns correlation and regression outputs.
#' @seealso \code{\link[WGCNA]{cor}}
#' @seealso \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#' @examples
#' data(DJEanlz)
#' SF <- system.file("extdata", "SF.expr.rds", package = "DJExpress")
#' SF.exp <- readRDS(SF)
#' Group1 <- colnames(DJEanlz$v.norm$E)[grep("SRR", colnames(DJEanlz$v.norm$E))]
#' DT.out <- DJEvsTrait(analize.out = DJEanlz, Group1 = Group1,traitData = SF.exp,
#' coeff = 0.2,select.junctions = c("chr1:225688773:225695652:2",
#' "chr1:225692756:225695652:2","chr1:225688773:225692692:2"),
#' test.type = "Correlation", cor.method = "bicor")
#' @importFrom grDevices dev.off recordPlot
#' @importFrom WGCNA bicor
#' @importFrom MatrixEQTL modelANOVA modelLINEAR modelLINEAR_CROSS
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEvsTrait <-
  function (analize.out,
            Group1,
            traitData,
            nThreads = 2,
            select.junctions = NULL,
            coeff = 0.3,
            pvalue = 0.05,
            plot.heatmap = FALSE,
            test.type = c("Correlation", "Regression", "both"),
            cor.method = c("bicor", "pearson", "kendall", "spearman"),
            useModel = c(modelANOVA, modelLINEAR, modelLINEAR_CROSS),
            covariates = NULL)
  {
    if(test.type=="Correlation" | test.type=="both"){
      options(stringsAsFactors = FALSE)
      datExpr0 = analize.out$v.norm$E
      datExpr = datExpr0[, -c(match(Group1, colnames(datExpr0)))]
      datExpr = as.data.frame(t(datExpr))
      WGCNA::enableWGCNAThreads(nThreads = nThreads)
      if (cor.method == "bicor"){
        TraitCor = WGCNA::bicor(datExpr, traitData, use = "p")
      }else{
        TraitCor = WGCNA::cor(datExpr, traitData, method = cor.method, use = "p")
      }
      TraitPvalue = WGCNA::corPvalueStudent(TraitCor, nrow(datExpr))
      if (plot.heatmap) {
        WGCNA::sizeGrWindow(10, 6)
        textMatrix = paste(signif(TraitCor, 2),
                           "\n(",
                           signif(TraitPvalue, 1),
                           ")",
                           sep = "")
        dim(textMatrix) = dim(TraitCor)
        p <- WGCNA::labeledHeatmap(
          Matrix = TraitCor,
          xLabels = colnames(traitData),
          yLabels = colnames(datExpr),
          ySymbols = colnames(datExpr),
          colorLabels = FALSE,
          colors = WGCNA::blueWhiteRed(50),
          textMatrix = textMatrix,
          setStdMargins = FALSE,
          cex.text = 0.5,
          zlim = c(-1, 1),
          main = paste("Junction-trait relationships")
        )
      }
      if(test.type=="Correlation"){
        if (!is.null(select.junctions)) {
          if (length(select.junctions) >= 5)
            stop("select.junctions in DJEvsTrait() should be a vector of length < 5")
          print(paste0("subsetting associations to: ", select.junctions))
          TraitCor.sub <-
            TraitCor[match(select.junctions, rownames(TraitCor)), , drop = FALSE]
          TraitPvalue.sub <-
            TraitPvalue[match(select.junctions, rownames(TraitPvalue)), , drop = FALSE]
          sig.cors <- list()
          for (i in 1:nrow(TraitCor.sub)) {
            sig.cor <-
              which(abs(round(TraitCor.sub[i, ], 1)) >= coeff &
                      TraitPvalue.sub[i, ] <= pvalue)
            sig.traits <-
              as.data.frame(cbind(gene = colnames(TraitCor.sub)[sig.cor], coef = TraitCor.sub[i, sig.cor]))
            sig.cors[[i]] <- sig.traits
          }
          names(sig.cors) <- rownames(TraitCor.sub)
          final.gene <- unique(do.call("rbind", sig.cors)[, 1])
            if (length(final.gene) > 0){
              if (plot.heatmap) {
                return(
                  list(
                    heatmap = p,
                    TraitCor = TraitCor,
                    TraitPvalue = TraitPvalue,
                    sig.cor = sig.cors
                  )
                )
              } else{
                return(list(
                  TraitCor = TraitCor,
                  TraitPvalue = TraitPvalue,
                  sig.cor = sig.cors
                ))
              }
            }
            if (length(final.gene) == 0)
              warning("no significant junction-trait correlations found")
            if (plot.heatmap) {
              return(
                list(
                  heatmap = p,
                  TraitCor = TraitCor,
                  TraitPvalue = TraitPvalue,
                  sig.cor = NULL
                )
              )
            } else{
              return(list(
                TraitCor = TraitCor,
                TraitPvalue = TraitPvalue,
                sig.cor = NULL
              ))
            }
        }else{
          warning("is.null(select.junctions): no junctions selected for SpliceRadar")
          if (plot.heatmap) {
            return(list(
              heatmap = p,
              TraitCor = TraitCor,
              TraitPvalue = TraitPvalue,
              sig.cor = NULL
            ))
          } else{
            return(list(
              TraitCor = TraitCor,
              TraitPvalue = TraitPvalue,
              sig.cor = NULL
            ))
          }
        }
      }
    }
    if(test.type=="Regression" | test.type=="both"){
      if(length(useModel)!=1)
        stop("Regression test.type selected, but any useModel specified")
      options(stringsAsFactors = FALSE)
      datExpr0 = analize.out$v.norm$E
      datExpr0 = datExpr0[, -c(match(Group1, colnames(datExpr0)))]
      datExpr <- as.data.frame(matrix(as.numeric(as.character(unlist(datExpr0))),nrow = nrow(datExpr0)))
      datExpr = as.data.frame(cbind(junctionID=rownames(datExpr0),
                                    datExpr))
      rownames(datExpr) <- NULL
      tData0 <- t(traitData)
      tData <- as.data.frame(matrix(as.numeric(as.character(unlist(tData0))),nrow = nrow(tData0)))
      tData = as.data.frame(cbind(traitID=rownames(tData0),
                                  tData))
      rownames(tData) <- NULL
      if(!is.null(covariates)){
        covData0 <- t(covariates)
        covData <- as.data.frame(matrix(as.numeric(as.character(unlist(covData0))),nrow = nrow(covData0)))
        covData = as.data.frame(cbind(traitID=rownames(covData0),
                                      covData))
        rownames(covData) <- NULL
        write.table(datExpr, file = "datExpr.txt", sep = "\t", row.names = FALSE)
        write.table(tData, file = "tData.txt", sep = "\t", row.names = FALSE)
        write.table(covData, file = "covData.txt", sep = "\t", row.names = FALSE)
        tData2 <- paste0(getwd(), "/", "tData.txt")
        datExpr2 <- paste0(getwd(), "/",  "datExpr.txt")
        covData2 <- paste0(getwd(), "/",  "covData.txt")
        output_file_name = tempfile()
        ## Load trait data
        snps = MatrixEQTL::SlicedData$new()
        snps$fileDelimiter = "\t"      # the TAB character
        snps$fileOmitCharacters = "NA" # denote missing values;
        snps$fileSkipRows = 1          # one row of column labels
        snps$fileSkipColumns = 1       # one column of row labels
        snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        snps$LoadFile( tData2 )
        ## Load junction expression data
        gene = MatrixEQTL::SlicedData$new()
        gene$fileDelimiter = "\t"      # the TAB character
        gene$fileOmitCharacters = "NA" # denote missing values;
        gene$fileSkipRows = 1          # one row of column labels
        gene$fileSkipColumns = 1       # one column of row labels
        gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        gene$LoadFile( datExpr2 )
        ## Load covariates
        cvrt = MatrixEQTL::SlicedData$new()
        cvrt$fileDelimiter = "\t"      # the TAB character
        cvrt$fileOmitCharacters = "NA" # denote missing values;
        cvrt$fileSkipRows = 1          # one row of column labels
        cvrt$fileSkipColumns = 1       # one column of row labels
        cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        cvrt$LoadFile( covData2 )
        errorCovariance = numeric()
        me = MatrixEQTL::Matrix_eQTL_engine(
          snps = snps,
          gene = gene,
          output_file_name = output_file_name,
          pvOutputThreshold = pvalue,
          useModel = useModel,
          errorCovariance = errorCovariance,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = FALSE,
          noFDRsaveMemory = FALSE)
        colnames(me$all$eqtls) <- c("trait", "junction", "statistic", "pvalue", "FDR", "beta")
        if(nrow(me$all$eqtls)==0)
          warning("no significant junction-trait associations by regression found")
      }else{
        write.table(datExpr, file = "datExpr.txt", sep = "\t", row.names = FALSE)
        write.table(tData, file = "tData.txt", sep = "\t", row.names = FALSE)
        tData2 <- paste0(getwd(), "/", "tData.txt")
        datExpr2 <- paste0(getwd(), "/",  "datExpr.txt")
        output_file_name = tempfile()
        ## Load trait data
        snps = MatrixEQTL::SlicedData$new()
        snps$fileDelimiter = "\t"      # the TAB character
        snps$fileOmitCharacters = "NA" # denote missing values;
        snps$fileSkipRows = 1          # one row of column labels
        snps$fileSkipColumns = 1       # one column of row labels
        snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        snps$LoadFile( tData2 )
        ## Load junction expression data
        gene = MatrixEQTL::SlicedData$new()
        gene$fileDelimiter = "\t"      # the TAB character
        gene$fileOmitCharacters = "NA" # denote missing values;
        gene$fileSkipRows = 1          # one row of column labels
        gene$fileSkipColumns = 1       # one column of row labels
        gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
        gene$LoadFile( datExpr2 )
        errorCovariance = numeric()
        me = MatrixEQTL::Matrix_eQTL_engine(
          snps = snps,
          gene = gene,
          output_file_name = output_file_name,
          pvOutputThreshold = pvalue,
          useModel = useModel,
          errorCovariance = errorCovariance,
          verbose = TRUE,
          pvalue.hist = TRUE,
          min.pv.by.genesnp = FALSE,
          noFDRsaveMemory = FALSE)
        colnames(me$all$eqtls) <- c("trait", "junction", "statistic", "pvalue", "FDR", "beta")
        if(nrow(me$all$eqtls)==0)
          warning("no significant junction-trait associations by regression found")
      }
      if(test.type=="Regression"){
        return(list(
          Sig.reg.out = me$all$eqtls
        ))
      }
      if(test.type=="both"){
        if (!is.null(select.junctions)) {
          if (length(select.junctions) >= 5)
            stop("select.junctions in DJEvsTrait() should be a vector of length < 5")
          print(paste0("subsetting associations to: ", select.junctions))
          TraitCor.sub <-
            TraitCor[match(select.junctions, rownames(TraitCor)), , drop = FALSE]
          TraitPvalue.sub <-
            TraitPvalue[match(select.junctions, rownames(TraitPvalue)), , drop = FALSE]
          sig.cors <- list()
          for (i in 1:nrow(TraitCor.sub)) {
            sig.cor <-
              which(abs(round(TraitCor.sub[i, ], 1)) >= coeff &
                      TraitPvalue.sub[i, ] <= pvalue)
            sig.traits <-
              as.data.frame(cbind(gene = colnames(TraitCor.sub)[sig.cor], coef = TraitCor.sub[i, sig.cor]))
            sig.cors[[i]] <- sig.traits
          }
          names(sig.cors) <- rownames(TraitCor.sub)
          final.gene <- unique(do.call("rbind", sig.cors)[, 1])
          if (length(final.gene) > 0){
            if (plot.heatmap) {
              return(
                list(
                  heatmap = p,
                  TraitCor = TraitCor,
                  TraitPvalue = TraitPvalue,
                  sig.cor = sig.cors,
                  Sig.reg.out = me$all$eqtls
                )
              )
            } else{
              return(list(
                TraitCor = TraitCor,
                TraitPvalue = TraitPvalue,
                sig.cor = sig.cors,
                Sig.reg.out = me$all$eqtls
              ))
            }
          }
          if (length(final.gene) == 0)
            warning("no significant junction-trait correlations found")
          if (plot.heatmap) {
            return(
              list(
                heatmap = p,
                TraitCor = TraitCor,
                TraitPvalue = TraitPvalue,
                sig.cor = NULL,
                Sig.reg.out = me$all$eqtls
              )
            )
          } else{
            return(list(
              TraitCor = TraitCor,
              TraitPvalue = TraitPvalue,
              sig.cor = NULL,
              Sig.reg.out = me$all$eqtls
            ))
          }
        }else{
          warning("is.null(select.junctions): no junctions selected for SpliceRadar")
          if (plot.heatmap) {
            return(list(
              heatmap = p,
              TraitCor = TraitCor,
              TraitPvalue = TraitPvalue,
              sig.cor = NULL,
              Sig.reg.out = me$all$eqtls
            ))
          } else{
            return(list(
              TraitCor = TraitCor,
              TraitPvalue = TraitPvalue,
              sig.cor = NULL,
              Sig.reg.out = me$all$eqtls
            ))
          }
        }
      }
    }
  }
