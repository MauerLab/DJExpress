#' DJEplotSplice: Differential junction usage plot
#'
#' Plots relative log-fold changes by junctions for the specified gene and highlights the significantly spliced junctions. This is a wrapper function of limma plotSplice.
#' @param analize.out output object from DJEanalize()
#' @param geneID gene name for which plot is generated
#' @param gtf Reference transcriptome in genecode gtf file format. Used to define gene coordinates
#' @param target.junction Junction IDs to highlight within gene model plot. Maximum 5 junctions.
#' @param genecolname column name of fit$genes containing gene IDs. Defaults to fit$genecolname
#' @param rank integer, if geneid=NULL then this ranked gene will be plotted
#' @param logFC numeric: log2-fold-change value cutoff in junction usage to define significance in plot
#' @param meanExp numeric: mean expression value cutoff in tested condition group to define significance in plot
#' @param FDR numeric: adjusted p-value cutoff in the linear model fit to define significance
#'
#' @return interactive highcharts plot showing relative logFC across junctions in target gene
#' @examples
#' data(DJEanlz)
#' iPlot.out <- DJEplotSplice(DJEanlz, geneID="ENAH", logFC = 0.8, FDR = 0.05)
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEplotSplice <-
  function (analize.out,
            geneID = NULL,
            gtf = NULL,
            target.junction = NULL,
            genecolname = NULL,
            rank = 1L,
            logFC = 2,
            meanExp = 10,
            FDR = 0.05)
  {
    junctionID <- NULL
    fit = analize.out$ex.norm
    coef = ncol(fit)
    diffSout = analize.out$dje.out
    SigTopTable = analize.out$dje.sig
    geneid = analize.out$dje.out$GeneID[match(geneID, analize.out$dje.out$GeneID)]
    if (is.null(genecolname))
      genecolname <- fit$genecolname
    else
      genecolname <- as.character(genecolname)
    if (is.null(geneid)) {
      if (rank == 1L)
        i <- which.min(fit$gene.F.p.value[, coef])
      else
        i <- order(fit$gene.F.p.value[, coef])[rank]
      geneid <- paste(fit$gene.genes[i, genecolname], collapse = ".")
    }
    else {
      geneid <- as.character(geneid)
      i <- which(fit$gene.genes[, genecolname] == geneid)[1]
      if (!length(i))
        stop(paste("geneid", geneid, "not found"))
    }
    j <- fit$gene.firstjunction[i]:fit$gene.lastjunction[i]
    junctioncolname <- fit$junctioncolname
    strcol <-
      grepl("strand", colnames(fit$gene.genes), ignore.case = TRUE)
    if (any(strcol))
      geneid <- paste0(geneid, " (", as.character(fit$gene.genes[i,
                                                                 strcol])[1], ")")
    if (is.null(junctioncolname)) {
      plot(
        fit$coefficients[j, coef],
        xlab = "",
        ylab = "logFC (this junction vs average-rest)",
        main = geneid,
        type = "b"
      )
    }
    else {
      junction.id <- fit$genes[j, junctioncolname]
      tjunct <- diffSout[match(junction.id, diffSout$junctionID), ]
      total.junction.ids <- tjunct$junctionID
      tjunct$FDR[which(is.na(tjunct$FDR))] = 1
      tjunct$logFC[which(is.na(tjunct$logFC))] = 0
      tjunct$meanExp.group2[which(is.na(tjunct$meanExp.group2))] = 0
      tjunct$meanExp.group2 <-
        as.numeric(as.character(tjunct$meanExp.group2))
    }
    tjunct$FDR[which(tjunct$FDR == 0)] = 1.0e-322
    fdr <- tjunct$FDR
    sig <- fdr < FDR
    tumor.withcounts.mean <-
      paste(tjunct$group2WithCounts,
            "(",
            round(tjunct$meanExp.group2, digits = 2) ,
            ")")
    mean.sig <- tjunct$meanExp.group2 >= meanExp
    NOmean.sig <- tjunct$meanExp.group2 < meanExp
    logFCs <- tjunct$logFC
    pos.logFC <- logFCs > logFC
    neg.logFC <- logFCs < -logFC
    sig.and.NOmean <- sig & NOmean.sig
    sig.and.mean <- sig & mean.sig
    sig.and.neg <- sig & neg.logFC
    sig.and.pos <- sig & pos.logFC
    sigmean.and.pos <- sig.and.mean & pos.logFC
    sigmean.and.neg <- sig.and.neg
    logFCs.ebayes <- tjunct$logFC.ebayes
    sig.logFC.e <- logFCs.ebayes < -logFC | logFCs.ebayes > logFC
    color = rep("#000000", length(total.junction.ids))
    color[which(sigmean.and.neg == "TRUE")] = "#005BA2"
    color[which(sigmean.and.pos == "TRUE")] = "#D44F4F"
    color[which(sigmean.and.pos == "TRUE" &
                  sig.logFC.e == "FALSE")] = "#D3D3D3"
    color[which(sigmean.and.pos == "TRUE" &
                  sig.logFC.e == "FALSE")] = "#D3D3D3"
    color[which(sigmean.and.neg == "TRUE" &
                  sig.logFC.e == "FALSE")] = "#D3D3D3"
    color[which(sigmean.and.neg == "TRUE" &
                  sig.logFC.e == "FALSE")] = "#D3D3D3"
    tjunct$color = color
    sigmean.and.pos <- sigmean.and.pos & sig.logFC.e
    sigmean.and.neg <- sigmean.and.neg & sig.logFC.e
    cex = rep(1, length(total.junction.ids))
    if (any(sigmean.and.pos)) {
      cex[which(sigmean.and.pos == "TRUE")] <- 3
    }
    if (any(sigmean.and.neg)) {
      cex[which(sigmean.and.neg == "TRUE")] <- 3
    }
    tjunct$cex <- cex
    sig.junct <- tjunct$junctionID[which(tjunct$color != "#000000")]
    tjunct$GrouplogFC <- NA
    tjunct$GrouplogFC[which(tjunct$color != "#000000")] <-
      SigTopTable$group[match(sig.junct, SigTopTable$junctionID)]
    tjunct$GrouplogFC[which(is.na(tjunct$GrouplogFC))] <-
      "no significant"
    tjunct$group2WithCounts <-
      round(as.numeric(as.character(tjunct$group2WithCounts)), digits = 2)
    tjunct$meanExp.group2 <-
      round(as.numeric(as.character(tjunct$meanExp.group2)), digits = 2)
    tjunct$logFC.ebayes <-
      round(as.numeric(as.character(tjunct$logFC.ebayes)), digits = 2)
    colnames(tjunct)[which(colnames(tjunct) == "meanExp.group2")] <-
      "meanExp_group2"
    colnames(tjunct)[which(colnames(tjunct) == "logFC.ebayes")] <-
      "logFC_absolute"
    p <- highcharter::highchart() %>%
      highcharter::hc_add_series(
        tjunct$logFC,
        type = "line",
        showInLegend = F,
        color = "grey"
      ) %>%
      highcharter::hc_add_series(
        data = tjunct,
        highcharter::hcaes(
          x = junctionID,
          y = logFC,
          color = color,
          size = cex
        ),
        type = "scatter",
        showInLegend = F,
        maxSize = "5%"
      ) %>%
      highcharter::hc_xAxis(title = list(text = paste("junctions in", geneid)), categories = total.junction.ids) %>%
      highcharter::hc_yAxis_multiples(
        list(
          title = list(text = "logFC (diffSplice)"),
          lineWidth = 3,
          minorGridLineWidth = 0,
          gridLineWidth = 0,
          plotBands = list(
            list(
              from = -logFC,
              to = 0,
              color = "lightgrey"
            ),
            list(from = 0, to = logFC, color = "lightgrey")
          ),
          plotLines = list(list(
            value = 0,
            color = "black",
            width = 3
          ))
        )
      ) %>%
      highcharter::hc_tooltip(pointFormat = "junction: {point.junctionID} <br> logFC diffSplice: {point.y} <br> logFC absolute: {point.logFC_absolute} <br> FDR: {point.FDR} <br> group: {point.GrouplogFC} <br> mean raw counts: {point.meanExp_group2} <br> number of samples with counts: {point.group2WithCounts} ",
                              crosshairs = TRUE)
    tjunct <- tjunct[,-c(match("color", colnames(tjunct)))]
    rownames(tjunct) <- NULL
    if(is.null(gtf)|is.null(target.junction)){
      warning("Path to gtf file not provided and/or target junction not selected. Junction-to-gene model plot won't be displayed")
      return(
        list(
          plot = p,
          junction.info = tjunct
        )
      )
    }else{
      if (unlist(base::strsplit(gtf, "[.]"))[length(unlist(base::strsplit(gtf, "[.]")))] ==
          "gz") {
        ggtf <-
          utils::read.delim(gzfile(gtf), header = FALSE, comment.char = "#")
      } else{
        ggtf <- utils::read.delim(gtf, header = FALSE, comment.char = "#")
      }
      domtogen <- system.file("extdata", "domain.to.genome.csv.gz", package = "DJExpress")
      domain.to.genome <-
        utils::read.delim(gzfile(domtogen), header = TRUE, sep=";")
      gID <- tjunct$GeneID[1]
      if(grepl("gene_id ENSG", as.character(ggtf$V9[1]))=="TRUE"){
        ggtf.gID <- ggtf[grep(gID, ggtf$V9),]
        target <- stringr::str_match(as.character(ggtf.gID$V9), "gene_name \\s*(.*?)\\s*; ")
        ggtf.gID <- ggtf.gID[which(target[,2] ==gID),]
        ggtf.gID.cds <- ggtf.gID[c(grep("CDS", ggtf.gID$V3)),]
        ggtf.gID <- ggtf.gID[c(grep("exon", ggtf.gID$V3)),]
        ggtf.gID.gene <- as.data.frame(cbind(ggtf.gID[,c(4,5)], y1=2, y2=3))
        window.range <- max(ggtf.gID.gene$V4)-min(ggtf.gID.gene$V4)
        arrow.dir <- list()
        targetjunct <- list()
        junct_col <- list()
        for (i in 1:length(target.junction)){
          junction <- target.junction[i]
          targetjunct[[i]] <- data.frame(x1 = as.numeric(as.character(unlist(strsplit(as.character(junction), "[:]"))[2])),
                                         x2 = as.numeric(as.character(unlist(strsplit(as.character(junction), "[:]"))[3])),
                                         y1 = 3, y2 = 3)
          if(as.character(unlist(strsplit(as.character(junction), "[:]"))[4]) == "1"){
            arrow.dir[[i]] <- c(min(ggtf.gID.gene$V4), min(ggtf.gID.gene$V4)+window.range/10)
          }
          if(as.character(unlist(strsplit(as.character(junction), "[:]"))[4]) == "2"){
            arrow.dir[[i]] <- c(max(ggtf.gID.gene$V4), max(ggtf.gID.gene$V4)-window.range/10)
          }
          if(as.numeric(as.character(tjunct$logFC[match(target.junction[i], tjunct$junctionID)])) < 0){
            junct_col[[i]] <- "#005BA2"
          }
          if(as.numeric(as.character(tjunct$logFC[match(target.junction[i], tjunct$junctionID)])) > 0){
            junct_col[[i]] <- "#D44F4F"
          }
        }
        ggtf.gID.gene.2 <- c()
        for(i in 1:nrow(ggtf.gID.gene)){
          if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE &
             dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE){

            ggtf.gID.gene.b <- as.data.frame(cbind(V4 = ggtf.gID.gene$V4[i],
                                                   V5 = ggtf.gID.gene$V5[i],
                                                   y1 = 2,
                                                   y2 = 3))
          }
          if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE &
             dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE){
            ggtf.gID.gene.b <- as.data.frame(cbind(V4 = ggtf.gID.gene$V4[i],
                                                   V5 = ggtf.gID.gene$V5[i],
                                                   y1 = 2.2,
                                                   y2 = 2.8))
          }
          if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE &
             dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE){
            ggtf.gID.gene.b1 <- as.data.frame(cbind(V4= ggtf.gID.gene$V4[i],
                                                    V5= min(ggtf.gID.cds$V4),
                                                    y1 = 2.2,
                                                    y2 = 2.8))
            ggtf.gID.gene.b2 <- as.data.frame(cbind(V4= min(ggtf.gID.cds$V4),
                                                    V5= ggtf.gID.gene$V5[i],
                                                    y1 = 2,
                                                    y2 = 3))
            ggtf.gID.gene.b <- as.data.frame(rbind(ggtf.gID.gene.b1, ggtf.gID.gene.b2))
          }
          if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE &
             dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE){
            ggtf.gID.gene.b1 <- as.data.frame(cbind(V4= ggtf.gID.gene$V4[i],
                                                    V5= max(ggtf.gID.cds$V5),
                                                    y1 = 2,
                                                    y2 = 3))
            ggtf.gID.gene.b2 <- as.data.frame(cbind(V4= max(ggtf.gID.cds$V5),
                                                    V5= ggtf.gID.gene$V5[i],
                                                    y1 = 2.2,
                                                    y2 = 2.8))

            ggtf.gID.gene.b <- as.data.frame(rbind(ggtf.gID.gene.b1, ggtf.gID.gene.b2))
          }
          ggtf.gID.gene.2 <- rbind(ggtf.gID.gene.2, ggtf.gID.gene.b)
        }
        domain.to.genome.geneID <- domain.to.genome[which(domain.to.genome$gene==geneID),c(14,15,7)]
        domain.to.genome.geneID <- as.data.frame(cbind(domain.to.genome.geneID, y1=2, y2=3))
        domain.to.genome.geneID$chr_start <- as.numeric(as.character(domain.to.genome.geneID$chr_start))
        domain.to.genome.geneID$chr_end <- as.numeric(as.character(domain.to.genome.geneID$chr_end))
        domain.to.genome.geneID$y1 <- as.numeric(as.character(domain.to.genome.geneID$y1))
        domain.to.genome.geneID$y2 <- as.numeric(as.character(domain.to.genome.geneID$y2))
        colnames(domain.to.genome.geneID)[3] <- "Domain"

        nb.cols <- length(levels(factor(domain.to.genome.geneID$Domain)))
        mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nb.cols)
        if(length(target.junction) ==1){
          JGplot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste(gID, "gene model")) +
            ggplot2::scale_x_continuous(name="")+
            ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
            ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
            ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
            ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
            ggplot2::scale_fill_manual(values = mycolors)+
            ggplot2::labs(fill = "Domains and PTMs")+
            ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::theme_classic()+
            ggplot2::theme(
              title = ggplot2::element_text(size=14),
              plot.title = ggplot2::element_text(hjust = 0.5),
              axis.text.y = ggplot2::element_blank(),
              line = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
              axis.ticks = ggplot2::element_blank())
        }
        if(length(target.junction) ==2){
          JGplot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste(gID, "gene model")) +
            ggplot2::scale_x_continuous(name="")+
            ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
            ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
            ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
            ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
            ggplot2::scale_fill_manual(values = mycolors)+
            ggplot2::labs(fill = "Domains and PTMs")+
            ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::theme_classic()+
            ggplot2::theme(
              title = ggplot2::element_text(size=14),
              plot.title = ggplot2::element_text(hjust = 0.5),
              axis.text.y = ggplot2::element_blank(),
              line = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
              axis.ticks = ggplot2::element_blank())
        }
        if(length(target.junction) ==3){
          JGplot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste(gID, "gene model")) +
            ggplot2::scale_x_continuous(name="")+
            ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
            ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
            ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
            ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
            ggplot2::scale_fill_manual(values = mycolors)+
            ggplot2::labs(fill = "Domains and PTMs")+
            ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::theme_classic()+
            ggplot2::theme(
              title = ggplot2::element_text(size=14),
              plot.title = ggplot2::element_text(hjust = 0.5),
              axis.text.y = ggplot2::element_blank(),
              line = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
              axis.ticks = ggplot2::element_blank())
        }
        if(length(target.junction) ==4){
          JGplot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste(gID, "gene model")) +
            ggplot2::scale_x_continuous(name="")+
            ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
            ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
            ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
            ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
            ggplot2::scale_fill_manual(values = mycolors)+
            ggplot2::labs(fill = "Domains and PTMs")+
            ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[4]], color=junct_col[[4]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::theme_classic()+
            ggplot2::theme(
              title = ggplot2::element_text(size=14),
              plot.title = ggplot2::element_text(hjust = 0.5),
              axis.text.y = ggplot2::element_blank(),
              line = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
              axis.ticks = ggplot2::element_blank())
        }
        if(length(target.junction) ==5){
          JGplot <- ggplot2::ggplot() +
            ggplot2::ggtitle(paste(gID, "gene model")) +
            ggplot2::scale_x_continuous(name="")+
            ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
            ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
            ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
            ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
            ggplot2::scale_fill_manual(values = mycolors)+
            ggplot2::labs(fill = "Domains and PTMs")+
            ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                  arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[4]], color=junct_col[[4]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                data = targetjunct[[5]], color=junct_col[[5]], curvature = -1, ncp = 10, linetype=2)+
            ggplot2::theme_classic()+
            ggplot2::theme(
              title = ggplot2::element_text(size=14),
              plot.title = ggplot2::element_text(hjust = 0.5),
              axis.text.y = ggplot2::element_blank(),
              line = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
              axis.ticks = ggplot2::element_blank())
        }
      }else{
        if(grepl("gene_id ENSG", as.character(ggtf$V9[1]))=="FALSE"){
          ggtf.gID <- ggtf[grep(gID, ggtf$V9),]
          target <- stringr::str_match(as.character(ggtf.gID$V9), "gene_id \\s*(.*?)\\s*; ")
          ggtf.gID <- ggtf.gID[which(target[,2] ==gID),]
          ggtf.gID.cds <- ggtf.gID[c(grep("CDS", ggtf.gID$V3)),]
          ggtf.gID <- ggtf.gID[c(grep("exon", ggtf.gID$V3)),]
          ggtf.gID.gene <- as.data.frame(cbind(ggtf.gID[,c(4,5)], y1=2, y2=3))
          window.range <- max(ggtf.gID.gene$V4)-min(ggtf.gID.gene$V4)
          arrow.dir <- c()
          if(as.character(unlist(strsplit(as.character(target.junction[1]), "[:]"))[4]) == "1"){
            arrow.dir <- c(min(ggtf.gID.gene$V4), min(ggtf.gID.gene$V4)+window.range/10)
          }
          if(as.character(unlist(strsplit(as.character(target.junction[1]), "[:]"))[4]) == "2"){
            arrow.dir <- c(max(ggtf.gID.gene$V4), max(ggtf.gID.gene$V4)-window.range/10)
          }
          targetjunct <- list()
          junct_col <- list()
          for (i in 1:length(target.junction)){
            junction <- target.junction[i]
            targetjunct[[i]] <- data.frame(x1 = as.numeric(as.character(unlist(strsplit(as.character(junction), "[:]"))[2])),
                                           x2 = as.numeric(as.character(unlist(strsplit(as.character(junction), "[:]"))[3])),
                                           y1 = 3, y2 = 3)
            if(as.numeric(as.character(tjunct$logFC[match(target.junction[i], tjunct$junctionID)])) < 0){
              junct_col[[i]] <- "#005BA2"
            }
            if(as.numeric(as.character(tjunct$logFC[match(target.junction[i], tjunct$junctionID)])) > 0){
              junct_col[[i]] <- "#D44F4F"
            }
          }
          ggtf.gID.gene.2 <- c()
          for(i in 1:nrow(ggtf.gID.gene)){
            if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE &
               dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE){

              ggtf.gID.gene.b <- as.data.frame(cbind(V4 = ggtf.gID.gene$V4[i],
                                                     V5 = ggtf.gID.gene$V5[i],
                                                     y1 = 2,
                                                     y2 = 3))
            }
            if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE &
               dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE){
              ggtf.gID.gene.b <- as.data.frame(cbind(V4 = ggtf.gID.gene$V4[i],
                                                     V5 = ggtf.gID.gene$V5[i],
                                                     y1 = 2.2,
                                                     y2 = 2.8))
            }
            if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE &
               dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE){
              ggtf.gID.gene.b1 <- as.data.frame(cbind(V4= ggtf.gID.gene$V4[i],
                                                      V5= min(ggtf.gID.cds$V4),
                                                      y1 = 2.2,
                                                      y2 = 2.8))
              ggtf.gID.gene.b2 <- as.data.frame(cbind(V4= min(ggtf.gID.cds$V4),
                                                      V5= ggtf.gID.gene$V5[i],
                                                      y1 = 2,
                                                      y2 = 3))
              ggtf.gID.gene.b <- as.data.frame(rbind(ggtf.gID.gene.b1, ggtf.gID.gene.b2))
            }
            if(dplyr::between(ggtf.gID.gene$V4[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==TRUE &
               dplyr::between(ggtf.gID.gene$V5[i], min(ggtf.gID.cds$V4), max(ggtf.gID.cds$V5))==FALSE){
              ggtf.gID.gene.b1 <- as.data.frame(cbind(V4= ggtf.gID.gene$V4[i],
                                                      V5= max(ggtf.gID.cds$V5),
                                                      y1 = 2,
                                                      y2 = 3))
              ggtf.gID.gene.b2 <- as.data.frame(cbind(V4= max(ggtf.gID.cds$V5),
                                                      V5= ggtf.gID.gene$V5[i],
                                                      y1 = 2.2,
                                                      y2 = 2.8))

              ggtf.gID.gene.b <- as.data.frame(rbind(ggtf.gID.gene.b1, ggtf.gID.gene.b2))
            }
            ggtf.gID.gene.2 <- rbind(ggtf.gID.gene.2, ggtf.gID.gene.b)
          }
          domain.to.genome.geneID <- domain.to.genome[which(domain.to.genome$gene==geneID),c(14,15,7)]
          domain.to.genome.geneID <- as.data.frame(cbind(domain.to.genome.geneID, y1=2, y2=3))
          domain.to.genome.geneID$chr_start <- as.numeric(as.character(domain.to.genome.geneID$chr_start))
          domain.to.genome.geneID$chr_end <- as.numeric(as.character(domain.to.genome.geneID$chr_end))
          domain.to.genome.geneID$y1 <- as.numeric(as.character(domain.to.genome.geneID$y1))
          domain.to.genome.geneID$y2 <- as.numeric(as.character(domain.to.genome.geneID$y2))
          colnames(domain.to.genome.geneID)[3] <- "Domain"

          nb.cols <- length(levels(factor(domain.to.genome.geneID$Domain)))
          mycolors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(nb.cols)
          if(length(target.junction) ==1){
            JGplot <- ggplot2::ggplot() +
              ggplot2::ggtitle(paste(gID, "gene model")) +
              ggplot2::scale_x_continuous(name="")+
              ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
              ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
              ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
              ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
              ggplot2::scale_fill_manual(values = mycolors)+
              ggplot2::labs(fill = "Domains and PTMs")+
              ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                    arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::theme_classic()+
              ggplot2::theme(
                title = ggplot2::element_text(size=14),
                plot.title = ggplot2::element_text(hjust = 0.5),
                axis.text.y = ggplot2::element_blank(),
                line = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
                axis.ticks = ggplot2::element_blank())
          }
          if(length(target.junction) ==2){
            JGplot <- ggplot2::ggplot() +
              ggplot2::ggtitle(paste(gID, "gene model")) +
              ggplot2::scale_x_continuous(name="")+
              ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
              ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
              ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
              ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
              ggplot2::scale_fill_manual(values = mycolors)+
              ggplot2::labs(fill = "Domains and PTMs")+
              ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                    arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::theme_classic()+
              ggplot2::theme(
                title = ggplot2::element_text(size=14),
                plot.title = ggplot2::element_text(hjust = 0.5),
                axis.text.y = ggplot2::element_blank(),
                line = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
                axis.ticks = ggplot2::element_blank())
          }
          if(length(target.junction) ==3){
            JGplot <- ggplot2::ggplot() +
              ggplot2::ggtitle(paste(gID, "gene model")) +
              ggplot2::scale_x_continuous(name="")+
              ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
              ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
              ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
              ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
              ggplot2::scale_fill_manual(values = mycolors)+
              ggplot2::labs(fill = "Domains and PTMs")+
              ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                    arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::theme_classic()+
              ggplot2::theme(
                title = ggplot2::element_text(size=14),
                plot.title = ggplot2::element_text(hjust = 0.5),
                axis.text.y = ggplot2::element_blank(),
                line = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
                axis.ticks = ggplot2::element_blank())
          }
          if(length(target.junction) ==4){
            JGplot <- ggplot2::ggplot() +
              ggplot2::ggtitle(paste(gID, "gene model")) +
              ggplot2::scale_x_continuous(name="")+
              ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
              ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
              ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
              ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
              ggplot2::scale_fill_manual(values = mycolors)+
              ggplot2::labs(fill = "Domains and PTMs")+
              ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                    arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[4]], color=junct_col[[4]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::theme_classic()+
              ggplot2::theme(
                title = ggplot2::element_text(size=14),
                plot.title = ggplot2::element_text(hjust = 0.5),
                axis.text.y = ggplot2::element_blank(),
                line = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
                axis.ticks = ggplot2::element_blank())
          }
          if(length(target.junction) ==5){
            JGplot <- ggplot2::ggplot() +
              ggplot2::ggtitle(paste(gID, "gene model")) +
              ggplot2::scale_x_continuous(name="")+
              ggplot2::scale_y_continuous(name=as.character(unlist(strsplit(as.character(junction), "[:]"))[1]),limits=c(1.5, 3.3))+
              ggplot2::geom_segment(ggplot2::aes(x = min(ggtf.gID.gene.2$V4), y = 2.5, xend = max(ggtf.gID.gene.2$V4), yend = 2.5)) +
              ggplot2::geom_rect(data=ggtf.gID.gene.2, mapping=ggplot2::aes(xmin=V4, xmax=V5, ymin=y1, ymax=y2))+
              ggplot2::geom_rect(data=domain.to.genome.geneID, mapping=ggplot2::aes(xmin=chr_start, xmax=chr_end, ymin=y1, ymax=y2, fill=Domain))+
              ggplot2::scale_fill_manual(values = mycolors)+
              ggplot2::labs(fill = "Domains and PTMs")+
              ggplot2::geom_segment(ggplot2::aes(x = arrow.dir[1], y = 1.9, xend = arrow.dir[2], yend = 1.9),
                                    arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"),type = "closed"), color="black", alpha=1) +
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[1]], color=junct_col[[1]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[2]], color=junct_col[[2]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[3]], color=junct_col[[3]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[4]], color=junct_col[[4]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::geom_curve(ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
                                  data = targetjunct[[5]], color=junct_col[[5]], curvature = -1, ncp = 10, linetype=2)+
              ggplot2::theme_classic()+
              ggplot2::theme(
                title = ggplot2::element_text(size=14),
                plot.title = ggplot2::element_text(hjust = 0.5),
                axis.text.y = ggplot2::element_blank(),
                line = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_text(angle = 0, vjust = -0.032, hjust = 0),
                axis.ticks = ggplot2::element_blank())
          }
        }
      }
      return(
        list(
          plot = p,
          JunctionToGene= JGplot,
          junction.info = tjunct
        )
      )

    }
  }
