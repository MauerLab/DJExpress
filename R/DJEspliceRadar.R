#' DJEspliceRadar: Radar plot of junction expression vs trait correlation
#'
#' Interactive plot with top correlated junction expression vs sample traits obtained with DJEvsTrait function.
#' @param DJEvsTrait.out output object from DJEvsTrait()
#' @param base.junction vector of junction IDs to be used for trait association
#' @param topTraits numeric: number of top junction-associated traits to be shown in the plot
#' @param select.traits vector of trait names to be shown in the plot
#' @param ordered.junction junction ID for which correlation coefficient values are sorted in descending order in the plot
#' @return Interactive highcharts-based radar plot showing top associations between selected junction expression and traits based on correlation coefficients
#' @examples
#' DvT <- system.file("extdata", "DvT.rds", package = "DJExpress")
#' DvT.out <- readRDS(DvT)
#' Sr.out <- DJEspliceRadar(DvT.out, ordered.junction = "chr3:39450216:39452244:1")
#' @import magrittr
#' @importFrom grDevices dev.off recordPlot
#' @importFrom graphics abline par text
#' @importFrom methods new
#' @importFrom stats cor.test dist lm median p.adjust pf predict pt var
#' @export
DJEspliceRadar <-
  function (DJEvsTrait.out,
            base.junction = NULL,
            topTraits = NULL,
            select.traits = NULL,
            ordered.junction=NULL)
  {
    if (is.null(DJEvsTrait.out$sig.cor))
      stop("is.null(sig.cor): select.junctions set to NULL in DJEvsTrait.out()")
    sig.cors <- DJEvsTrait.out$sig.cor
    TraitCor <- DJEvsTrait.out$TraitCor
    final.gene <- unique(do.call("rbind", sig.cors)[, 1])
    final.TraitCor <- TraitCor[match(names(sig.cors), rownames(TraitCor)), match(final.gene, colnames(TraitCor))]
    sig.cors.2 <-
      data.frame(junction = rownames(final.TraitCor), final.TraitCor)
    table_radar <- tibble::as_tibble(sig.cors.2)
    if (nrow(table_radar) >= 5)
      stop("select.junctions in DJEvsTrait.out() should be a vector of length < 5")
    df.2 <-
      matrix(as.numeric(as.character(unlist(t(
        table_radar[, -c(1)]
      )))), nrow = nrow(t(table_radar[, -c(1)])))
    df.2 <- data.frame(colnames(table_radar)[-c(1)] , df.2)
    junct.list <-
      paste("junction", seq(1, length(table_radar$junction)), sep = "_")
    colnames(df.2) <- c("trait", junct.list)
    junct.names <- c("trait", table_radar$junction)

    if (is.null(select.traits)) {
      if (!is.null(base.junction)) {
        base.j <- which(junct.names == base.junction)
        df.2 <- df.2[order(df.2[, base.j]), ]
        if (!is.null(topTraits)) {
          top10.base <-
            df.2[order(abs(df.2[, base.j]), decreasing = T), base.j][1:topTraits]
          df.2 <- df.2[match(top10.base, df.2[, base.j]), ]
        }

      }
      if (!is.null(ordered.junction)){
        df.2 <- df.2[order(df.2[,match(ordered.junction, table_radar$junction)]),]
      }
      if (ncol(df.2) == 2) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            )
          )
      }
      if (ncol(df.2) == 3) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
          )
      }
      if (ncol(df.2) == 4) {
        p <- highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
            ,
            list(
              name = paste0("correlation to junction ", junct.names[4]),
              data = df.2$junction_3,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            )
          )
      }
      if (ncol(df.2) == 5) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
            ,
            list(
              name = paste0("correlation to junction ", junct.names[4]),
              data = df.2$junction_3,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[5]),
              data = df.2$junction_4,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            )
          )
      }
    } else{
      if (!is.null(base.junction))
        stop(
          "!is.null(base.junction): select.traits and base.junction cannot be simultaneously assigned "
        )
      base.j <- which(junct.names == base.junction)
      df.2 <- df.2[order(df.2[, base.j]), ]
      if (!is.null(topTraits))
        stop(
          "!is.null(base.junction): select.traits and topTraits cannot be simultaneously assigned "
        )
      if (!is.null(ordered.junction)){
        df.2 <- df.2[order(df.2[,match(ordered.junction, table_radar$junction)]),]
      }
      if (ncol(df.2) == 2) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            )
          )
      }
      if (ncol(df.2) == 3) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
          )
      }
      if (ncol(df.2) == 4) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
            ,
            list(
              name = paste0("correlation to junction ", junct.names[4]),
              data = df.2$junction_3,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            )
          )
      }
      if (ncol(df.2) == 5) {
        p <-   highcharter::highchart() %>%
          highcharter::hc_chart(polar = TRUE) %>%
          highcharter::hc_title(text = "correlation junction expression vs trait") %>%
          highcharter::hc_xAxis(
            categories = df.2$trait,
            tickmarkPlacement = "on",
            lineWidth = 0
          ) %>%
          highcharter::hc_yAxis(
            gridLineInterpolation = "circle",
            lineWidth = 0,
            min=-1,max=1
          ) %>%
          highcharter::hc_series(
            list(
              name = paste0("correlation to junction ", junct.names[2]),
              data = df.2$junction_1,
              pointPlacement = "on",
              type = "line",
              color = "darkred"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[3]),
              data = df.2$junction_2,
              pointPlacement = "on",
              type = "line",
              color = "red"
            )
            ,
            list(
              name = paste0("correlation to junction ", junct.names[4]),
              data = df.2$junction_3,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            ),
            list(
              name = paste0("correlation to junction ", junct.names[5]),
              data = df.2$junction_4,
              pointPlacement = "on",
              type = "line",
              color = "blue"
            )
          )
      }
    }
    return(p)
  }
