% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DJEannotate.R
\name{DJEannotate}
\alias{DJEannotate}
\title{DJEannotate: Gene annotation for junctions used in differential expression analysis}
\usage{
DJEannotate(import.out, gtf)
}
\arguments{
\item{import.out}{output object from DJEimport()}

\item{gtf}{Reference transcriptome in genecode gtf file format. Used to define gene ID for junctions}
}
\value{
List object containing junction expression and assigned gene IDs
}
\description{
Matches junction ID in expression matrix with respective gene name and returns output object for DJEprepare.
}
\examples{
DJEimp <- system.file("extdata", "DJEimp.rds", package = "DJExpress")
imp.out <- readRDS(DJEimp)
gtf0 <- system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
ann.out <- DJEannotate(imp.out, gtf0)
}
