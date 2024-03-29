% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8-plot-scores.R
\name{prepareDataForPlot}
\alias{prepareDataForPlot}
\title{Prepare Data For Plot}
\usage{
prepareDataForPlot(data, scores = NULL, runs = NULL)
}
\arguments{
\item{data}{data frame, output of \code{\link{computeSigScores}}}

\item{scores}{(optional) character vector, indicating
the summary score(s) to plot from \code{data}}

\item{runs}{(optional) numeric vector, indicating the
repeats to plot from \code{data}}
}
\value{
A data frame in long format with 4 columns:
\describe{
  \item{\code{run}}{contains an integer run IDs; \code{0}
  corresponds to the signature, IDs greater than \code{0}
  correspond to the random signatures}
  \item{\code{sampleID}}{contains the sample ID}
  \item{\code{score}}{contains the summary score for each sample}
  \item{\code{summaryScore}}{contains summary score name}
}
}
\description{
This function reshape the input data frame
to create a long format version usable by ggplot2.
}
\author{
Alessandro Barberis
}
\keyword{internal}
