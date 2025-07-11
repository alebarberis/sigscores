% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8-plot-scores.R
\name{heatmapSigScores}
\alias{heatmapSigScores}
\title{Plot Summary Scores as Heatmap}
\usage{
heatmapSigScores(
  data,
  scores = NULL,
  name = "summary\\nscores",
  column_title = "Summary Scores",
  ...
)
}
\arguments{
\item{data}{data frame, output of \code{\link{computeSigScores}}}

\item{scores}{(optional) character vector, indicating
the summary score(s) to plot from \code{data}}

\item{name}{Name of the heatmap. By default the heatmap name is used as the title of the heatmap legend.}

\item{column_title}{Title on the column.}

\item{...}{further arguments to internal function \code{\link{heatMap}}}
}
\value{
A \code{ComplexHeatmap} object.
}
\description{
This function generates a
\code{ComplexHeatmap} object.
}
\examples{
#Set seed for reproducibility
set.seed(seed = 5381L)

#Define row/col size
n = 10

#Create input matrix
x = matrix(
 data = stats::runif(n = n*n, min = 0, max = 1000),
 nrow = n,
 ncol = n,
 dimnames = list(
    paste0("g",seq(n)),
    paste0("S",seq(n))
 )
)

#Compute Summary Scores
x = computeSigScores(
 x = x,
 i = rownames(x)
)

#Plot scores
heatmapSigScores(data = x)

}
\seealso{
\code{\link{heatmapCorSigScores}}
}
\author{
Alessandro Barberis
}
