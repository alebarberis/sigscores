% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8-plot-scores.R
\name{scatterplotSigScores}
\alias{scatterplotSigScores}
\title{Plot Summary Scores as Scatter Plots}
\usage{
scatterplotSigScores(
  data,
  scores = NULL,
  runs = NULL,
  add.points = T,
  point.shape = 16,
  point.jitter.w = 0.2,
  point.jitter.h = NULL,
  labs.title = "Summary Scores",
  labs.x = "Sample IDs",
  labs.y = "Scores",
  labs.col = "Summary\\nScores",
  axis.text.x.angle = 90,
  axis.text.x.size = 9,
  axis.text.y.angle = 0,
  axis.text.y.size = 9,
  ...
)
}
\arguments{
\item{data}{data frame, output of \code{\link{computeSigScores}}}

\item{scores}{(optional) character vector, indicating
the summary score(s) to plot from \code{data}}

\item{runs}{(optional) numeric vector, indicating the
repeats to plot from \code{data}}

\item{add.points}{logical, whether to add the computed scores
of the individual summary measures as points in the plot}

\item{point.shape}{the shape to use to plot the scores. It can take
five types of values:
\itemize{
  \item An integer in `[0, 25]`
  \item The name of the shape
  \item A single character, used as a plotting symbol
  \item A . to draw the smallest rectangle that is visible,
  usually 1 pixel
  \item An \code{NA}, to draw nothing
}

See \code{vignette("ggplot2-specs")} for further details}

\item{point.jitter.w, point.jitter.h}{Amount of vertical and
horizontal jitter. The jitter is added in both positive and negative
directions, so the total spread is twice the value specified here.
See \code{\link[ggplot2]{position_jitter}} for further details}

\item{labs.title}{The text for the title}

\item{labs.x}{The title of the x axis}

\item{labs.y}{The title of the y axis}

\item{labs.col}{The title of the legend}

\item{axis.text.x.angle, axis.text.y.angle}{Specify the x and y
axis tick labels angles (in \code{[0, 360]})
See \code{\link[ggplot2]{element_text}} for further details}

\item{axis.text.x.size, axis.text.y.size}{Specify the x and y
axis tick labels size in pts.
See \code{\link[ggplot2]{element_text}} for further details}

\item{...}{further arguments to \code{\link{ggPlot}}}
}
\value{
A \code{ggplot} object.
}
\description{
This function generates a
\code{ggplot} object.

\code{scatterplotSigScores} uses the internal function
\code{\link{ggPlot}} which contains different function
calls to \code{ggplot2} functions in order to create a
pre-defined plot.
}
\examples{
#Set seed for reproducibility
set.seed(seed = 5381L)

#Define row/col size
n = 10

#Create input matrix
x = matrix(
 data = stats::runif(n = n*n, min = 0, max = 100),
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
scatterplotSigScores(data = x)

#Plot scatter plots per summary score
scatterplotSigScores(
 data       = x,
 scores     = c("mean", "median"),
 facet.rows = "summaryScore"
)

}
\seealso{
\code{\link{ggPlot}}
}
\author{
Alessandro Barberis
}
\keyword{internal}
