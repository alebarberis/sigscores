% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5-scorers.R
\name{midrangeScorer}
\alias{midrangeScorer}
\title{Midrange Scorer}
\usage{
midrangeScorer(
  x,
  i = NULL,
  na.rm = TRUE,
  transform.fun = NULL,
  transform.args = list(),
  transform.sub = F
)
}
\arguments{
\item{x}{a numerical matrix features-by-samples}

\item{i}{(optional) numerical vector giving the (row) position
in \code{x} or character vector matching the (row) names in \code{x}.
If \code{missing} or \code{i = NULL}, the entire \code{x} is
considered for the computation of the score}

\item{na.rm}{logical, whether to remove \code{NA}
values from \code{x} before computation}

\item{transform.fun}{function to transform the data.
If provided, \code{x} is transformed
(\code{x = transform.fun(x, transform.args)})
before the computation of the scores.
See \code{\link{getDataTransformer}} for further details
about built-in options}

\item{transform.args}{list of parameters to the
data transformation function}

\item{transform.sub}{logical, whether to transform
\code{x} after it is subset for \code{i} (used to
speedup computation). Default is \code{FALSE},
meaning the transformation would be applied directly
to \code{x} provided in input}
}
\value{
A numerical value representing the computed score.
A default \code{NA} value is returned if the score can't be
computed, or if \code{i} values are not present in \code{x}.
}
\description{
This scorer returns the *midrange* score(s).
See the **Details** section below for further information.
}
\details{
The midrange score is the average of the lowest
and highest values in a set of data.
It is a measure of central tendency like
the mean, median, and mode. However, it is more prone to bias than
these other measures because it relies solely upon the two most
extreme scores, which could potentially be outliers.

This function calculates the range of values in each
column of \code{x} by calling \code{\link[matrixStats]{colRanges}}.
Then, it uses the min and max values to compute the score.
}
\seealso{
\code{\link[matrixStats]{colRanges}}
}
\author{
Alessandro Barberis
}
