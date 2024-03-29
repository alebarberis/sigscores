% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3-stats-functions.R
\name{IQM}
\alias{IQM}
\title{Interquartile Mean}
\usage{
IQM(x, i = NULL, na.rm = T)
}
\arguments{
\item{x}{(named) numerical vector}

\item{i}{(optional) numerical vector giving the position in \code{x}
or character vector matching the names in \code{x}.
If \code{missing} or \code{i = NULL}, the entire \code{x} is
considered for the computation of the score}

\item{na.rm}{unused argument, provided for consistency
with other functions.
\code{NA} values are always removed from \code{x} before
computation}
}
\value{
A numerical value representing the computed measure.
A default \code{NA} value is returned if the score can't be
computed.
}
\description{
This function computes the *interquartile mean*
of the \code{x} values.
See the **Details** section below for further information.
}
\details{
The *interquartile mean* is a statistical measure
of central tendency based on the truncated mean of the
interquartile range. It is computed as:

\deqn{IQM(x) = \frac{2}{n}\sum_{i=\frac{n}{4}+1}^{\frac{3n}{4}}  x_{i}}

where \eqn{x_{i}} is the \eqn{i}-th element of the ordered vector.

Like the median it is insensitive to outliers.
}
\examples{
#Dataset size divisible by four
x = c(5,8,4,38,8,6,9,7,7,3,1,6)
IQM(x)#6.5

#Dataset size not divisible by four
x = c(1,2,3,4,5)
IQM(x = x)#3

x = c(1,3,5,7,9,11,13,15,17)
IQM(x)#9

}
\references{
https://en.wikipedia.org/wiki/Interquartile_mean
}
\author{
Alessandro Barberis
}
