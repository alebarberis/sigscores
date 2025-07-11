% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7-significance.R
\name{computeSignificance}
\alias{computeSignificance}
\title{Compute Significance of a Summary Statistic}
\usage{
computeSignificance(
  obs,
  rnd,
  type = c("asl", "ci", "se"),
  ci = c("standard", "percentile"),
  conf.level = 0.95,
  ...
)
}
\arguments{
\item{obs}{a numeric scalar representing the observed value of the summary statistic.}

\item{rnd}{numerical vector containing the permutation/bootstrap
statistics.}

\item{type}{character string indicating the significance to
compute. Three options are available:
\describe{
\item{\code{asl}}{achieved significance level}
\item{\code{ci}}{confidence interval}
\item{\code{se}}{standard error}
}}

\item{ci}{character string indicating the type of
confidence interval. Two options are available:
\describe{
\item{\code{standard}}{the standard confidence interval}
\item{\code{percentile}}{the percentile confidence interval}
}}

\item{conf.level}{numeric scalar, the desired confidence level. Used if \code{type = "ci"}.}

\item{...}{further arguments passed to \code{\link{computeASL}},
\code{\link{computeCI}}, or \code{\link{computeSE}}. For example,
\code{alternative} for \code{type = "asl"}, \code{distribution} or \code{n}
for \code{type = "ci"} with \code{ci = "standard"}.}
}
\description{
Compute Significance of a Summary Statistic
}
\seealso{
\code{\link{computeASL}},
\code{\link{computeCI}},
\code{\link{computeSE}}
}
\author{
Alessandro Barberis
}
