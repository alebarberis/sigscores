% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-parallel-functions.R
\name{forLoop}
\alias{forLoop}
\title{Looper}
\usage{
forLoop(
  n.iter,
  cores = 1L,
  .inorder = TRUE,
  fun,
  ...,
  logger = NULL,
  verbose = TRUE
)
}
\arguments{
\item{n.iter}{number of iterations of the loop}

\item{cores}{number of cores to use for parallel execution.}

\item{.inorder}{logical flag indicating whether the \code{.combine}
function requires the task results to be combined in the same order
that they were submitted.  If the order is not important, then it
setting \code{.inorder} to \code{FALSE} can give improved performance.
The default value is `TRUE.}

\item{fun}{A function to be executed inside the loop. It should contain an \code{iloop} parameter.}

\item{...}{further arguments to \code{fun}}

\item{logger}{a \code{\link{Logger}}}

\item{verbose}{logical, indicating whether R should report extra information on progress}
}
\description{
Common interface to run a for loop sequentially or in parallel.
The core of the loop is defined by the function \code{fun}.
}
\author{
Alessandro Barberis
}
\keyword{internal}
