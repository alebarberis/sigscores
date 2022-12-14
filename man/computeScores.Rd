% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5-scorers.R
\name{computeScores}
\alias{computeScores}
\title{Compute Scores}
\usage{
computeScores(
  x,
  i = NULL,
  na.rm = TRUE,
  scorers = NULL,
  scores = NULL,
  args = NULL,
  sample.id = TRUE,
  cores = 1L,
  logger = NULL,
  filename = "sigscores",
  outdir = NULL
)
}
\arguments{
\item{x}{features-by-samples matrix}

\item{i}{(optional) numerical vector giving the rows
in \code{x} or character vector matching the row
names in \code{x}
If \code{missing} or \code{i = NULL}, all the rows
in \code{x} are considered for the computation of
the scores}

\item{na.rm}{logical, whether to remove \code{NA}
values before computation}

\item{scorers}{named list of scoring functions.
If provided, \code{scores} is not considered.
Each function must accept some specific arguments,
i.e. \code{x}, \code{i}, \code{na.rm}, \code{...}
and is expected to compute a score for each column
in \code{x}}

\item{scores}{(optional) character vector, indicating the
summary score(s) to compute}

\item{args}{named list, where the names must match the
\code{scores} or the names of \code{scorers}.
Each element in the list is another list
containing the arguments to pass to the function used
for computing the named score. For example,
\code{args = list(trimmedMean = list(trim = 0.4))}
indicates to use \code{trim = 0.4} when computing the
trimmed mean scores (\code{scores = "trimmedMean"} or
\code{scorers = list(trimmedMean = getScorer("trimmedMean"))})}

\item{sample.id}{logical, whether to report the
sample ID as a column in the output data frame}

\item{cores}{number of cores to use for parallel execution.}

\item{logger}{(optional) a \code{\link{Logger}} object.
If provided, it will be used to report extra information
on progress. To create a Logger use \code{\link{createLogger}}}

\item{filename}{(optional) character string, a name without
extension for the output file}

\item{outdir}{(optional) character string, path to
the output directory. If provided the returned data
will be stored}
}
\value{
A data frame containing the computed
score(s) for each sample. Each row corresponds to
a different sample.
}
\description{
This function computes summary
score(s) of the signature \code{i} in input
considering each column vector in the input matrix
\code{x}.
}
\details{
\code{\link{computeScores}} uses the provided scorers
to calculate the measures.

If \code{scorers = NULL} and \code{scores} are provided,
the scorers are retrieved by using the function
\code{\link{getScorers}}:

\describe{
\item{\code{"sum"         }}{\code{\link{sumScorer}}}
\item{\code{"weightedSum" }}{\code{\link{weightedSumScorer}}}
\item{\code{"mean"        }}{\code{\link{meanScorer}}}
\item{\code{"trimmedMean" }}{\code{\link{trimmedMeanScorer}}}
\item{\code{"weightedMean"}}{\code{\link{weightedMeanScorer}}}
\item{\code{"median"      }}{\code{\link{medianScorer}}}
\item{\code{"mode"        }}{\code{\link{modeScorer}}}
\item{\code{"midrange"    }}{\code{\link{midrangeScorer}}}
\item{\code{"midhinge"    }}{\code{\link{midhingeScorer}}}
\item{\code{"trimean"     }}{\code{\link{trimeanScorer}}}
\item{\code{"iqr"         }}{\code{\link{iqrScorer}}}
\item{\code{"iqm"         }}{\code{\link{iqmScorer}}}
\item{\code{"mad"         }}{\code{\link{madScorer}}}
\item{\code{"aad"         }}{\code{\link{aadScorer}}}
\item{\code{"ssgsea"      }}{\code{\link{ssgseaScorer}}}
\item{\code{"gsva"        }}{\code{\link{gsvaScorer}}}
\item{\code{"plage"       }}{\code{\link{plageScorer}}}
\item{\code{"zscore"      }}{\code{\link{zscoreScorer}}}
}

Look at the different functions to know which specific
arguments they accept (arguments can be passed via the
\code{args} parameter).

Scorers also accepts a transformation function
via the \code{transform.fun} argument, which
is used to transform the data before the computation
of the scores so that:
\code{x = transform.fun(x = x, transform.args)},
where \code{transform.args} is a list of parameters passed
to the transformation function.
Look at the different functions for further details.
}
\examples{
\dontrun{
#set seed for reproducibility
set.seed(seed = 5381L)

#Define row/col size
nr = 20
nc = 10

#Create input matrix
x = matrix(
 data = stats::runif(n = nr*nc, min = 0, max = 1000),
 nrow = nr,
 ncol = nc,
 dimnames = list(
   paste0("g",seq(nr)),
   paste0("S",seq(nc))
 )
)

#Compute all scores
computeScores(
 x = x,
 i = rownames(x)[1:10]
)

#Compute all scores using 'scorers' argument
computeScores(
 x = x,
 i = rownames(x)[1:10],
 scorers = list(
    'score1' = getScorer('sum'),
    'score2' = getScorer('mean')
 )
)

#Pass different parameters to the same scorer
computeScores(
 x = x,
 i = rownames(x)[1:10],
 scorers = list(
  'score1' = getScorer('trimmedMean'),
  'score2' = getScorer('trimmedMean')
 ),
 args = list(
  'score1' = list(trim = 0),
  'score2' = list(trim = 0.3)
 )
)

#Transform data and compute the scores
computeScores(
 x = x,
 i = rownames(x)[1:10],
 scorers = list(
  'score1' = getScorer('weightedSum'),
  'score2' = getScorer('mean')
 ),
 args = list(
  'score1' = list(transform.fun = getDataTransformer('quantile')),
  'score2' = list(
     transform.fun = getDataTransformer('stepFunction'),
     transform.args = list(
       method = 'median',
       by = 'rows'
     )
   )
 )
)
}

}
\author{
Alessandro Barberis
}
\keyword{internal}
