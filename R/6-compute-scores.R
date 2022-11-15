#'@include 0-utility-functions.R 5-scorers.R
NULL

#'Compute the Summary Scores
#'
#'@description This function computes summary
#'score(s) of the signature \code{i} in input
#'considering each column vector in the input matrix
#'\code{x}.
#'
#'A parallel execution to speed up the computation
#'on a multi-core machine can be run by setting
#'the argument \code{cores} with a number greater
#'than \code{1}.
#'
#'@inheritParams computeScores
#'@param sampling character string, indicating whether
#'to compute the scores using the provided data
#'(\code{sampling = "none"}, default), whether
#'to sample the data (\code{sampling = "permutation"}
#'and \code{sampling = "bootstrap"}),
#'or whether to generate random signatures, i.e.
#'vectors the same size of \code{i} with values
#'randomly assigned from the possible values in
#'\code{x}.
#'
#'Five options are available:
#' \describe{
#'   \item{\code{none}}{use \code{x} as it is}
#'   \item{\code{permutation}}{random sampling without replacement from row elements of \code{x}}
#'   \item{\code{bootstrap}}{random sampling with replacement from row elements of \code{x}}
#'   \item{\code{rndsig}}{random signatures of same length of \code{i} generated
#'   from all possible values in \code{x}}
#'   \item{\code{rndsigsub}}{random signatures of same length of \code{i} generated
#'   from all possible values in \code{x} after removing \code{i} values}
#' }
#'
#'See \code{\link{sampleData}} and \code{\link{randomSignatures}}
#'for further details
#'@param n.repeat integer, number of repeated samples
#'to generate
#'
#'@return A data frame containing the computed
#'score(s) for each sample. Each row corresponds to
#'a different sample.
#'
#'If \code{sampling = "random"}, \code{sampling = "bootstrap"},
#'\code{sampling = "rndsig"} or \code{sampling = "rndsigsub"},
#'the data frame contains a column with the run information.
#'
#'The two columns containing the run/sample information are:
#'
#'\describe{
#'   \item{sampleID}{the name of the sample}
#'   \item{run}{integer indicating in which run -
#'       out of the \code{n.repeat} - was computed the score}
#'}
#'
#'@author Alessandro Barberis
#'
#'@details
#'\code{\link{computeSigScores}} uses internally
#'\code{\link{computeScores}} to handle the computation of
#'the scores.
#'
#'The scoring functions are created via the generic scorer
#'\code{\link{computeScore}}:
#'\describe{
#' \item{\code{"sum"         }}{\code{\link{sumScore}         }}
#' \item{\code{"weightedSum" }}{\code{\link{weightedSumScore} }}
#' \item{\code{"mean"        }}{\code{\link{meanScore}        }}
#' \item{\code{"trimmedMean" }}{\code{\link{trimmedMeanScore} }}
#' \item{\code{"weightedMean"}}{\code{\link{weightedMeanScore}}}
#' \item{\code{"median"      }}{\code{\link{medianScore}      }}
#' \item{\code{"mode"        }}{\code{\link{modeScore}        }}
#' \item{\code{"midrange"    }}{\code{\link{midrangeScore}    }}
#' \item{\code{"midhinge"    }}{\code{\link{midhingeScore}    }}
#' \item{\code{"trimean"     }}{\code{\link{trimeanScore}     }}
#' \item{\code{"iqr"         }}{\code{\link{iqrScore}         }}
#' \item{\code{"iqm"         }}{\code{\link{iqmScore}         }}
#' \item{\code{"mad"         }}{\code{\link{madScore}         }}
#' \item{\code{"aad"         }}{\code{\link{aadScore}         }}
#' \item{\code{"ssgsea"      }}{\code{\link{ssgseaScore}      }}
#' \item{\code{"gsva"        }}{\code{\link{gsvaScore}        }}
#' \item{\code{"plage"       }}{\code{\link{plageScore}       }}
#' \item{\code{"zscore"      }}{\code{\link{zscoreScore}      }}
#'}
#'
#'\code{\link{computeScore}} handles an input matrix
#'by calling \code{\link{computeColMeasures}}.
#'Internally, \code{\link{computeColMeasures}} uses these
#'functions to compute the measures:
#'
#'\describe{
#' \item{\code{"sum"         }}{\code{\link{colSummations}}}
#' \item{\code{"weightedSum" }}{\code{\link{colWeightedSums}}}
#' \item{\code{"mean"        }}{\code{\link{colArithmeticMeans}}}
#' \item{\code{"trimmedMean" }}{\code{\link{colTrimmedMeans}}}
#' \item{\code{"weightedMean"}}{\code{\link{colWeightedArithmeticMeans}}}
#' \item{\code{"median"      }}{\code{\link{colMidpoints}}}
#' \item{\code{"mode"        }}{\code{\link{colModes}}}
#' \item{\code{"midrange"    }}{\code{\link{colMidranges}}}
#' \item{\code{"midhinge"    }}{\code{\link{colMidhinges}}}
#' \item{\code{"trimean"     }}{\code{\link{colTrimeans}}}
#' \item{\code{"iqr"         }}{\code{\link{colIQRs}}}
#' \item{\code{"iqm"         }}{\code{\link{colIQMs}}}
#' \item{\code{"mad"         }}{\code{\link{colMADs}}}
#' \item{\code{"aad"         }}{\code{\link{colAADs}}}
#' \item{\code{"ssgsea"      }}{\code{\link{colSsgsea}}}
#' \item{\code{"gsva"        }}{\code{\link{colGsva}}}
#' \item{\code{"plage"       }}{\code{\link{colPlage}}}
#' \item{\code{"zscore"      }}{\code{\link{colZscore}}}
#'}
#'Look at the different functions to know which specific
#'arguments they accept (arguments can be passed via the
#'\code{args} parameter).
#'
#'\code{\link{computeColMeasures}} also accepts a
#'transformation function via the \code{transform.fun}
#'argument, which is used to transform the data
#'before the computation of the scores so that:
#'\code{x = transform.fun(x = x, transform.args)},
#'where \code{transform.args} is a list of parameters passed
#'to the transformation function.
#'See \code{\link{computeColMeasures}} for further details.
#'A transformation function and related arguments can be
#'passed via the \code{args} parameter (see **Examples**).
#'
#'
#'The functions used for random sampling are:
#'\describe{
#' \item{\code{"permutation"}}{\code{\link{sampleData}}}
#' \item{\code{"bootstrap"  }}{\code{\link{sampleData}}}
#' \item{\code{"rndsig"     }}{\code{\link{randomSignatures}}}
#' \item{\code{"rndsigsub"  }}{\code{\link{randomSignatures}}}
#' }
#'
#'@seealso
#'Use \code{\link{getAvailableScores}} to list the available
#'built-in scores.
#'
#'Use \code{\link{getAvailableDataTransformers}} to list the available
#'built-in data transformers
#'
#'
#'@examples
#'\dontrun{
#'#set seed for reproducibility
#'set.seed(seed = 5381L)
#'
#'#Define row/col size
#'nr = 20
#'nc = 10
#'
#'#Create input matrix
#'x = matrix(
#'  data = stats::runif(n = nr*nc, min = 0, max = 1000),
#'  nrow = nr,
#'  ncol = nc,
#'  dimnames = list(
#'    paste0("g",seq(nr)),
#'    paste0("S",seq(nc))
#'  )
#')
#'
#'#Compute all scores
#'computeSigScores(
#'  x = x,
#'  i = rownames(x)[1:10]
#')
#'
#'#Compute one score
#'computeSigScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scores = 'mean'
#')
#'
#'#Compute one score passing an argument
#'computeSigScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scores = 'trimmedMean',
#'  args = list(trimmedMean = list(trim = 0.2))
#')
#'
#'#Transform data and compute the scores
#'computeSigScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'   'score1' = getScorer('weightedSum'),
#'   'score2' = getScorer('trimmedMean')
#'  ),
#'  args = list(
#'   'score1' = list(transform.fun = getDataTransformer('quantile')),
#'   'score2' = list(
#'      trim = 0.2,
#'      transform.fun = getDataTransformer('stepFunction'),
#'      transform.args = list(
#'        method = 'median',
#'        by = 'rows'
#'      )
#'    )
#'  )
#')
#'
#'#Compute scores with permutation
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "permutation",
#'  n.repeat = 10
#')
#'
#'#Compute scores with bootstrap
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "bootstrap",
#'  n.repeat = 10
#')
#'
#'#Compute scores with random signatures
#'#(elements of i are possible)
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "rndsig",
#'  n.repeat = 10
#')
#'
#'#Compute scores with random signatures
#'#(elements of i are excluded)
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "rndsigsub",
#'  n.repeat = 10
#')
#'
#'}
#'
#'
#'@export
computeSigScores <- function(
    x,
    i      = NULL,
    na.rm  = TRUE,
    scores = c("sum", "weightedSum",
               "mean", "trimmedMean", "weightedMean",
               "median", "mode", "midrange", "midhinge",
               "trimean", "iqr", "iqm", "mad", "aad",
               "ssgsea", "gsva", "plage", "zscore"),
    scorers   = NULL,
    args      = NULL,
    sampling  = c("none", "permutation", "bootstrap", "rndsig", "rndsigsub"),
    n.repeat  = 10L,
    cores     = 1L
){

  #check input ------------------------------------------------
   ##sampling
  sampling = match.arg(sampling)

  #compute scores ------------------------------------------------
  if(isTRUE(identical(sampling, "none"))){
    ##compute scores
    out = computeScores(
      x         = x,
      i         = i,
      na.rm     = na.rm,
      scores    = scores,
      scorers   = scorers,
      sample.id = T,
      cores     = cores,
      args      = args
    )

  } else {

    if(isTRUE(sampling %in% c("permutation", "bootstrap"))){
      ##random  data ------------------------------------------------
      ##sample
      rnd = sampleData(
        x      = x,
        n      = n.repeat,
        method = sampling
      )

      ##compute scores
      out = forLoop(
        n.iter   = n.repeat,
        cores    = cores,
        .inorder = F,
        fun = function(iloop, x, rnd, i, args, ...) {
          #get rnames
          rnames = rownames(x)
          #permute
          x = x[rnd[[iloop]],,drop=F]
          #set rnames
          rownames(x) = rnames
          #compute
          o = computeScores(
            x         = x,
            i         = i,
            sample.id = TRUE,
            args      = args,
            ...
          )
          #add run
          o = cbind(run = iloop, o)
          #return
          return(o)
        },
        x       = x,
        rnd     = rnd,
        na.rm   = na.rm,
        scores  = scores,
        scorers = scorers,
        i       = i,
        args    = args
      )
    } else if(isTRUE(sampling %in% c("rndsig", "rndsigsub"))){
      ##random signatures  ------------------------------------------
      exclude = switch(
        sampling,
        "rndsig"    = FALSE,
        "rndsigsub" = TRUE
      )

      ##sample
      rnd = randomSignatures(
        x       = x,
        i       = i,
        n       = n.repeat,
        exclude = exclude
      )

      ##compute scores
      out = forLoop(
        n.iter   = n.repeat,
        cores    = cores,
        .inorder = F,
        fun = function(iloop, x, rnd, args, ...) {
          #get rnames
          rnames = rownames(x)
          #set rnames
          rownames(x) = rnames
          #compute
          o = computeScores(
            x         = x,
            i         = rnd[[iloop]],
            sample.id = TRUE,
            args      = args,
            ...
          )
          #add run
          o = cbind(run = iloop, o)
          #return
          return(o)
        },
        x       = x,
        rnd     = rnd,
        na.rm   = na.rm,
        scores  = scores,
        scorers = scorers,
        args    = args
      )

    } else {
      stop("Error: 'sampling' is not valid. Please, check your input.\n")
    }

    ##shape as data frame
    out = data.table::setDF(x = data.table::rbindlist(l = out, use.names = TRUE, idcol = FALSE))
  }

  #return         ---------------------------------------------------
  return(out)
}

