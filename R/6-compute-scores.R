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
#'@param x numerical matrix, features-by-samples
#'@param i (optional) numerical vector giving the rows
#'in \code{x} or character vector matching the row
#'names in \code{x}
#'If \code{missing} or \code{i = NULL}, all the rows
#'in \code{x} are considered for the computation of
#'the score
#'@param na.rm logical, whether to remove \code{NA}
#'values before computation
#'@param scores character vector, indicating the
#'summary score(s) to compute
#'@inheritParams forLoop
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
#'@param n.repeat integer, number of repeated samples
#'to generate
#@param ... further arguments to the function used
#for computing the scores
#'@param args named list, where the names must match the
#'\code{scores}. Each element in the list is another list
#'containing the arguments to pass to the function used
#'for computing the named score. For example,
#'\code{args = list(trimmedMean = list(trim = 0.4))}
#'indicates to use \code{trim = 0.4} when computing the
#'trimmed mean score
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
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "rndsig",
#'  n.repeat = 10
#')
#'computeSigScores(
#'  x        = x,
#'  i        = rownames(x)[1:10],
#'  sampling = "rndsigsub",
#'  n.repeat = 10
#')
#'
#'}
#'
#' @seealso
#'\code{\link{getAvailableScores}},
#'\code{\link{computeScores}}
#'
#'The generic scorer:
#'\code{\link{computeScore}}
#'
#'The built-in functions used to compute the measures:
#'\code{\link{colSummations}},
#'\code{\link{colWeightedSums}},
#'\code{\link{colArithmeticMeans}},
#'\code{\link{colTrimmedMeans}},
#'\code{\link{colWeightedArithmeticMeans}},
#'\code{\link{colMidpoints}},
#'\code{\link{colModes}},
#'\code{\link{colMidranges}},
#'\code{\link{colMidhinges}},
#'\code{\link{colTrimeans}},
#'\code{\link{colIQRs}},
#'\code{\link{colIQMs}},
#'\code{\link{colMADs}},
#'\code{\link{colAADs}},
#'\code{\link{colSsgsea}},
#'\code{\link{colGsva}},
#'\code{\link{colPlage}},
#'\code{\link{colZscore}}
#'
#'The functions used for random sampling:
#'\code{\link{sampleData}},
#'\code{\link{randomSignatures}}
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
      scorers   = NULL,
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


