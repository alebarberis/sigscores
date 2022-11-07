#'@include utility-functions.R 1-matrix-scores.R 0-vector-scores.R
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
#'(\code{sampling = "none"}, default) or whether
#'to sample the data.
#'
#'Three options are available:
#' \describe{
#'   \item{\code{none}}{use \code{x} as it is}
#'   \item{\code{permutation}}{random sampling without replacement}
#'   \item{\code{bootstrap}}{random sampling with replacement}
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
#'If \code{sampling = "random"} or \code{sampling = "bootstrap"},
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
#'}
#'
#' @seealso
#'\code{\link{getAvailableScores}},
#'\code{\link{meanScores}},
#'\code{\link{trimmedMeanScores}},
#'\code{\link{weightedMeanScores}},
#'\code{\link{interquartileMeanScores}},
#'\code{\link{medianScores}},
#'\code{\link{modeScores}},
#'\code{\link{midrangeScores}},
#'\code{\link{midhingeScores}},
#'\code{\link{trimeanScores}},
#'\code{\link{bristowScores}},
#'\code{\link{reviewedBristowScores}},
#'\code{\link{weightedSumScores}},
#'\code{\link{ssGseaScores}},
#'\code{\link{gsvaScores}},
#'\code{\link{plageScores}},
#'\code{\link{zScores}},
#'\code{\link{sampleData}}
computeSigScores <- function(
    x,
    i      = NULL,
    na.rm  = TRUE,
    scores = c("mean", "trimmedMean", "weightedMean", "median",
               "mode", "midrange", "midhinge",
               "trimean", "bristow", "reviewedBristow",
               "IQM", "weightedSum",
               "ssGSEA", "gsva", "plage", "zscore"),
    sampling  = c("none", "permutation", "bootstrap"),
    n.repeat  = 10L,
    cores     = 1L,
    args      = NULL
){

  #check input ------------------------------------------------
  ##x
  if(isFALSE(is.matrix(x))){
    x = matrix(data = x, ncol = 1, dimnames = list(names(x)))
  }
  ##sample names
  colnames(x) = colnames(x = x, do.NULL = FALSE, prefix = "S")
  ##gene names
  rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")
  ##scores
  if(isTRUE(missing(scores) | is.null(scores) | is.na(scores))){
    scores = getAvailableScores()$id
  } else {
    ###intersection
    scores = intersect(scores, getAvailableScores()$id)
    ###check match
    if(isFALSE(length(scores)>0)){
      stop("Error: provided 'scores' not valid. Please, check available scores using 'getAvailableScores()'.\n")
    }
  }

  ##signature
  if(isTRUE(missing(i) | is.null(i))){
    i = rownames(x)
  }
  ###get size
  n = length(i)

  ##sampling
  sampling = match.arg(sampling)

  #compute scores ------------------------------------------------
  if(isTRUE(identical(sampling, "none"))){
    ##compute scores
    out = computeMatrixScores(
      x         = x,
      i         = i,
      na.rm     = na.rm,
      scores    = scores,
      sample.id = T,
      cores     = cores,
      args      = args
    )

  } else {

    ##create random  ------------------------------------------------
    ##sample
    rnd = sampleData(
      x      = x,
      n      = n.repeat,
      method = sampling
    )

    #compute sampling scores ---------------------------------------
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
        o = computeMatrixScores(
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
      x      = x,
      rnd    = rnd,
      na.rm  = na.rm,
      scores = scores,
      i      = i,
      args   = args
    )

    ##shape as data frame
    out = data.table::setDF(x = data.table::rbindlist(l = out, use.names = TRUE, idcol = FALSE))
  }

  #return         ------------------------------------------------
  return(out)
}


