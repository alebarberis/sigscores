#'@include 0-utility-functions.R 1-log-functions.R 2-parallel-functions.R 3-stats-functions.R 4-data-preprocessing.R
NULL

# Generic Scorer ----------------------------------------------------------


#'Generic Scorer Function
#'
#'@description This function computes a summary score
#'from an input vector or from each column vector
#'of an input matrix.
#'
#'@param x (named) numerical vector or matrix
#'@param i (optional) numerical vector giving the (row) position
#'in \code{x} or character vector matching the (row) names in \code{x}.
#'If \code{missing} or \code{i = NULL}, the entire \code{x} is
#'considered for the computation of the score
#'@param w numerical vector of weights the same length as
#'\code{i}. It is used for the computation of the weighted scores
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'@param score character string indicating the summary score to compute
#'@param ... further arguments to \code{score}
#@param transform.fun (optional) data transformation function.
#If provided, it must accept \code{x} in input.
#The data is transformed before the computation
#of the score
#'@param transform string, the transformation method.
#'Available options are:
#'\describe{
#' \item{\code{'none'}}{\code{x} is returned}
#' \item{\code{'stepFunction'}}{step function}
#' \item{\code{'quantile'}}{quantile normalisation}
#'}
#'See \code{\link{transformData}} for further details
#'@param transform.args list of parameters to the
#'data transformation function
#'@param transform.sub logical, whether to transform
#'\code{x} after it is subset for \code{i} (used to
#'speedup computation). Default is \code{FALSE},
#'meaning the transformation would be applied directly
#'to \code{x} provided in input
#'
#'@return A numerical value representing the computed score.
#'A default \code{NA} value is returned if the score can't be
#'computed, or if \code{i} values are not present in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{computeVectorScore}},
#'\code{\link{computeColScores}}
#'
#'@keywords internal
computeScore <- function(
    x,
    i     = NULL,
    w     = NULL,
    na.rm = T,
    score  = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqr", "iqm", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore"),
    transform = c('none', 'stepFunction', 'quantile'),
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){
  score = match.arg(score)
  transform = match.arg(transform)

  if(isTRUE(is.vector(x))){
    out = computeVectorScore(x=x,i=i,w=w,na.rm=na.rm,score=score,transform=transform,transform.args=transform.args,transform.sub=transform.sub,...)
  } else if(isTRUE(is.matrix(x))){
    out = computeColScores(x=x,i=i,w=w,na.rm=na.rm,score=score,transform=transform,transform.args=transform.args,transform.sub=transform.sub,...)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}

# #'@export
# computeScore.vector <- function(x, ...){computeVectorScore(x, ...)}
# #'@export
# computeScore.matrix  <- function(x, ...){computeColScores(x, ...)}





#'Compute Vector Score
#'
#'@description This function computes a summary score
#'from an input vector
#'
#'@param x (named) numerical vector
#'@param i (optional) numerical vector giving the position in \code{x}
#'or character vector matching the names in \code{x}.
#'If \code{missing} or \code{i = NULL}, the entire \code{x} is
#'considered for the computation of the score
#'@inheritParams computeScore
#'
#'@return A numerical value representing the computed score.
#'A default \code{NA} value is returned if the score can't be
#'computed, or if \code{i} values are not present in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeVectorScore <- function(
    x,
    i,
    w     = NULL,
    na.rm = T,
    score = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqr", "iqm", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore"),
    transform = c('none', 'stepFunction', 'quantile'),
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){

  #match
  score = match.arg(score)
  transform = match.arg(transform)

  #check input
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = length(x))}
  i = updateSig(x = x, i = i)
  if(isTRUE(is.null(i))){   return(getDefaultNaValue())}

  #update
  if(isTRUE(identical(transform, "none"))){
    x = subsetData(x = x, i = i)
  } else {
    if(isTRUE(transform.sub)){
      x = subsetData(x = x, i = i)
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
    } else {
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
      x = subsetData(x = x, i = i)
    }
  }

  #compute
  if(isTRUE(score %in% c("weightedSum", "weightedMean"))){
    if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(i))}
    #update
    w = w[getIndex(i)]
    #compute
    out = computeMeasure(score = score, x = x, w = w, na.rm = na.rm, ...)
  } else {
    out = computeMeasure(score = score, x = x, na.rm = na.rm, ...)
  }


  #return
  return(out)
}




#'Compute Score for Each Column of x
#'
#'@description This function computes a summary score for
#'each column in an input matrix.
#'
#'@param x (named) numerical vector
#'@param i (optional) numerical vector giving the position in \code{x}
#'or character vector matching the names in \code{x}.
#'If \code{missing} or \code{i = NULL}, the entire \code{x} is
#'considered for the computation of the score
#'@inheritParams computeScore
#'
#'@return A numerical value representing the computed score.
#'A default \code{NA} value is returned if the score can't be
#'computed, or if \code{i} values are not present in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeColScores <- function(
  x,
  i,
  w,
  na.rm = T,
  score = c("sum", "weightedSum",
            "mean", "trimmedMean", "weightedMean",
            "median", "mode", "midrange", "midhinge",
            "trimean", "iqr", "iqm", "mad", "aad",
            "ssgsea", "gsva", "plage", "zscore"),
  transform = c('none', 'stepFunction', 'quantile'),
  transform.args = list(),
  transform.sub  = FALSE,
  ...
  ){

  #get dim
  numCol = ncol(x)

  #match
  score     = match.arg(score)
  transform = match.arg(transform)

  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = switch(
    score,
    "sum"             =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "weightedSum"     =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "mean"            =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "trimmedMean"     =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "weightedMean"    =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "median"          =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "mode"            =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "midrange"        =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "midhinge"        =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "trimean"         =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "bristow"         =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "iqr"             =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "iqm"             =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "mad"             =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "aad"             =  computeColStatScores(score = score, x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "ssgsea"          =  computeColEnrichmentScores(score = score, x = x, i = i, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "gsva"            =  computeColEnrichmentScores(score = score, x = x, i = i, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "plage"           =  computeColEnrichmentScores(score = score, x = x, i = i, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...),
    "zscore"          =  computeColEnrichmentScores(score = score, x = x, i = i, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...)
  )

  #return
  return(out)
}

computeColStatScores <- function(
    x,
    i,
    w,
    na.rm = T,
    score = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqr", "iqm", "mad", "aad"),
    transform = c('none', 'stepFunction', 'quantile'),
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){


  #get dim
  numCol = ncol(x)
  numRow = nrow(x)

  #match
  score = match.arg(score)
  transform = match.arg(transform)

  #check i
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = numRow)}
  i = updateSig(x = x, i = i)
  if(isTRUE(is.null(i))){ return(rep(x = getDefaultNaValue(), times = numCol)) }

  #update
  if(isTRUE(identical(transform, "none"))){
    x = subsetData(x = x, i = i, rm.dupli = T)
  } else {
    if(isTRUE(transform.sub)){
      x = subsetData(x = x, i = i, rm.dupli = T)
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
    } else {
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
      x = subsetData(x = x, i = i, rm.dupli = T)
    }
  }

  #compute
  if(isTRUE(score %in% c("weightedSum", "weightedMean"))){
    #check w
    if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(i))}

    #update
    w = w[getIndex(i)]
    out = computeColMeasures(score = score, x = x, w = w, na.rm = na.rm, ...)
  } else {
    out = computeColMeasures(score = score, x = x, na.rm = na.rm, ...)
  }


  #return
  return(out)
}


computeColEnrichmentScores <- function(
    x,
    i,
    na.rm = T,
    score = c("ssgsea", "gsva", "plage", "zscore"),
    transform = c('none', 'stepFunction', 'quantile'),
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){
  #get dim
  numRow = nrow(x)
  numCol = ncol(x)

  #match
  score = match.arg(score)
  transform = match.arg(transform)

  #check i
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = numRow)}
  i = updateSig(x = x, i = i)
  if(isTRUE(is.null(i))){ return(rep(x = getDefaultNaValue(), times = numCol)) }

  #update
  if(isFALSE(identical(transform, "none"))){
    x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
  }

  #compute
  out = computeColMeasures(score = score, x = x, rows = i, na.rm = na.rm, ...)

  #return
  return(out)
}


# Scorers -----------------------------------------------------------------


sumScore          <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "sum"         , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
weightedSumScore  <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "weightedSum" , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
meanScore         <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "mean"        , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
trimmedMeanScore  <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "trimmedMean" , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
weightedMeanScore <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "weightedMean", x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
medianScore       <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "median"      , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
modeScore         <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "mode"        , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
midrangeScore     <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "midrange"    , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
midhingeScore     <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "midhinge"    , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
trimeanScore      <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "trimean"     , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
iqrScore          <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "iqr"         , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
iqmScore          <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "iqm"         , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
madScore          <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "mad"         , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
aadScore          <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "aad"         , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
ssgseaScore       <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "ssgsea"      , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
gsvaScore         <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "gsva"        , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
plageScore        <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "plage"       , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}
zscoreScore       <- function(x, i = NULL, w = NULL, na.rm = TRUE, transform = 'none', transform.args = list(), transform.sub = F, ...){return(computeScore(score = "zscore"      , x = x, i = i, w = w, na.rm = na.rm, transform = transform, transform.args = transform.args, transform.sub = transform.sub, ...))}

# Get Scorer(s) -----------------------------------------------------------

#'Get Scorer
#'
#'@description This function is a dispatcher for the score
#'function selected in input.
#'
#'@param score character string, one of the supported measures
#'@return A function.
#'
#'@author Alessandro Barberis
#'
#'@export
getScorer <- function(
  score = c("sum", "weightedSum",
            "mean", "trimmedMean", "weightedMean",
            "median", "mode", "midrange", "midhinge",
            "trimean", "iqr", "iqm", "mad", "aad",
            "ssgsea", "gsva", "plage", "zscore")
  ){
  score = match.arg(score)

  #compute
  out = switch(
    score,
    "sum"             = sumScore         ,
    "weightedSum"     = weightedSumScore ,
    "mean"            = meanScore        ,
    "trimmedMean"     = trimmedMeanScore ,
    "weightedMean"    = weightedMeanScore,
    "median"          = medianScore      ,
    "mode"            = modeScore        ,
    "midrange"        = midrangeScore    ,
    "midhinge"        = midhingeScore    ,
    "trimean"         = trimeanScore     ,
    "iqr"             = iqrScore         ,
    "iqm"             = iqmScore         ,
    "mad"             = madScore         ,
    "aad"             = aadScore         ,
    "ssgsea"          = ssgseaScore      ,
    "gsva"            = gsvaScore        ,
    "plage"           = plageScore       ,
    "zscore"          = zscoreScore
  )

  return(out)
}


#'Get Scorers
#'
#'@description This function is a dispatcher for the score
#'functions selected in input.
#'
#'@param scores character vector, the selected supported measures
#'@return A list containing the scoring functions.
#'
#'@author Alessandro Barberis
#'
#'@seealso The generic scorer:
#'\code{\link{computeScore}}
#'
#'And the built-in functions used to compute the measures:
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
#'@export
getScorers <- function(
    scores = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqr", "iqm", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore")){

  #check
  isId = scores %in% getAvailableScores()$id
  scores = scores[isId]
  if(isFALSE(all(isId))){
    if(isTRUE(length(scores)>0)){
      warning("Some score ids are not valid. Please, check your input.\n")
    } else {
      stop("Error: 'scores' is not valid. Please, check your input.\n")
    }
  }

  out = lapply(X = scores, FUN = getScorer)

  names(out) = scores

  return(out)
}

# Compute Scores ----------------------------------------------------------

#'Compute Scores
#'
#'@description This function computes summary
#'score(s) of the signature \code{i} in input
#'considering each column vector in the input matrix
#'\code{x}.
#'
#'@param x features-by-samples matrix
#'@param i (optional) numerical vector giving the rows
#'in \code{x} or character vector matching the row
#'names in \code{x}
#'If \code{missing} or \code{i = NULL}, all the rows
#'in \code{x} are considered for the computation of
#'the scores
#@param w numerical vector of weights the same length as
#\code{i}
#'@param na.rm logical, whether to remove \code{NA}
#'values before computation
#'@param scores character vector, indicating the
#'summary score(s) to compute
#'@param scorers named list of scoring functions.
#'If provided, \code{scores} is not considered.
#'Each function must accept some specific arguments,
#'i.e. \code{x}, \code{i}, \code{na.rm}, \code{...}
#'and is expected to compute a score for each column
#'in \code{x}
#'@param args named list, where the names must match the
#'\code{scores} or the names of \code{scorers}.
#'Each element in the list is another list
#'containing the arguments to pass to the function used
#'for computing the named score. For example,
#'\code{args = list(trimmedMean = list(trim = 0.4))}
#'indicates to use \code{trim = 0.4} when computing the
#'trimmed mean scores
#'@param sample.id logical, whether to report the
#'sample ID as a column in the output data frame
#'@inheritParams forLoop
#'@param logger a \code{\link{Logger}}
#'
#'@return A numerical vector containing the computed
#'score for each sample.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
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
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10]
#')
#'
#'#Compute all scores using 'scorers' argument
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'     'score1' = getScorer('sum'),
#'     'score2' = getScorer('mean')
#'  )
#')
#'
#'#Pass different parameters to the same scorer
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'   'score1' = getScorer('trimmedMean'),
#'   'score2' = getScorer('trimmedMean')
#'  ),
#'  args = list(
#'   'score1' = list(trim = 0),
#'   'score2' = list(trim = 0.3)
#'  )
#')
#'
#'#Transform data and compute the scores
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'   'score1' = getScorer('weightedSum'),
#'   'score2' = getScorer('mean')
#'  ),
#'  args = list(
#'   'score1' = list(transform = 'quantile'),
#'   'score2' = list(
#'      transform = 'stepFunction',
#'      method = 'median',
#'      by = 'rows'
#'    )
#'  )
#')
#'}
computeScores <- function(
    x,
    i         = NULL,
    na.rm     = TRUE,
    scorers   = NULL,
    scores    = NULL,
    args      = NULL,
    sample.id = TRUE,
    cores     = 1L,
    logger    = NULL
  ){

  #check input ------------------------------------------------
  ##x
  if(isFALSE(is.matrix(x))){
    x = matrix(data = x, ncol = 1, dimnames = list(names(x)))
  }
  ##sample names
  if(isTRUE(is.null(colnames(x)))){colnames(x) = colnames(x = x, do.NULL = FALSE, prefix = "S")}
  ##gene names
  if(isTRUE(is.null(rownames(x)))){rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")}

  ##scorers
  if(isTRUE(missing(scorers) | is.null(scorers))){
    if(isTRUE(missing(scores) | is.null(scores))){
      scores = getAvailableScores()$id
    } else {
      ###intersection
      scores = intersect(scores, getAvailableScores()$id)
      ###check match
      if(isFALSE(length(scores)>0)){
        stop("Error: provided 'scores' not valid. Please, check available scores using 'getAvailableScores()'.\n")
      }
    }
    scorers = getScorers(scores = scores)
  } else {
    if(isTRUE(is.list(scorers))){
      if(isFALSE(all(sapply(scorers, FUN = is.function)))){
        stop("Error: provided 'scorers' not valid. Each element of the list must be a function computing a score.\n")
      }
    } else {
      stop("Error: provided 'scorers' not valid. It should be a named list, where each element is a function computing a score.\n")
    }
  }
  ##args
  if(isFALSE(missing(args) | is.null(args))){
    if(isTRUE(is.list(args) && !is.null(names(args)) && (any(names(args) %in% names(scorers))))){
      #update
      args = args[(names(args) %in% names(scorers))]
    } else{
      stop("Error: provided 'args' not valid. It must be a named list, with names matching the provided 'scorers'.\n")
    }
  }

  ##signature
  if(isTRUE(missing(i) | is.null(i))){
    i = rownames(x)
  }
  ##get size
  n = length(i)

  #compute scores ------------------------------------------------
  ##loop over scores
  ###number of iteration
  n.iter = length(scorers)
  ###compute scores
  out = forLoop(
    n.iter   = n.iter,
    cores    = cores,
    .inorder = T,
    fun = function(iloop, scorers, args, ...) {
      #compute
      o = do.call(
        what = scorers[[iloop]],
        args = c(
          list(...),
          args[[names(scorers)[[iloop]]]]
        )
      )
      #return
      return(o)
    },
    x      = x,
    na.rm  = na.rm,
    scorers = scorers,
    i      = i,
    logger = logger,
    # w      = w,
    args   = args
  )

  ##set names
  names(out) = names(scorers)

  ##shape as data frame
  out = data.frame(
    out,
    row.names = colnames(x),
    stringsAsFactors = F
  )

  #sample ID
  if(isTRUE(sample.id)){
    #add column
    out$sampleID = rownames(out)

    #re-order columns
    out = out[,unique(c("sampleID", colnames(out))), drop=F]
  }

  #return         ------------------------------------------------
  return(out)
}
