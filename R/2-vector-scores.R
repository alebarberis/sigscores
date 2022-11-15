#'@include 0-utility-functions.R 3-stats-functions.R
NULL

#'Compute Score
#'
#'@description This function computes a summary score
#'from an input vector
#'
#'@param x (named) numerical vector
#'@param i (optional) numerical vector giving the position in \code{x}
#'or character vector matching the names in \code{x}.
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
#'@keywords internal
computeScore <- function(
    x,
    i,
    w,
    na.rm = T,
    score = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqm", "iqr", "mad", "aad",
              "ssGSEA", "gsva", "plage", "zscore"),
    transform = c('none', 'stepFunction', 'quantile'),
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){

  #match
  score = match.arg(score)
  transform = match.arg(transform)

  #subset
  if(isTRUE(score %in% c("weightedSum", "weightedMean"))){
    if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = length(x))}
    if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(i))}

    #subset input
    i = updateSig(x = x, i = i)

    #check input
    if(isTRUE(is.null(i))){   return(getDefaultNaValue())}

    #update
    w = w[getIndex(i)]
  }

  #update
  if(isTRUE(identical(transform, "none"))){
    x = subsetVector(x = x, i = i)
  } else {
    if(isTRUE(transform.sub)){
      x = subsetVector(x = x, i = i)
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
    } else {
      x = do.call(what = transformData, args = c(list(x = x, f = transform), transform.args))
      x = subsetVector(x = x, i = i)
    }
  }

  #compute
  out = computeMeasure(score = score, x = x, w = w, na.rm = na.rm, ...)

  #return
  return(out)
}



#'Sum Score
#'
#'@description This function computes the *sum* of
#'values in \code{x}.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[base]{sum}}
#'
sumScore <- function(
    x,
    i,
    na.rm = T,
    transform = 'none',
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){

  #compute
  x = computeScore(
    score = "sum",
    x = x,
    i = i,
    na.rm = na.rm,
    transform = transform,
    transform.args = transform.args,
    transform.sub  = transform.sub,
    ...
  )

  #return
  return(x)
}


#'Weighted Sum Score
#'
#'@description This function computes the *weighted sum* of
#'values in \code{x}.
#'
#'@inheritParams computeScore
#'@param w numerical vector of weights the same length as
#'\code{i}
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[base]{sum}}
#'
weightedSumScore <- function(
    x,
    i,
    w,
    na.rm = T,
    transform = 'none',
    transform.args = list(),
    transform.sub  = FALSE,
    ...
){

  #compute
  x = computeScore(
    score = "weightedSum",
    x = x,
    w = w,
    i = i,
    na.rm = na.rm,
    transform = transform,
    transform.args = transform.args,
    transform.sub  = transform.sub,
    ...
  )

  #return
  return(x)
}

#'Arithmetic Mean Score
#'
#'@description This function computes the *arithmetic mean*.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[base]{mean}}
#'
meanScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  x = computeMeasure(score = "mean", x = x, na.rm = na.rm, trim = 0)

  #return
  return(x)
}



#'Trimmed Arithmetic Mean Score
#'
#'@description This function computes the *trimmed arithmetic mean*.
#'
#'@inheritParams computeScore
#'@inheritParams base::mean
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[base]{mean}}
#'
#'
trimmedMeanScore <- function(
    x,
    i,
    na.rm = TRUE,
    trim  = 0
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  x = computeMeasure(score = "mean", x = x, na.rm = na.rm, trim = trim)

  #return
  return(x)
}




#'Weighted Arithmetic Mean Score
#'
#'@description This function computes the *weighted arithmetic mean*.
#'
#'@inheritParams computeScore
#'@param w numerical vector of weights the same length as
#'\code{i}
#@inheritParams stats::weighted.mean
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[stats]{weighted.mean}}
weightedMeanScore <- function(
    x,
    i,
    w,
    na.rm = TRUE
){

  #check input
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = length(x))}
  if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(i))}

  #subset input
  i = updateSig(x = x, i = i)

  #check input
  if(isTRUE(is.null(i))){   return(getDefaultNaValue())}

  #update
  x = subsetVector(x = x, i = i)
  w = w[getIndex(i)]

  #compute
  x = computeMeasure(score = "weightedMean", x = x, na.rm = na.rm, w = w)

  #return
  return(out)
}





#'Median Scores
#'
#'@description The *median* is the value separating the lower half
#'of a set of measures from the higher half.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[stats]{median}}
#'
#'
medianScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  out = computeMeasure(score = "median", x = x, na.rm = na.rm)

  #return
  return(out)
}


#'Mode Score
#'
#'@description Compute the *mode* of a vector, i.e. the value
#' that has highest number of occurrences. If different values have
#' the same number of occurrences, the first one is reported.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'
modeScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  out = computeMeasure(score = "mode", x = x, na.rm = na.rm)

  #return
  return(out)
}


#'Midrange Score
#'
#'@description This function computes the *midrange* score, i.e.
#'the average of the lowest and highest values in a set of data.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@inherit midrange details
#'
#'@inherit midrange references
#'
#'@seealso
#'\code{\link{meanScore}},
#'\code{\link{medianScore}},
#'\code{\link{modeScore}}
midrangeScore <- function(
  x,
  i,
  na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  out = computeMeasure(score = "midrange", x = x, na.rm = na.rm)

  #return
  return(out)
}




#'Midhinge Score
#'
#'@description The *midhinge* of a set of values
#'is the mean of the first and third quartiles.
#'
#'@inheritParams computeScore
#'@param na.rm unused argument, provided for consistency
#'with other functions.
#'Missing values and NaN's are allowed
#'only if \code{na.rm = TRUE} in \code{\link[stats]{quantile}}.
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
midhingeScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  out = computeMeasure(score = "midhinge", x = x, na.rm = na.rm)

  #return
  return(out)
}



#'Trimean Score
#'
#'@description The *trimean* contains information about the
#'center and some of the position of the data.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
trimeanScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #compute
  out = computeMeasure(score = "trimean", x = x, na.rm = na.rm)

  #return
  return(out)
}



#'Interquartile Mean (IQM) Score
#'
#'@description The *interquartile mean* is a statistical measure
#'of central tendency based on the truncated mean of the
#'interquartile range. It is computed as:
#'
#'\deqn{IQM(x) = \frac{2}{n}\sum_{i=\frac{n}{4}+1}^{\frac{3n}{4}}  x_{i}}
#'
#'where \eqn{x_{i}} is the \eqn{i}-th element of the ordered vector.
#'
#'Like the median it is insensitive to outliers.
#'
#'@inheritParams computeScore
#'@param na.rm unused argument, provided for consistency
#'with other functions.
#'
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@references See https://en.wikipedia.org/wiki/Interquartile_mean
interquartileMeanScore <- function(
    x,
    i,
    na.rm   = TRUE
  ){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  out = computeMeasure(score = "iqm", x = x, na.rm = na.rm)

  #return
  return(out)
}

madScore <- function()


#'Bristow Score
#'
#'@description This function computes the *Bristow score*.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeScore
#'@param m numerical vector of length equivalent to \code{i},
#'containing the median abundance for each element of
#'\code{i} in the cohort. If \code{m = NULL}, the median
#'of \code{x} is instead used
#'
#'@inherit computeScore return
#'
#'@details To compute the score, each element in \code{i} is
#'evaluated against the median abundance for the same element
#'in the cohort (if \code{m} is provided) or against the
#'median abundance of \code{x} (when \code{m = NULL}).
#'If an element \eqn{x_{j}} is greater than the median, it
#'receives a score of \code{+1}; if it is lower than the median,
#'it receives a score of \code{-1}. The final score is the sum
#'of the scores for each element in \code{i}:
#'
#'\deqn{bristowScore = \sum_{j=1}^{n}  f(x) =
#'\left\{
#'\begin{array}{l}
#'+1, \quad \textrm{if} \quad x_{j} > median_{j} \\
#'-1, \quad \textrm{if} \quad x_{j} < median_{j}
#'\end{array}
#'\right\}
#'}
#'
#'where \eqn{n} is the number of elements in \code{i};
#'\eqn{x_{j}} is the \eqn{j}-th element of vector \code{x};
#'\eqn{median_{j}} is the \eqn{j}-th element of vector \code{m} if
#'provided, or the median of \code{x} when \code{m = NULL}.
#'
#'@inherit computeScore author
#'
#'@references
#'Lalonde, E. et al., *Tumour genomic and microenvironmental heterogeneity for integrated prediction of 5-year biochemical recurrence of prostate cancer: a retrospective cohort study*,
#'The Lancet Oncology (2014), appendix
bristowScore <- function(
    x,
    i,
    m     = NULL,
    na.rm = TRUE
){

  #check input
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = length(x))}
  if(isTRUE(missing(m) || is.null(m))) {
    #compute
    m = median(x = x, na.rm = na.rm)
  }
  if(isTRUE(length(m)==1)){m = rep(x = m, times = length(i))}

  if(isFALSE(length(m)==length(i))){
    stop("Error: 'm' must be of length 1 or of the same length of 'i'. Please, check your input.\n")
  }

  #subset input
  i = updateSig(x = x, i = i)
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(i))){   return(getDefaultNaValue())}
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #update w
  m = m[getIndex(i)]

  #check NAs
  if(isTRUE(na.rm)){
    keep = !is.na(x)
    #update
    x = x[keep]
    m = m[keep]
  }

  #compute
  ##out
  out = vector(mode = "numeric", length = length(x))
  ##update
  out[x > m] = 1
  out[x < m] = -1
  ##score
  out = sum(out, na.rm = na.rm)

  #return
  return(out)
}


#'Reviewd Bristow Score
#'
#'@description This function computes the *reviewed Bristow score*.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeScore
#'@inheritParams bristowScore
#'
#'@inherit computeScore return
#'
#'@details To compute the score, each element in \code{i} is
#'evaluated against the median abundance for the same element
#'in the cohort (if \code{m} is provided) or against the
#'median abundance of \code{x} (when \code{m = NULL}).
#'If an element \eqn{x_{j}} is greater than the median, it
#'receives a score of \code{+1}; if it is lower than the median,
#'it receives a score of \code{-1}. The final score is the sum
#'of the scores for each element in \code{i}.
#'
#'Summarising, the score is computed in 2 steps:
#'\enumerate{
#'  \item compute the Bristow score
#'  \item divide the score by the number of elements in the signature
#'}
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link{bristowScore}}
reviewedBristowScore <- function(
    x,
    i,
    m     = NULL,
    na.rm = TRUE
){

  #compute
  ##bristowScore
  out = bristowScore(x = x, i = i, m = m, na.rm = na.rm)
  #subset input
  x = subsetVector(x = x, i = i)
  #check NAs
  if(isTRUE(na.rm)){ x = x[!is.na(x)]}
  ##reviewed score
  out = out / length(x)

  #return
  return(out)
}




#'Single Sample Gene Set Enrichment Analysis (ssGSEA) Score
#'
#'@description This function computes the *ssGSEA* score.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@details This function always returns NA.
#'
#'@inherit computeScore author
ssGseaScore <- function(
    x,
    i,
    na.rm   = TRUE,
    ...
){

  #default
  out = getDefaultNaValue()

  #return
  return(out)
}

#'Gene Set Variation Analysis (GSVA) Score
#'
#'@description This function computes the *GSVA* score.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@details This function always returns NA.
#'
#'@inherit computeScore author
gsvaScore <- function(
    x,
    i,
    na.rm   = TRUE,
    ...
){

  #default
  out = getDefaultNaValue()

  #return
  return(out)
}


#'Pathway Level Analysis of Gene Expression (PLAGE) Score
#'
#'@description This function computes the *PLAGE* score.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@details This function always returns NA.
#'
#'@inherit computeScore author
plageScore <- function(
    x,
    i,
    na.rm   = TRUE,
    ...
){

  #default
  out = getDefaultNaValue()

  #return
  return(out)
}

#'Z Score
#'
#'@description This function computes the *z-score*.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@details This function always returns NA.
#'
#'@inherit computeScore author
zScore <- function(
    x,
    i,
    na.rm   = TRUE,
    ...
){

  #default
  out = getDefaultNaValue()

  #return
  return(out)
}

#'Sum of Step Function (SSF) Score
#'
#'@description This function computes the *Sum of Step Function score*.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeScore
#'@param thr threshold, numerical vector of length equivalent to \code{i},
#'containing the measure of central tendency for each element of
#'\code{i} in the cohort. If \code{threshold = NULL}, the measure
#'of central tendency of \code{x} is instead used
#'@param method character string, the measure of central tendency to use
#'@param scale logical, whether to scale the score by the number of
#'elements in \code{i}
#'@inherit computeScore return
#'
#'@details To compute the score, each element in \code{i} is
#'evaluated against the median abundance for the same element
#'in the cohort (if \code{m} is provided) or against the
#'median abundance of \code{x} (when \code{m = NULL}).
#'If an element \eqn{x_{j}} is greater than the median, it
#'receives a score of \code{+1}; if it is lower than the median,
#'it receives a score of \code{-1}. The final score is the sum
#'of the scores for each element in \code{i}:
#'
#'\deqn{SSF = \sum_{j=1}^{n}  f(x) =
#'\left\{
#'\begin{array}{l}
#'+1, \quad \textrm{if} \quad x_{j} > threshold_{j} \\
#'-1, \quad \textrm{if} \quad x_{j} < threshold_{j}
#'\end{array}
#'\right\}
#'}
#'
#'where \eqn{n} is the number of elements in \code{i};
#'\eqn{x_{j}} is the \eqn{j}-th element of vector \code{x};
#'\eqn{median_{j}} is the \eqn{j}-th element of vector \code{m} if
#'provided, or the median of \code{x} when \code{m = NULL}.
#'
#'
#'@inherit computeScore author
ssfScore <- function(
    x,
    i,
    w,
    thr    = NULL,
    method = c("median", "mean", "mode", "midrange", "trimean", "iqm"),
    scale  = TRUE,
    na.rm  = TRUE
  ){

  #check input
  if(isTRUE(missing(i) || is.null(i))) {i = seq_len(length.out = length(x))}
  if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(i))}
  if(isTRUE(missing(thr) || is.null(thr))) {
    #match
    method = match.arg(method)

    #compute
    thr = computeScore(x = x, na.rm = na.rm, score = method)
  }
  if(isTRUE(length(thr)==1)){m = rep(x = thr, times = length(i))}

  if(isFALSE(length(thr)==length(i))){
    stop("Error: 'thr' must be of length 1 or of the same length of 'i'. Please, check your input.\n")
  }

  #subset input
  i = updateSig(x = x, i = i)
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(i))){   return(getDefaultNaValue())}
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #update
  thr = thr[getIndex(i)]
  w   =   w[getIndex(i)]

  #check NAs
  if(isTRUE(na.rm)){
    keep = !is.na(x)
    #update
    x = x[keep]
    w = w[keep]
    thr = thr[keep]
  }

  #compute
  ##out
  out = vector(mode = "numeric", length = length(x))
  ##update
  out[x > thr] = 1
  out[x < thr] = -1
  ##score
  if(isTRUE(scale)){
    out = stats::weighted.mean(x = out, w = w)
  } else {
    #weights
    out = out * w
    #score
    out = sum(out, na.rm = na.rm)
  }

  #return
  return(out)
}




