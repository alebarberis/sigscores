#'@include utility-functions.R
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
#'considered for the computation of the score.
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'@param score character string indicating the summary score to compute
#'@param ... further arguments to \code{score}
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
    na.rm,
    score = c("mean", "trimmedMean", "weightedMean", "median",
              "mode", "midrange", "midhinge",
              "trimean", "briskow", "reviewedBriskow",
              "IQM", "weightedSum",
              "ssGSEA", "gsva", "plage", "zscore"),
    ...
  ){

  #match
  score = match.arg(score)

  #compute
  out = switch(
    score,
    "mean"            =              meanScore(x = x, i = i, na.rm = na.rm, ...),
    "trimmedMean"     =       trimmedMeanScore(x = x, i = i, na.rm = na.rm, ...),
    "weightedMean"    =      weightedMeanScore(x = x, i = i, na.rm = na.rm, ...),
    "median"          =            medianScore(x = x, i = i, na.rm = na.rm, ...),
    "mode"            =              modeScore(x = x, i = i, na.rm = na.rm, ...),
    "midrange"        =          midrangeScore(x = x, i = i, na.rm = na.rm, ...),
    "midhinge"        =          midhingeScore(x = x, i = i, na.rm = na.rm, ...),
    "trimean"         =           trimeanScore(x = x, i = i, na.rm = na.rm, ...),
    "briskow"         =           briskowScore(x = x, i = i, na.rm = na.rm, ...),
    "reviewedBriskow" =   reviewedBriskowScore(x = x, i = i, na.rm = na.rm, ...),
    "IQM"             = interquartileMeanScore(x = x, i = i, na.rm = na.rm, ...),
    "weightedSum"     =       weightedSumScore(x = x, i = i, na.rm = na.rm, ...),
    "ssGSEA"          =            ssGseaScore(x = x, i = i, na.rm = na.rm, ...),
    "gsva"            =              gsvaScore(x = x, i = i, na.rm = na.rm, ...),
    "plage"           =             plageScore(x = x, i = i, na.rm = na.rm, ...),
    "zscore"          =                 zScore(x = x, i = i, na.rm = na.rm, ...)
  )

  #return
  return(out)
}


#'Arithmetic Mean Scores
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso [base::mean()]
#'
#'
meanScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##get value
  out = mean(x = x, na.rm = na.rm, trim = 0)

  #return
  return(out)
}



#'Trimmed Arithmetic Mean Score
#'
#'@inheritParams computeScore
#'@inheritParams base::mean
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso [base::mean()]
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

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##get value
  out = mean(x = x, na.rm = na.rm, trim = trim)

  #return
  return(out)
}




#'Weighted Arithmetic Mean Score
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
#'@seealso [stats::weighted.mean]
weightedMeanScore <- function(
    x,
    i,
    # w = rep(x = 1, times = length(x)),
    w,
    na.rm = TRUE
){

  #check input
  if(missing(i) || is.null(i)) {i = seq_len(length.out = length(x))}
  if(missing(w) || is.null(w)) {w = rep(x = 1, times = length(i))}

  #subset input
  i = updateSig(x = x, i = i)
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #update w
  w = w[getIndex(i)]

  #compute
  out = stats::weighted.mean(x = x, w = w, na.rm = na.rm)

  #return
  return(out)
}





#'Median Scores
#'
#'@description The median is the value separating the lower half
#'of a set of measures from the higher half.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso [stats::median()]
#'
#'
medianScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##get value
  out = median(x = x, na.rm = na.rm)

  #return
  return(out)
}


#'Mode Score
#'
#'@description Compute the mode of a vector, i.e. the value
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

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #check input
  useNA = if(isTRUE(na.rm)){"no"}else{"ifany"}
  isNum = is.numeric(x)

  # #has name
  # hasNames = isTRUE(!is.null(names(x)))

  #compute
  ##get the frequency of elements
  out = table(x, useNA = useNA)
  ##get index of max freq
  index = which.max(out)
  ##get value
  out = names(out)[index]

  # update
  # if(isFALSE(hasNames)){names(out) = NULL}
  if(isTRUE(isNum)){out = as.numeric(out)}

  #return
  return(out)
}




#'Midrange Score
#'
#'@description This function computes the midrange score, i.e.
#'the average of the lowest and highest values in a set of data.
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@details The midrange score is a measure of central tendency like
#'the mean, median, and mode. However, it is more prone to bias than
#'these other measures because it relies solely upon the two most
#'extreme scores, which could potentially be outliers.
#'
#'@inherit computeScore author
#'
#'@references https://dictionary.apa.org/midrange
#'
#'@seealso [meanScore, medianScore, modeScore]
midrangeScore <- function(
  x,
  i,
  na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##max
  maxValue = max(x = x, na.rm = na.rm)
  ##min
  minValue = min(x = x, na.rm = na.rm)
  ##score
  out = (maxValue + minValue) / 2

  #return
  return(out)
}



#'Midhinge Score
#'
#'@description The midhinge of a set of values
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

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##quartile
  quartileValues = stats::quantile(x = x, probs = c(0.25, 0.75), na.rm = TRUE)
  q1 = quartileValues[1]
  q2 = quartileValues[2]
  ##score
  out = stats::setNames(object = (q1 + q2) / 2, nm = NULL)

  #check output


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

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##median
  medianValue = median(x = x, na.rm = na.rm)
  ##midhinge
  midhingeValue = midhingeScore(x = x, na.rm = na.rm)
  ##score
  out = (medianValue + midhingeValue) / 2

  #return
  return(out)
}



#'Briskow Score
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
briskowScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##median
  medianValue = median(x = x, na.rm = na.rm)
  ##out
  out = vector(mode = "numeric", length = length(x))
  ##update
  out[x > medianValue] = 1
  out[x < medianValue] = -1
  ##score
  out = sum(out, na.rm = na.rm)

  #return
  return(out)
}


#'Reviewd Briskow Score
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
reviewedBriskowScore <- function(
    x,
    i,
    na.rm = TRUE
){

  #subset input
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #compute
  ##median
  out = briskowScore(x = x, na.rm = na.rm)
  ##score
  out = out / length(x)

  #return
  return(out)
}




#'Interquartile Mean (IQM) Score
#'
#'@description The interquartile mean is a statistical measure
#'of central tendency based on the truncated mean of the
#'interquartile range. It is computed as:
#'
#'\deqn{IQM(x) = \frac{2}{n}\sum_{i=\frac{n}{4}+1}^{\frac{3n}{4}}  x_{i}}
#'
#'where \eqn{x_{i}} is the \eqn{i}-th element of the ordered vector.
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
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}


  # #check NAs
  # if(isTRUE(na.rm)){
  #   x = x[!is.na(x)]
  # }
  #rm NAs
  x = x[!is.na(x)]

  #length
  n = length(x)

  #check length
  if(isTRUE(n > 3 & !all(is.na(x)))){

    #sort
    x = sort(x = x, decreasing = F)

    #quartile size
    qsize = trunc(x = n / 4)

    # start
    # istart = qsize + 1

    #subset
    x = x[-((n - qsize + 1):n)]
    x = x[-(1:qsize)]

    #check data is divisible by 4
    isDivBy4 = isTRUE((n %% 4) == 0)
    if(isDivBy4){
      # #end
      # iend   = 3 * qsize
      # #subset
      # x = x[istart:iend]
      #compute score
      out = sum(x) * 2 / n
    } else {
      #weight
      w = 1 - ((n / 4) - qsize)
      #elem
      x1 = x[1]
      xN = x[length(x)]
      #update
      x = x[-c(1, length(x))]
      #compute score
      out = (sum(x) + w*x1 + w*xN) * (2 / n)
    }

  } else {
    out = getDefaultNaValue()
  }

  #return
  return(out)
}



#'Weighted Sum Score
#'
#'@inheritParams computeScore
#@inheritParams stats::weighted.mean
#'@param w numerical vector of weights the same length as
#'\code{i}
#'@param normalisation the normalisation method to perform before
#'computing the weighted sum
#'@param ... further arguments to \code{normalisation}
#'
#'
#'@inherit computeScore return
#'
#'@inherit computeScore author
#'
#'@seealso
#'\code{\link[stats]{weighted.mean}},
#'\code{\link{quantileNormalization}}
weightedSumScore <- function(
    x,
    i,
    # w = rep(x = 1, times = length(x)),
    w,
    na.rm = TRUE,
    normalisation = "quantile",
    ...
){

  #default
  if(missing(i) || is.null(i)) {i = seq_len(length.out = length(x))}
  if(missing(w) || is.null(w)) {w = rep(x = 1, times = length(i))}

  #subset input
  i = updateSig(x = x, i = i)
  x = subsetVector(x = x, i = i)

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  #update w
  w = w[getIndex(i)]

  #match arg
  normalisation = match.arg(normalisation)

  #-----------------------------------------------------------#
  # Check NAs
  if(isTRUE(na.rm)){
    #check attribute
    oldna = stats::na.action(object = x)
    #remove
    x = stats::na.omit(object = x)
    #get index of removed elements
    rmv = stats::na.action(object = x)
    #update w
    if(isTRUE(is.null(oldna) & !is.null(rmv))){
      w = w[-rmv]
    }
    #clean
    rm(oldna, rmv)
    #check
    if(isFALSE(length(w)==length(x))){
      stop("Error after removing NAs: length of 'w' not matching length of 'x'.\n")
    }
  }

  #-----------------------------------------------------------#
  #0) convert to 1-column matrix
  x = matrix(data = x, ncol = 1)

  #1) Normalise data
  x =   x = normaliseData(x = x, method = normalisation, na.rm = na.rm, ...)

  #compute
  out = apply(X = x, MARGIN = 2, FUN = function(xi, w, na.rm){
    xi = xi * w
    out = sum(xi, na.rm = na.rm)
    return(out)
  }, w = w, na.rm = na.rm)

  if(isTRUE(is.na(out))){out = getDefaultNaValue()}

  #-----------------------------------------------------------#
  #return
  return(out)

}


#'Single Sample Gene Set Enrichment Analysis (ssGSEA) Score
#'
#'@inheritParams computeScore
#'
#'@inherit computeScore return
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
#'@inheritParams computeScore
#'
#'@inherit computeScore return
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
#'@inheritParams computeScore
#'
#'@inherit computeScore return
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
#'@inheritParams computeScore
#'
#'@inherit computeScore return
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
