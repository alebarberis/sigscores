#'@include 0-utility-functions.R
NULL

# Generic Function for Statistical Measures -------------------------------

#'Compute Measure
#'
#'@description This function computes statistical measure(s)
#'from an input vector or matrix.
#'It is an internal function, and it is not intended to be
#'called by users.
#'
#'@param x (named) numerical vector or a matrix features-by-samples
#'@param score character string indicating the statistical measure
#' to compute
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'@param ... further arguments to \code{score}
#'
#'@return A numerical vector representing the computed measure(s).
#'
#'@details If \code{x} is a vector, the output is a single value,
#'if \code{x} is a matrix, the output is a vector of length
#'\code{ncol(x)}.
#'
#'A default \code{NA} value is returned if the score can't be
#'computed, or if \code{i} values are not present in \code{x}.
#'If \code{x} is a matrix, the output is a vector of \code{NA}
#'values.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{computeVectorMeasure}},
#'\code{\link{computeColMeasures}}
#'
#'@keywords internal
computeMeasure <- function(
    x,
    score = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqm", "iqr", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore"),
    na.rm = TRUE,
    ...){

  if(isTRUE(is.vector(x))){
    out = computeVectorMeasure(x=x,score=score,na.rm=na.rm,...)
  } else if(isTRUE(is.matrix(x))){
    out = computeColMeasures(x=x,score=score,na.rm=na.rm,...)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}


# Statistical measures for vector ----------------------------------------

#'Compute Measure
#'
#'@description This function computes a statistical measure
#'from an input vector. It is an internal function, and
#'it is not intended to be called by users.
#'
#'@param x (named) numerical vector
#'@param i (optional) numerical vector giving the position in \code{x}
#'or character vector matching the names in \code{x}.
#'If \code{missing} or \code{i = NULL}, the entire \code{x} is
#'considered for the computation of the score
#'@inheritParams computeMeasure
#'
#'@return A numerical value representing the computed measure.
#'A default \code{NA} value is returned if the score can't be
#'computed.
#A length-one numerical vector.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{summation}},
#'\code{\link{weightedSum}},
#'\code{\link{arithmeticMean}},
#'\code{\link{trimmedMean}},
#'\code{\link{weightedMean}},
#'\code{\link{midpoint}},
#'\code{\link{modalValue}},
#'\code{\link{midrange}},
#'\code{\link{midhinge}},
#'\code{\link{trimean}},
#'\code{\link{iqr}},
#'\code{\link{IQM}},
#'\code{\link{MAD}},
#'\code{\link{AAD}},
#'\code{\link{ssgsea}},
#'\code{\link{gsva}},
#'\code{\link{plage}},
#'\code{\link{zscore}}
#'@keywords internal
computeVectorMeasure <- function(
  x,
  i     = NULL,
  score = c("sum", "weightedSum",
    "mean", "trimmedMean", "weightedMean",
    "median", "mode", "midrange", "midhinge",
    "trimean", "iqm", "iqr", "mad", "aad",
    "ssgsea", "gsva", "plage", "zscore"),
  na.rm = TRUE,
  ...){

  #check input
  if(isTRUE(is.null(x))){   return(getDefaultNaValue())}
  if(isTRUE(all(is.na(x)))){return(getDefaultNaValue())}

  score = match.arg(score)

  #compute
  out = switch(
    score,
    "sum"             =            summation(x = x, i = i, na.rm = na.rm, ...),
    "weightedSum"     =          weightedSum(x = x, i = i, na.rm = na.rm, ...),
    "mean"            =       arithmeticMean(x = x, i = i, na.rm = na.rm, ...),
    "trimmedMean"     =          trimmedMean(x = x, i = i, na.rm = na.rm, ...),
    "weightedMean"    =         weightedMean(x = x, i = i, na.rm = na.rm, ...),
    "median"          =             midpoint(x = x, i = i, na.rm = na.rm, ...),
    "mode"            =           modalValue(x = x, i = i, na.rm = na.rm, ...),
    "midrange"        =             midrange(x = x, i = i, na.rm = na.rm, ...),
    "midhinge"        =             midhinge(x = x, i = i, na.rm = na.rm, ...),
    "trimean"         =              trimean(x = x, i = i, na.rm = na.rm, ...),
    "iqr"             =                  iqr(x = x, i = i, na.rm = na.rm, ...),
    "iqm"             =                  IQM(x = x, i = i, na.rm = na.rm, ...),
    "mad"             =                  MAD(x = x, i = i, na.rm = na.rm, ...),
    "aad"             =                  AAD(x = x, i = i, na.rm = na.rm, ...),
    "ssgsea"          =               ssgsea(x = x, i = i, na.rm = na.rm, ...),
    "gsva"            =                 gsva(x = x, i = i, na.rm = na.rm, ...),
    "plage"           =                plage(x = x, i = i, na.rm = na.rm, ...),
    "zscore"          =               zscore(x = x, i = i, na.rm = na.rm, ...)
  )

  #check output
  if(isTRUE(is.na(out))){out = getDefaultNaValue()}

  #return
  return(out)
}


#'Sum
#'
#'@description This function computes the sum of all
#'the values in \code{x}.
#'
#'@inheritParams computeVectorMeasure
#'@param ... further arguments to \code{\link[base]{sum}}
#'
#'@details It is a wrapper to \code{\link[base]{sum}}
#'function.
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'summation(x = c(1,2))
#'
#'@export
summation <- function(
    x,
    i     = NULL,
    na.rm = T,
    ...){
  #update
  if(!is.null(i)){
    x = x[i]
    w = w[i]
  }
  #score
  out = sum(x = x, na.rm = na.rm, ...)

  return(out)
}


#'Weighted Sum
#'
#'@description This function computes the sum of the
#'elements of a vector (\code{x}) multiplied by the elements
#'of another vector (\code{w}).
#'
#'@param x (named) numerical vector
#'@param w numerical vector of weights the same length as
#'\code{x}
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'x = c(1,1,2,3,3,3,NA)
#'w = c(1,1,8,1,1,1,1)
#'weightedSum(x = x, w = w, na.rm = TRUE)
#'
#'@export
weightedSum <- function(
    x,
    w,
    i     = NULL,
    na.rm = T){
  #check
  if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(x))}
  #update
  if(!is.null(i)){
    x = x[i]
    w = w[i]
  }
  #weights
  out = x * w
  #score
  out = sum(out, na.rm = na.rm)

  return(out)
}



#'Arithmetic Mean
#'
#'@description This function computes the *arithmetic mean*.
#'
#'@inheritParams computeVectorMeasure
#'@param ... further arguments to \code{\link[base]{mean}}
#'
#'@details It is a wrapper to \code{\link[base]{mean}}
#'function.
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'x = c(1,1,2,3,3,NA)
#'arithmeticMean(x = x, na.rm = TRUE)
#'
#'@export
arithmeticMean <- function(
    x,
    i     = NULL,
    na.rm = T,
    ...){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  out = mean(x, na.rm = na.rm, trim = 0, ...)
  return(out)
}


#'Weighted Arithmetic Mean
#'
#'@description This function computes the *weighted arithmetic mean*.
#'
#'@inheritParams computeVectorMeasure
#'@param ... further arguments to \code{\link[stats]{weighted.mean}}
#'
#'@details It is a wrapper to \code{\link[stats]{weighted.mean}}
#'function.
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'x = c(1,1,2,3,3,3,NA)
#'w = c(1,1,8,1,1,1,1)
#'weightedMean(x = x, w = w, na.rm = TRUE)
#'
#'@export
weightedMean <- function(
  x,
  w,
  i     = NULL,
  na.rm = TRUE,
  ...){
  #check
  if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(x))}
  #update
  if(!is.null(i)){
    x = x[i]
    w = w[i]
  }
  #compute
  out = stats::weighted.mean(x = x, w = w, na.rm = na.rm, ...)
  return(out)
}


#'Trimmed Arithmetic Mean Score
#'
#'@description This function computes the *trimmed arithmetic mean*
#'using \code{\link[base]{mean}}.
#'
#'@inheritParams base::mean
#'
#'@inherit base::mean return
#'
#'@author Alessandro Barberis
#'
#'@details It is a wrapper to \code{\link[base]{mean}}
#'function.
#'
#'@seealso
#'\code{\link[base]{mean}}
#'
#'@examples
#'x = c(seq_len(length.out = 10), NA)
#'trimmedMean(x = x, na.rm = TRUE, trim = 0.1)
#'
#'@export
trimmedMean <- function(
  x,
  trim  = 0,
  i     = NULL,
  na.rm = TRUE,
  ...
  ){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  out = mean(x = x, na.rm = na.rm, trim = trim, ...)
  return(out)
}


#'Median
#'
#'@description This function computes the *median*.
#'
#'@inheritParams computeVectorMeasure
#'@param ... further arguments to \code{\link[stats]{median}}
#'
#'@details It is a wrapper to \code{\link[stats]{median}}
#'function.
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'x = c(seq_len(length.out = 10), NA)
#'midpoint(x = x, na.rm = TRUE)
#'
#'@export
midpoint <- function(
  x,
  i     = NULL,
  na.rm = T,
  ...){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  out = stats::median(x = x, na.rm = na.rm, ...)
  return(out)
}

#'Mode
#'
#'@description Compute the *mode* of a vector, i.e. the value
#' that has highest number of occurrences. If different values have
#' the same number of occurrences, the first one is reported.
#'
#'@inheritParams computeVectorMeasure
#'
#'@return A length-one numerical vector.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'x = c(1,1,2,3,3,3,NA)
#'modalValue(x = x, na.rm = TRUE)
#'
#'@export
modalValue <- function(
    x,
    i     = NULL,
    na.rm = TRUE){

  #update
  if(!is.null(i)){x = x[i]}

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

#'Midrange
#'
#'@description This function computes the *midrange* score, i.e.
#'the average of the lowest and highest values in a set of data.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details The midrange score is a measure of central tendency like
#'the mean, median, and mode. However, it is more prone to bias than
#'these other measures because it relies solely upon the two most
#'extreme scores, which could potentially be outliers.
#'
#'@inherit computeVectorMeasure author
#'
#'@references https://dictionary.apa.org/midrange
#'
#'@examples
#'x = c(0,1,1,2,3,3,3,100,NA)
#'midrange(x = x, na.rm = TRUE)
#'
#'@export
midrange <- function(
    x,
    i     = NULL,
    na.rm = T
){
  #update
  if(!is.null(i)){x = x[i]}

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


#'Midhinge
#'
#'@description The *midhinge* of a set of values
#'is the mean of the first and third quartiles.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeVectorMeasure
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@details This measure of location is computed as:
#'
#'\deqn{MH(x) = \frac{Q_{1}(x) + Q_{3}(x)}{2} = \frac{median(x) + midhinge(x)}{2}}
#'
#'It is related to the interquartile range (\eqn{IQR = Q_{3} - Q_{1}}),
#'which is a measure of dispersion.
#'
#'Since missing values and NaN's are allowed only
#'if \code{na.rm = TRUE} in \code{\link[stats]{quantile}},
#'if \code{x} has \code{NA} values and \code{na.rm = FALSE}
#'then \code{NA} is returned.
#'
#'@references https://en.wikipedia.org/wiki/Midhinge
#'
#'@examples
#'x = c(0,1,1,2,3,3,3,100,NA)
#'midhinge(x = x, na.rm = TRUE)
#'
#'@export
midhinge <- function(
    x,
    i     = NULL,
    na.rm = TRUE
){
  #update
  if(!is.null(i)){x = x[i]}

  #compute
  if(isTRUE(!na.rm && any(is.na(x)))){
    out = getDefaultNaValue()
  } else {
    ##quartile
    quartileValues = stats::quantile(x = x, probs = c(0.25, 0.75), na.rm = TRUE)
    q1 = quartileValues[1]
    q3 = quartileValues[2]
    ##score
    out = stats::setNames(object = (q1 + q3) / 2, nm = NULL)
  }

  #return
  return(out)
}

#'Tukey's Trimean
#'
#'@description The *trimean* contains information about the
#'center and some of the position of the data.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@details The trimean is defined as a weighted average of
#'the distribution's median and its two quartiles:
#'
#'\deqn{TM(x) = \frac{Q_{1}(x) + 2Q_{2}(x) + Q_{3}(x)}{4} = \frac{median(x) + midhinge(x)}{2}}
#'
#'
#'@references https://en.wikipedia.org/wiki/Trimean
#'
#'@examples
#'x = c(0,1,1,2,3,3,3,100,NA)
#'trimean(x = x, na.rm = TRUE)
#'
#'@export
trimean <- function(
  x,
  i     = NULL,
  na.rm = TRUE
){
  #update
  if(!is.null(i)){x = x[i]}
  #median
  medianValue = median(x = x, na.rm = na.rm)
  #midhinge
  midhingeValue = midhinge(x = x, na.rm = na.rm)
  #score
  out = (medianValue + midhingeValue) / 2
  return(out)
}



#'Interquartile Range
#'
#'@description This function computes the *interquartile range*
#'of the \code{x} values.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeVectorMeasure
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@details The *interquartile range* (also called *midspread*
#'or *H-spread*) is a measure of statistical dispersion.
#'It is defined as the difference between the 75th and 25th
#'percentiles of the data:
#'
#'\deqn{IQR(x) = Q_{3}(x) - Q_{1}(x)}
#'
#'Since missing values and NaN's are allowed only
#'if \code{na.rm = TRUE} in \code{\link[stats]{quantile}},
#'if \code{x} has \code{NA} values and \code{na.rm = FALSE}
#'then \code{NA} is returned.
#'
#'@references https://en.wikipedia.org/wiki/Interquartile_range
#'
#'@examples
#'x = c(7,7,31,31,47,75,87,115,116,119,119,155,177,NA)
#'iqr(x = x)#88
#'
#'@export
iqr <- function(
  x,
  i     = NULL,
  na.rm = TRUE){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  if(isTRUE(!na.rm && any(is.na(x)))){
    out = getDefaultNaValue()
  } else {
    ##quartile
    quartileValues = stats::quantile(x = x, probs = c(0.25, 0.75), na.rm = TRUE)
    q1 = quartileValues[1]
    q3 = quartileValues[2]
    ##score
    out = stats::setNames(object = (q3 - q1), nm = NULL)
  }
  #return
  return(out)
}



#'Interquartile Mean
#'
#'@description This function computes the *interquartile mean*
#'of the \code{x} values.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeVectorMeasure
#'@param na.rm unused argument, provided for consistency
#'with other functions.
#'\code{NA} values are always removed from \code{x} before
#'computation
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@details The *interquartile mean* is a statistical measure
#'of central tendency based on the truncated mean of the
#'interquartile range. It is computed as:
#'
#'\deqn{IQM(x) = \frac{2}{n}\sum_{i=\frac{n}{4}+1}^{\frac{3n}{4}}  x_{i}}
#'
#'where \eqn{x_{i}} is the \eqn{i}-th element of the ordered vector.
#'
#'Like the median it is insensitive to outliers.
#'
#'@references https://en.wikipedia.org/wiki/Interquartile_mean
#'
#'@examples
#'#Dataset size divisible by four
#'x = c(5,8,4,38,8,6,9,7,7,3,1,6)
#'IQM(x)#6.5
#'
#'#Dataset size not divisible by four
#'x = c(1,2,3,4,5)
#'IQM(x = x)#3
#'
#'x = c(1,3,5,7,9,11,13,15,17)
#'IQM(x)#9
#'
#'@export
IQM <- function(
    x,
    i     = NULL,
    na.rm = T){
  #update
  if(!is.null(i)){x = x[i]}
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


#'Median Absolute Deviation
#'
#'@description This function computes the *median absolute deviation*.
#'
#'@inheritParams computeVectorMeasure
#'@param ... further arguments to \code{\link[stats]{mad}}
#'
#'@details It is a wrapper to \code{\link[stats]{mad}}
#'function.
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@references https://en.wikipedia.org/wiki/Median_absolute_deviation
#'
#'@examples
#'x = c(1,1,2,2,4,6,9)
#'MAD(x)#1.4826
#'
#'@export
MAD <- function(
  x,
  i     = NULL,
  na.rm = TRUE,
  ...){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  out = stats::mad(x = x, na.rm = na.rm, ...)
  return(out)
}

#'Average Absolute Deviation
#'
#'@description This function computes the *average absolute deviation*,
#'i.e. the mean of the absolute deviations from a central point
#'(median by default).
#'
#'@inheritParams computeVectorMeasure
#'@param center numerical value, the central point
#'
#'@inherit computeVectorMeasure return
#'
#'@inherit computeVectorMeasure author
#'
#'@references https://en.wikipedia.org/wiki/Average_absolute_deviation
#'
#'@examples
#'x = c(2,2,3,4,14)
#'AAD(x = x, center = median(x))#2.8
#'AAD(x = x, center = mean(x))#3.6
#'AAD(x = x, center = modalValue(x))#3.0
#'
#'@export
AAD <- function(
  x,
  i     = NULL,
  na.rm  = TRUE,
  center = NULL
  ){
  #update
  if(!is.null(i)){x = x[i]}
  #compute
  if(isTRUE(missing(center) || is.null(center))){
    center = median(x = x, na.rm = na.rm)
  }
  out = abs(x = x - center)
  out = mean(out, na.rm = na.rm)
  return(out)
}

#'Enrichment Score
#'
#'@description This function computes an enrichment score.
#'
#'@param method character string, the name of the method
#'to be used for computing the enrichment scores.
#'Current available options are:
#'\describe{
#' \item{\code{'ssgsea'}}{Single Sample Gene Set Enrichment Analysis}
#' \item{\code{'gsva'}}{Gene Set Variation Analysis}
#' \item{\code{'plage'}}{Pathway Level Analysis of Gene Expression}
#' \item{\code{'zscore'}}{Z-Score}
#'}
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details This function always returns NA.
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'enrichment()
#'
#'@keywords internal
enrichment <- function(
  method = c("ssgsea", "gsva", "plage", "zscore"),
  ...
){
  out = getDefaultNaValue()
  return(out)
}

#'Enrichment Score
#'
#'@description This function computes an enrichment score.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details This function always returns NA.
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'ssgsea()
#'
#'@export
ssgsea <- function(...){return(enrichment(method = "ssgsea", ...))}

#'Enrichment Score
#'
#'@description This function computes an enrichment score.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details This function always returns NA.
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'gsva()
#'
#'@export
gsva   <- function(...){return(enrichment(method = "gsva", ...))}

#'Pathway Level Analysis of Gene Expression (PLAGE) Score
#'
#'@description This function computes an enrichment score.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details This function always returns NA.
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'plage()
#'
#'@export
plage  <- function(...){return(enrichment(method = "plage", ...))}

#'Enrichment Score
#'
#'@description This function computes an enrichment score.
#'
#'@inheritParams computeVectorMeasure
#'
#'@inherit computeVectorMeasure return
#'
#'@details This function always returns NA.
#'
#'@inherit computeVectorMeasure author
#'
#'@examples
#'zscore()
#'
#'@export
zscore <- function(...){return(enrichment(method = "zscore", ...))}

# Statistical measures for matrix ----------------------------------------

#'Computes Measure for Each Column in a Matrix
#'
#'@description This function calculates statistical measure(s)
#'for each column in an input matrix.
#'It is an internal function, and it is not intended to be
#'called by users.
#'
#'@param x a numerical matrix features-by-samples
#'@param rows (optional) numerical vector giving the rows
#'in \code{x} or character vector matching the row
#'names in \code{x} to operate over.
#'If \code{missing} or \code{rows = NULL}, all the rows
#'in \code{x} are considered for the computation of
#'the measures
#'@inheritParams computeMeasure
#'
#'@return A numerical vector containing the computed
#'score for each column in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@details
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
#'\code{...} parameter).
#'
#'@keywords internal
computeColMeasures <- function(
    x,
    rows  = NULL,
    score = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqm", "iqr", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore"),
    na.rm = TRUE,
    ...){

  #match
  score = match.arg(score)

  #get dim
  numCol = ncol(x)

  #names
  cnames = colnames(x)

  #check input
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = switch(
    score,
    "sum"             =               colSummations(x = x, na.rm = na.rm, rows = rows, ...),
    "weightedSum"     =             colWeightedSums(x = x, na.rm = na.rm, rows = rows, ...),
    "mean"            =          colArithmeticMeans(x = x, na.rm = na.rm, rows = rows, ...),
    "trimmedMean"     =             colTrimmedMeans(x = x, na.rm = na.rm, rows = rows, ...),
    "weightedMean"    =  colWeightedArithmeticMeans(x = x, na.rm = na.rm, rows = rows, ...),
    "median"          =                colMidpoints(x = x, na.rm = na.rm, rows = rows, ...),
    "mode"            =                    colModes(x = x, na.rm = na.rm, rows = rows, ...),
    "midrange"        =                colMidranges(x = x, na.rm = na.rm, rows = rows, ...),
    "midhinge"        =                colMidhinges(x = x, na.rm = na.rm, rows = rows, ...),
    "trimean"         =                 colTrimeans(x = x, na.rm = na.rm, rows = rows, ...),
    "iqr"             =                     colIQRs(x = x, na.rm = na.rm, rows = rows, ...),
    "iqm"             =                     colIQMs(x = x, na.rm = na.rm, rows = rows, ...),
    "mad"             =                     colMADs(x = x, na.rm = na.rm, rows = rows, ...),
    "aad"             =                     colAADs(x = x, na.rm = na.rm, rows = rows, ...),
    "ssgsea"          =                   colSsgsea(x = x, na.rm = na.rm, rows = rows, ...),
    "gsva"            =                     colGsva(x = x, na.rm = na.rm, rows = rows, ...),
    "plage"           =                    colPlage(x = x, na.rm = na.rm, rows = rows, ...),
    "zscore"          =                   colZscore(x = x, na.rm = na.rm, rows = rows, ...)
  )

  #update
  names(out) = cnames

  #return
  return(out)
}

#'Column Sums
#'
#'@description This function computes the *sum*
#'for each column in an input matrix.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colSums2}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[matrixStats]{colSums2}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,3,1,4), ncol = 2)
#'colSummations(x)
#'
#'@export
colSummations <- function(
  x,
  rows  = NULL,
  na.rm = TRUE,
  ...
){

  out = matrixStats::colSums2(
    x     = x,
    na.rm = na.rm,
    rows  = rows,
    ...
    # dims = 1L
    # rows  = i
  )

  return(out)
}


#'Column Weighted Sums
#'
#'@description This function computes the *weighted sum*
#'for each column in an input matrix.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... currently not used
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link[matrixStats]{colSums2}}
#'to each column of the matrix resulting from \code{w * x}.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,3,1,4), ncol = 2)
#'w = c(5, 1)
#'colWeightedSums(x = x, w = w)
#'
#'@export
colWeightedSums <- function(
  x,
  w,
  rows  = NULL,
  na.rm = TRUE,
  ...){

  #check
  if(isTRUE(missing(w) || is.null(w))) {w = rep(x = 1, times = length(x))}

  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
    w = w[rows]
  }

  #compute
  ##multiply
  out = w * x
  ##sum
  out = colSummations(
    x = out,
    na.rm = na.rm
  )

  #update
  out = unlist(out)

  return(out)
}

#'Column Means
#'
#'@description This function computes the *mean*
#'for each column in an input matrix.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colMeans2}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[matrixStats]{colMeans2}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,3,1,4), ncol = 2)
#'colArithmeticMeans(x = x)
#'
#'@export
colArithmeticMeans <- function(
    x,
    na.rm = TRUE,
    rows  = NULL,
    ...){
  out = matrixStats::colMeans2(
    x     = x,
    na.rm = na.rm,
    rows  = rows,
    ...
  )
  return(out)
}


#'Column Weighted Means
#'
#'@description This function computes the *weighted mean*
#'for each column in an input matrix.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colWeightedMeans}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[matrixStats]{colWeightedMeans}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,3,1,4), ncol = 2)
#'w = c(1,1)
#'colWeightedArithmeticMeans(x = x, w = w)
#'
#'@export
colWeightedArithmeticMeans <- function(
    x,
    w,
    rows  = NULL,
    na.rm = TRUE,
    ...){
  out = matrixStats::colWeightedMeans(
    x     = x,
    w     = w,
    na.rm = na.rm,
    rows  = rows,
    ...
  )
  return(out)
}


#'Column Trimmed Means
#'
#'@description This function computes the *trimmed mean*
#'for each column in an input matrix.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@inheritParams trimmedMean
#'@param ... further arguments to \code{\link{trimmedMean}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{trimmedMean}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = sample(seq_len(10)), ncol = 2)
#'colTrimmedMeans(x = x, trim = 0)
#'colTrimmedMeans(x = x, trim = 0.2)
#'
#'@export
colTrimmedMeans <- function(
  x,
  trim  = 0,
  rows  = NULL,
  na.rm = TRUE,
  ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimmedMean, trim = trim, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}

#'Column Medians
#'
#'@description This function computes the *median*
#'for each column in an input matrix.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colMedians}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[matrixStats]{colMedians}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,2,3,1,2,3), ncol = 2)
#'colMidpoints(x = x)
#'
#'@export
colMidpoints <- function(
    x,
    rows  = NULL,
    na.rm = TRUE,
    ...){
  #compute
  out = matrixStats::colMedians(
    x     = x,
    na.rm = na.rm,
    rows  = rows,
    ...
    # dims = 1L
    # rows  = i
  )

  return(out)
}


#'Column Modes
#'
#'@description This function computes the *mode*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link{modalValue}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{modalValue}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,2,2,1,1,3), ncol = 2)
#'colModes(x = x)
#'
#'@export
colModes <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = modalValue, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}


#'Column Midranges
#'
#'@description This function computes the *midrange*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... currently not used
#'
#'@inherit computeColMeasures return
#'
#'@details The midrange score is the average of the lowest
#'and highest values in a set of data.
#'It is a measure of central tendency like
#'the mean, median, and mode. However, it is more prone to bias than
#'these other measures because it relies solely upon the two most
#'extreme scores, which could potentially be outliers.
#'
#'This function calculates the range of values in each
#'column of \code{x} by calling \code{\link[matrixStats]{colRanges}}.
#'Then, it uses the min and max values to compute the score.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,2,2,1,1,3), ncol = 2)
#'colMidranges(x = x)
#'
#'@export
colMidranges <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  ##ranges
  ranges = matrixStats::colRanges(
    x     = x,
    na.rm = na.rm)
  ##max
  maxValues = ranges[,2]
  ##min
  minValues = ranges[,1]
  ##score
  out = (maxValues + minValues) / 2

  return(out)
}

#'Column Trimeans
#'
#'@description This function computes the *trimean*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link{trimean}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{trimean}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,2,2,1,1,3), ncol = 2)
#'colTrimeans(x = x)
#'
#'@export
colTrimeans <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimean, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}

#'Column Midhinges
#'
#'@description This function computes the *midhinge*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link{midhinge}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{midhinge}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = matrix(data = c(1,2,2,1,1,3), ncol = 2)
#'colMidhinges(x = x)
#'
#'@export
colMidhinges <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = midhinge, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}


#'Column Interquartile Ranges
#'
#'@description This function computes the *interquartile range*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link{iqr}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{iqr}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = c(7,7,31,31,47,75,87,115,116,119,119,155,177,NA)
#'x = cbind(x,x)
#'colnames(x) = c("S1", "S2")
#'colIQRs(x = x)#88 88
#'
#'@export
colIQRs <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){
  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = iqr, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}

#'Column Interquartile Mean
#'
#'@description This function computes the *interquartile mean*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link{IQM}}
#'
#'@inherit computeColMeasures return
#'
#'@details This function applies \code{\link{IQM}}
#'to each column of the input matrix.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = c(1,2,3,4,5)
#'x = cbind(x,x)
#'colnames(x) = c("S1", "S2")
#'colIQMs(x = x)#3 3
#'
#'@export
colIQMs <- function(
  x,
  rows  = NULL,
  na.rm = T,
  ...
  ){

  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = IQM, na.rm = na.rm, ..., simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}

#'Column Median Absolute Deviation
#'
#'@description This function computes the *median absolute deviation*
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colMads}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[matrixStats]{colMads}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = c(1,1,2,2,4,6,9)
#'x = cbind(S1 = x, S2 = x)
#'colMADs(x)#1.4826 1.4826
#'
#'@export
colMADs <- function(
    x,
    rows  = NULL,
    na.rm = T,
    ...
){

  #compute
  out = matrixStats::colMads(
    x     = x,
    na.rm = na.rm,
    rows  = rows,
    ...
  )

  return(out)
}

#'Column Average Absolute Deviation
#'
#'@description This function computes the *average absolute deviation*
#'for each column vector in the input matrix \code{x}.
#'
#'@inheritParams computeColMeasures
#'@param ... further arguments to \code{\link[matrixStats]{colMeans2}}
#'
#'@inherit computeColMeasures return
#'
#'@inherit computeColMeasures author
#'
#'@examples
#'x = cbind(c(2,2,3,4,14), c(2,2,3,4,14))
#'colAADs(x)#2.8 2.8
#'
#'@export
colAADs <- function(
    x,
    center,
    rows  = NULL,
    na.rm = T,
    ...
){

  #update
  if(isTRUE(!is.null(rows))){
    x = x[rows,,drop=F]
  }

  #check
  if(isTRUE(missing(center) || is.null(center))) {
    center = colMidpoints(x = x, na.rm = na.rm)
  }
  if(isTRUE(length(center)==1)){center = rep(x = center, times = ncol(x))}

  #subtract
  out = sweep(x = x, MARGIN = 2, STATS = center, FUN = "-")

  #abs
  out = abs(out)

  #compute
  out = matrixStats::colMeans2(
    x     = out,
    na.rm = na.rm,
    ...
  )
  return(out)
}

#'Enrichment Scores
#'
#'@description This function computes an enrichment score
#'for each column vector in the input matrix \code{x}.
#'See the **Details** section below for further information.
#'
#'@inheritParams computeColMeasures
#'@param method character string, the name of the method
#'to be used for computing the enrichment scores.
#'Current available options are:
#'\describe{
#' \item{\code{'ssgsea'}}{Single Sample Gene Set Enrichment Analysis}
#' \item{\code{'gsva'}}{Gene Set Variation Analysis}
#' \item{\code{'plage'}}{Pathway Level Analysis of Gene Expression}
#' \item{\code{'zscore'}}{Z-Score}
#'}
#'
#'@param ... further arguments to \code{\link[GSVA]{gsva}}
#'
#'@inherit computeColMeasures return
#'
#'@details It is a wrapper to \code{\link[GSVA]{gsva}}
#'function.
#'
#'@inherit computeColMeasures author
#'
#'@keywords internal
colEnrichment <- function(
  x,
  rows   = NULL,
  na.rm  = TRUE,
  method = c("ssgsea", "gsva", "plage", "zscore"),
  verbose = FALSE,
  ...){

  #match
  method = match.arg(method)

  #remove
  if(isTRUE(na.rm)){
    x = stats::na.omit(object = x)
    #get index of removed elements
    rmv = stats::na.action(object = x)
    #update
    if(!is.null(rmv)){
      rows = rows[-rmv]
    }
    #clean
    rm(rmv)
  }

  #Remove genes with constant expression values
  keep = !areRowValuesEqual(x)
  if(isFALSE(all(keep))){
    x = x[keep,,drop=F]
    rows = rows[keep]
  }

  #get dim
  numRow = nrow(x)
  numCol = ncol(x)

  #check
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(numRow < 2)){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(rownames(x))) & is.character(rows)){stop("Error: 'x' must have row names.\n")}

  #check rownames
  if(isTRUE(is.null(rownames(x)))){
    rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")
  }

  #check i
  if(isTRUE(is.numeric(rows))){
    # #check dim
    # rows = rows[rows<=numRow]
    #get names
    rows = rownames(x)[rows]
    #bool
    hasElements = length(rows) > 0
  } else {
    hasElements = length(intersect(x = rows, y = rownames(x))) > 0
  }

  #compute
  ##Check if has matching elements to avoid Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
  if(hasElements){
    out = unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = rows), verbose=verbose, method = method, ...))[,1])
  } else {
    out = rep(x = getDefaultNaValue(), times = numCol)
  }

  #return
  return(out)
}

#'Single Sample Gene Set Enrichment Analysis (ssGSEA) Scores
#'
#'@inherit colEnrichment description
#'@inheritParams colEnrichment
#'@inherit colEnrichment return
#'@inherit colEnrichment details
#'@inherit colEnrichment author
#'
#'@export
colSsgsea <- function(
    x,
    rows   = NULL,
    na.rm  = TRUE,
    verbose = FALSE,
    ...){
  return(colEnrichment(method = "ssgsea", x = x, rows = rows, na.rm = na.rm, verbose = verbose, ...))
}

#'Gene Set Variation Analysis (GSVA) Scores
#'
#'@inherit colEnrichment description
#'@inheritParams colEnrichment
#'@inherit colEnrichment return
#'@inherit colEnrichment details
#'@inherit colEnrichment author
#'
#'@export
colGsva   <- function(
    x,
    rows   = NULL,
    na.rm  = TRUE,
    verbose = FALSE,
    ...){
  return(colEnrichment(method = "gsva"  , x = x, rows = rows, na.rm = na.rm, verbose = verbose, ...))
}

#'Pathway Level Analysis of Gene Expression (PLAGE) Scores
#'
#'@inherit colEnrichment description
#'@inheritParams colEnrichment
#'@inherit colEnrichment return
#'@inherit colEnrichment details
#'@inherit colEnrichment author
#'
#'@export
colPlage  <- function(
    x,
    rows   = NULL,
    na.rm  = TRUE,
    verbose = FALSE,
    ...){
  return(colEnrichment(method = "plage" , x = x, rows = rows, na.rm = na.rm, verbose = verbose, ...))
}

#'Z-Scores
#'
#'@inherit colEnrichment description
#'@inheritParams colEnrichment
#'@inherit colEnrichment return
#'@inherit colEnrichment details
#'@inherit colEnrichment author
#'
#'@export
colZscore <- function(
    x,
    rows   = NULL,
    na.rm  = TRUE,
    verbose = FALSE,
    ...){
  return(colEnrichment(method = "zscore", x = x, rows = rows, na.rm = na.rm, verbose = verbose, ...))
}

# Mathematical functions --------------------------------------------------

#'Step Function
#'
#'@description This function implements a *Step Function*.
#'See the **Details** section below for further information.
#'
#'@param x a (named) numerical vector or a matrix features-by-samples
#'@param thr numerical value, the threshold to use to
#'divide values in \code{x} in two intervals
#'@param y numerical vector of length 3, the values to
#'use for the piecewise function.
#'The first value is used for the left interval,
#'the second value is used when \code{x == threshold},
#'and the third value is used for the right interval
#'@param by character string, indicating whether to
#'apply the threshold to rows or columns of
#'\code{x}.
#'In the first case, the threshold for the \eqn{i}-th row
#'is the same across columns.
#'In the second case, the threshold is the same within
#'the same column vector.
#'It is used when \code{x} is a matrix
#'
#'@return A numerical vector or matrix, containing the output
#'of the step function for each element of \code{x}.
#'The default settings of \code{thr} and \code{y}
#'makes this function equivalent to a sign function.
#'
#'@details This function implements a piecewise function
#'that is defined as follows:
#'
#'\deqn{SF = f(x) =
#'\left\{
#'\begin{array}{l}
#'y_{1}, \quad \textrm{if} \quad x > threshold \\
#'y_{2}, \quad \textrm{if} \quad x = threshold \\
#'y_{3}, \quad \textrm{if} \quad x < threshold
#'\end{array}
#'\right\}
#'}
#'
#'@author Alessandro Barberis
#'
#'@examples
#'#set seed for reproducibility
#'set.seed(seed = 5381L)
#'
#'#x is a vector
#'x = sample(x = -10:10, size = 10)
#'stepFunction(x, thr = 0)
#'
#'#x is a matrix
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
#'#compute
#'stepFunction(x, thr = 500)
#'
#'@keywords internal
stepFunction <- function(
    x,
    y   = c(-1,0,1),
    thr = 0,
    by = c("rows", "cols")
){

  if(isTRUE(is.vector(x))){
    out = stepFunctionForVector(x=x, y=y, thr=thr)
  } else if(isTRUE(is.matrix(x))){
    out = stepFunctionForMatrix(x=x, y=y, thr=thr, by = by)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)

}

#'Step Function
#'
#'@inherit stepFunction  description
#'
#'@param x (named) numerical vector
#'@param thr numerical value, the threshold to use to
#'divide values in \code{x} in two intervals
#'@inheritParams stepFunction
#'@return A numerical vector, containing the output
#'of the step function for each element of \code{x}.
#'The default settings of \code{thr} and \code{y}
#'makes this function equivalent to a sign function.
#'
#'@inherit stepFunction details
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
stepFunctionForVector <- function(
    x,
    y   = c(-1,0,1),
    thr = 0){

  #out
  out = rep(x = NA, times = length(x))

  #update
  out[x < thr]  = y[1]
  out[x == thr] = y[2]
  out[x > thr]  = y[3]

  #names
  names(out) = names(x)

  #return
  return(out)
}


#'Step Function
#'
#'@inherit stepFunction description
#'
#'@param x numerical matrix, features-by-samples
#'@param thr numerical vector, the threshold(s) to use to
#'divide values in \code{x} in two intervals.
#'If it is a single value, the same threshold is applied
#'throughout \code{x}.
#'If \code{by = 'rows'} and \code{thr} has number
#'of elements equivalent to rows in \code{x} then
#'a different threshold is applied to each row.
#'If \code{by = 'cols'} and \code{thr} has number
#'of elements equivalent to columns in \code{x} then
#'a different threshold is applied to each column
#'@inheritParams stepFunction
#'@param by character string, indicating whether to
#'apply the threshold to rows or columns of
#'\code{x}.
#'In the first case, the threshold for the \eqn{i}-th row
#'is the same across columns.
#'In the second case, the threshold is the same within
#'the same column vector
#' @return A numerical matrix, containing the output
#' of the step function for each element of \code{x}.
#' The default settings of \code{thr} and \code{y}
#' makes this function equivalent to a sign function.
#'
#'@inherit stepFunction details
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
stepFunctionForMatrix <- function(
    x,
    y   = c(-1,0,1),
    thr = 0,
    by = c("rows", "cols")){

  #0) Match
  by = match.arg(by)

  #1) dim
  dimX = dim(x)
  len  = length(thr)

  #2) Check input
  if(isTRUE(len>1)){
    if(isTRUE(identical(by, "rows") && len!=dimX[1])){stop("Error: 'thr' must be a single value or a vector of values of length nrow(x). Please, check your input.\n")}
    if(isTRUE(identical(by, "cols") && len!=dimX[2])){stop("Error: 'thr' must be a single value or a vector of values of length ncol(x). Please, check your input.\n")}
  }

  #1) Store row and column names
  rnames = rownames(x)
  cnames = colnames(x)

  #2) Apply function
  if(isTRUE(identical(by, "rows"))){
    x = as.matrix(as.data.frame(apply(X = x, MARGIN = 2, FUN = stepFunctionForVector, y = y, thr = thr, simplify = F)))
  } else {
    x = as.matrix(as.data.frame(apply(X = x, MARGIN = 1, FUN = stepFunctionForVector, y = y, thr = thr, simplify = F)))
    x = t(x)
  }


  #3) Create output
  x = matrix(data = x, nrow = dimX[1], ncol = dimX[2])
  rownames(x) = rnames
  colnames(x) = cnames

  #return
  return(x)
}


