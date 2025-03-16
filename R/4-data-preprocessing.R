#'@include 0-utility-functions.R 3-stats-functions.R
NULL

#'Normalise Data
#'
#'@description This function is a dispatcher for the
#'selected normalisation method.
#'
#'@param x numerical matrix, features-by-samples.
#'@param method string, the normalisation method.
#'Available options are:
#'\describe{
#'\item{quantile}{quantile normalisation}
#'}
#'@param ... further arguments to the normalisation method.
#'
#'@return A numerical matrix containing the normalized
#'values.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{quantileNormalization}}
#'
#'@keywords internal
normaliseData <- function(
    x,
    method = c('quantile'),
    ...){

  #match arg
  method = match.arg(method)

  #Compute
  x = switch(
    method,
    'quantile' =  quantileNormalization(x = x, ...)
  )

  #return
  return(x)
}



#'Transform Data
#'
#'@description This function is a dispatcher for the
#'selected transformation method.
#'
#'@param x numerical matrix, features-by-samples
#'@param f string, the transformation method.
#'Available options are:
#'\describe{
#' \item{\code{'none'}}{\code{x} is returned}
#' \item{\code{'stepFunction'}}{step function}
#' \item{\code{'quantile'}}{quantile normalisation}
#' \item{\code{'z'}}{z-standardisation}
#'}
#'@param ... further arguments to the transformation method.
#'
#'@return A numerical matrix containing the transformed
#'values.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{stepFunctionTranformation}},
#'\code{\link{quantileNormalization}}
#'\code{\link{zStandardisation}}
#'
#'@keywords internal
transformData <- function(
    x,
    f = c('none', 'stepFunction', 'quantile', 'z'),
    ...){

  #match arg
  f = match.arg(f)

  #Compute
  x = switch(
    f,
    'none'         = x,
    'stepFunction' = stepFunctionTranformation(x = x, ...),
    'quantile'     = quantileNormalization(x = x, ...),
    'z'            = zStandardisation(x = x, ...)
  )

  #return
  return(x)
}


#'Get Data Transformation Function
#'
#'@description This function is a dispatcher for the
#'selected transformation method.
#'
#'@param f string, the transformation method.
#'Available options are:
#'\describe{
#' \item{\code{'stepFunction'}}{step function}
#' \item{\code{'quantile'}}{quantile normalisation}
#' \item{\code{'z'}}{z-standardisation}
#'}
#'
#'@return A data transformation function.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{stepFunctionTranformation}},
#'\code{\link{quantileNormalization}},
#'\code{\link{zStandardisation}}
#'
#'@export
getDataTransformer <- function(
    f = c('stepFunction', 'quantile', 'z')
  ){
  #match arg
  f = match.arg(f)

  #Compute
  x = switch(
    f,
    'stepFunction' = stepFunctionTranformation,
    'quantile'     = quantileNormalization,
    'z'            = zStandardisation
  )

  #return
  return(x)
}


#'Quantile Normalization
#'
#'@description This function normalises the input
#'data via quantile normalisation technique.
#'
#'@param x numerical matrix, features-by-samples
#'@inheritParams base::rank
#'
#'@inherit normaliseData return
#'
#'@inherit normaliseData author
#'
#'@references https://www.nature.com/articles/s41598-020-72664-6,
#'https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#'
#'@examples
#'\dontrun{
#'
#'x = as.data.frame(
#'matrix(
#'data = c(2,3,6,5,5,3,6,6,
#'4,3,5,3,3,5,4,5,
#'5,4,3,4,7,2,5,5,
#'3,5,4,4,5,6,3,4,
#'4,5,5,6,6,5,5,7),
#'nrow = 5,
#'byrow = TRUE))
#'
#'quantileNormalization(x)
#'
#'}
#'
#'@keywords internal
quantileNormalization <- function(
    x,
    na.rm       = TRUE,
    ties.method = c( "min", "max", "first", "last")
  ){

  if(isTRUE(is.vector(x))){
    out = quantileNormalizationForVector(x=x)
  } else if(isTRUE(is.matrix(x))){
    out = quantileNormalizationForMatrix(x=x, na.rm=na.rm, ties.method = ties.method)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}


#'Quantile Normalization
#'
#'@description This function normalises the input
#'data via quantile normalisation technique.
#'
#'@param x numerical matrix, features-by-samples
#'@inheritParams base::rank
#'
#'@inherit normaliseData return
#'
#'@inherit normaliseData author
#'
#'@references https://www.nature.com/articles/s41598-020-72664-6,
#'https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
#'
#'@examples
#'\dontrun{
#'
#'x = as.data.frame(
#'matrix(
#'data = c(2,3,6,5,5,3,6,6,
#'4,3,5,3,3,5,4,5,
#'5,4,3,4,7,2,5,5,
#'3,5,4,4,5,6,3,4,
#'4,5,5,6,6,5,5,7),
#'nrow = 5,
#'byrow = TRUE))
#'
#'quantileNormalization(x)
#'
#'}
#'@export
#'@keywords internal
quantileNormalizationForMatrix <- function(
    x,
    na.rm       = TRUE,
    ties.method = c( "min", "max", "first", "last")
){

  ties.method = match.arg(ties.method)

  #rm NAs
  if(isTRUE(na.rm)){
    x = stats::na.omit(object = x)
  }

  dimX = dim(x)

  #0) Store row and column names
  rnames = rownames(x)
  cnames = colnames(x)

  #1) Store Ranks
  xrank = as.matrix(as.data.frame(apply(X = x, MARGIN = 2, FUN = rank, ties.method = ties.method, simplify = F)))

  #2) Sort matrix
  x = as.matrix(as.data.frame(apply(X = x, MARGIN = 2, FUN = sort, decreasing = FALSE, simplify = F, na.last = TRUE)))

  #3) Compute the mean
  # rowM = rowMeans(x = x, na.rm = na.rm, dims = 1)
  rowM = matrixStats::rowMeans2(x = x, na.rm = na.rm, dim. = dimX)

  #4) Substitute
  x = apply(X = xrank, MARGIN = 2, FUN = function(v, rowM){
    return(rowM[v])
  }, rowM = rowM)

  #5) Create output
  x = matrix(data = x, nrow = dimX[1], ncol = dimX[2])
  rownames(x) = rnames
  colnames(x) = cnames

  #6) return
  return(x)
}


#'Quantile Normalization
#'
#'@description Currently this function just returns
#'\code{x} as it is.
#'
#'@param x a (named) numerical vector
#'
#'@return This function always returns \code{x}.
#'
#'@inherit normaliseData author
#'
#'@keywords internal
quantileNormalizationForVector <- function(x){return(x)}

# #'Class-Specific Quantile Normalization
# classSpecificQuantileNormalization <- function(
    #     x
# ){
#
# }
#
# #'Class-Specific Quantile Normalization
# discreteQuantileNormalization <- function(
    #     x
# ){
#
# }
#
# #'Smooth Quantile Normalisation
# smoothQuantileNormalisation <- function(
    #     x
# ){
#
# }





#'Step Function Transformation
#'
#'@description This function transform input data via
#'a *Step Function*.
#'See \code{\link{stepFunction}} for further information.
#'
#'@param x a (named) numerical vector or a matrix features-by-samples
#'@inheritParams stepFunction
#'@param method character string, indicating the score to be used as
#'threshold. This parameter is only used when \code{thr = NULL} or
#'not provided
#'@param by character string, indicating whether to compute the
#'threshold by applying \code{method} to rows or columns of
#'\code{x}.
#'In the first case, the threshold for the \eqn{i}-th row
#'is the same across columns.
#'In the second case, the threshold is the same within
#'the same column vector.
#'It is used when \code{x} is a matrix
#'@param na.rm logical, whether to remove \code{NA} values
#'before the computation of the threshold
#'
#'@return A numerical vector or matrix, containing the output
#'of the transformation.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'#set seed for reproducibility
#'set.seed(seed = 5381L)
#'
#'#x is a vector
#'x = sample(x = -10:10, size = 10)
#'stepFunctionTranformation(x, thr = 0)
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
#'stepFunctionTranformation(x, thr = 500)
#'
#'#use median as threshold
#'stepFunctionTranformation(x, method = 'median')
#'
#'#use median as threshold (by columns)
#'stepFunctionTranformation(x, method = 'median', by = 'cols')
#'
#'@export
stepFunctionTranformation <- function(
  x,
  y      = c(-1,0,1),
  thr    = NULL,
  method = c("median", "mean", "mode", "midrange", "trimean", "iqm", "iqr", "mad"),
  by     = c("rows", "cols"),
  na.rm  = TRUE
){

  #match
  method = match.arg(method)
  by     = match.arg(by)

  if(isTRUE(is.vector(x))){
    out = stepFunctionTranformationForVector(x=x,y=y,thr=thr,method=method,na.rm=na.rm)
  } else if(isTRUE(is.matrix(x))){
    out = stepFunctionTranformationForMatrix(x=x,y=y,thr=thr,method=method,by=by,na.rm=na.rm)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}

#'Step Function Transformation for Vector
#'
#'@inherit stepFunctionTranformation description
#'
#'@param x a (named) numerical vector
#'@inheritParams stepFunctionTranformation
#'
#'@return A vector with the transformed values.
#'
#'@inherit stepFunctionTranformation author
#'
#'@keywords internal
stepFunctionTranformationForVector <- function(
    x,
    y   = c(-1,0,1),
    thr,
    method = c("median", "mean", "mode", "midrange", "trimean", "iqm", "iqr", "mad"),
    na.rm  = T
){

  if(isTRUE(missing(thr) || is.null(thr))) {
    #match
    method = match.arg(method)

    #compute
    thr = computeMeasure(x = x, na.rm = na.rm, score = method)
  }

  #check
  if(isTRUE(length(thr)==1)){thr = rep(x = thr, times = length(x))}

  #check NAs
  # if(isTRUE(na.rm)){
  #   keep = !is.na(x)
  #   #update
  #   x   =   x[keep]
  #   thr = thr[keep]
  #   #clean
  #   rm(keep)
  # }

  #transform
  x = stepFunction(x = x, y = y, thr = thr)

  return(x)
}

#'Step Function Transformation for Matrix
#'
#'@inherit stepFunctionTranformation description
#'
#'@param x a matrix features-by-samples
#'@inheritParams stepFunctionTranformation
#'
#'@return A matrix with the transformed values.
#'
#'@keywords internal
stepFunctionTranformationForMatrix <- function(
    x,
    y   = c(-1,0,1),
    thr,
    method = c("median", "mean", "mode", "midrange", "trimean", "iqm", "iqr", "mad"),
    by     = c("rows", "cols"),
    na.rm = T
){

  #match
  method = match.arg(method)
  by     = match.arg(by)

  #check
  if(isTRUE(missing(thr) || is.null(thr))) {
    if(isTRUE(identical(by, "cols"))){
      #compute
      thr = apply(X = x, MARGIN = 2, FUN = computeMeasure, score = method, na.rm = na.rm, simplify = FALSE)
    } else {
      #compute
      thr = apply(X = x, MARGIN = 1, FUN = computeMeasure, score = method, na.rm = na.rm, simplify = FALSE)
    }
    #as vector
    thr = as.vector(unlist(thr))
  }

  #compute
  x = stepFunction(x = x, y = y, thr = thr, by = by)

  return(x)
}



#'Z-Standardisation
#'
#'@description This function transform the input
#'data via z-standardisation.
#'See the **Details** section below for further information.
#'
#'@param x a (named) numerical vector or a matrix features-by-samples.
#'@param by character string, indicating whether to compute the
#'standardisation by rows or columns of \code{x}.
#'It is used when \code{x} is a matrix.
#'@param robust logical, whether to compute a robust z-standardisation.
#'@param na.rm logical, whether to remove \code{NA} values
#'before the computation.
#'
#'@return A numerical vector or matrix, containing the output
#'of the transformation.
#'
#'@details The z-standardisation (or z-score normalisation) is a method for
#'transforming the data so that it has a mean of zero and a standard deviation
#'of one. It is computed as:
#'
#'\deqn{z(x) = \frac{x - \mu}{\sigma}}
#'
#'where \eqn{\mu} and \eqn{\sigma} are the mean and standard deviation of the
#'population, respectively.
#'
#'Since z-scores can be affected by unusually large or small data values,
#'we can also compute a more robust modified version as:
#'
#'\deqn{z(x) = \frac{x - \mathrm{median}}{\mathrm{MAD}}}
#'
#'where \eqn{\mathrm{MAD}} is the median absolute deviation of the population.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{zStandardisationForVector}},
#'\code{\link{zStandardisationForMatrix}}
#'
#'
#'@examples
#'#set seed for reproducibility
#'set.seed(seed = 5381L)
#'
#'#x is a vector
#'x = sample(x = -10:10, size = 10)
#'#Standard z-score transformation
#'zStandardisation(x)
#'#Robust z-score transformation
#'zStandardisation(x=x,robust=TRUE)
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
#'#Standardise the features
#'zStandardisation(x=x,by="rows")
#'
#'@export
zStandardisation <- function(
  x,
  robust = FALSE,
  by     = c("rows", "cols"),
  na.rm  = TRUE
){
  # Match
  by = match.arg(by)

  if(isTRUE(is.vector(x))){
    out = zStandardisationForVector(x=x,robust=robust,na.rm=na.rm)
  } else if(isTRUE(is.matrix(x))){
    out = zStandardisationForMatrix(x=x,robust=robust,by=by,na.rm=na.rm)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}

#'Z-Standardisation for Vector
#'
#'@description This function transform the input
#'data via z-standardisation.
#'See the **Details** section below for further information.
#'@param x a (named) numerical vector.
#'@inheritParams zStandardisation
#'
#'@return A vector with the transformed values.
#'
#'@details The z-standardisation (or z-score normalisation) is a method for
#'transforming the data so that it has a mean of zero and a standard deviation
#'of one. It is computed as:
#'
#'\deqn{z(x) = \frac{x - \mu}{\sigma}}
#'
#'where \eqn{\mu} and \eqn{\sigma} are the mean and standard deviation of the
#'population, respectively.
#'
#'Since z-scores can be affected by unusually large or small data values,
#'we can also compute a more robust modified version as:
#'
#'\deqn{z(x) = \frac{x - \mathrm{median}}{\mathrm{MAD}}}
#'
#'where \eqn{\mathrm{MAD}} is the median absolute deviation of the population.
#'
#'@inherit zStandardisation author
#'
#'@references https://en.wikipedia.org/wiki/Standard_score
#'
#'@keywords internal
zStandardisationForVector <- function(
    x,
    robust = FALSE,
    na.rm  = TRUE
){

  # Check
  if(isTRUE(robust)){
    # Center
    center = stats::median(x = x, na.rm = na.rm)
    # Scale
    scale = stats::mad(x = x, center = center, na.rm = na.rm)
  } else {
    # Center
    center = mean(x = x, trim = 0, na.rm = na.rm)
    # Scale
    scale = stats::sd(x = x, na.rm = na.rm)
  }
  # Compute
  x = (x - center) / scale

  return(x)

}

#'Z-Standardisation for Matrix
#'
#'@description This function transform the input
#'data via z-standardisation.
#'See the **Details** section below for further information.
#'
#'@param x numerical matrix, features-by-samples.
#'@param robust logical, whether to compute a robust z-standardisation.
#'@param by character string, indicating whether to compute the
#'standardisation by rows or columns of \code{x}.
#'@param na.rm logical, whether to remove \code{NA} values
#'before the computation.
#'
#'@inherit zStandardisationForVector details
#'
#'@inherit zStandardisation author
#'
#'@references https://en.wikipedia.org/wiki/Standard_score
#'
#'@keywords internal
zStandardisationForMatrix <- function(
    x,
    robust = FALSE,
    by     = c("rows", "cols"),
    na.rm  = TRUE
){
  # Match arg
  by = match.arg(by)

  # Compute
  if(identical(by, "rows")){
    if(isTRUE(robust)){
      # Center
      center = matrixStats::rowMedians(x = x, na.rm = na.rm)
      # Scale
      scale = matrixStats::rowMads(x = x, center = center, na.rm = na.rm)
      # Standardise
      x = (x - center) / scale
    } else {
      # Center
      center = matrixStats::rowMeans2(x = x, na.rm = na.rm)
      # Scale
      scale = matrixStats::rowSds(x = x, center = center, na.rm = na.rm)
    }
    # Standardise
    x = (x - center) / scale
  } else {
    if(isTRUE(robust)){
      # Center
      center = matrixStats::colMedians(x = x, na.rm = na.rm)
      # Scale
      scale = matrixStats::colMads(x = x, center = center, na.rm = na.rm)
      # Standardise
      x = (x - center) / scale
    } else {
      # Center
      center = matrixStats::colMeans2(x = x, na.rm = na.rm)
      # Scale
      scale = matrixStats::colSds(x = x, center = center, na.rm = na.rm)
    }
    # Standardise
    x = sweep(x = x, MARGIN = 2L, STATS = center, FUN = "-") / scale
  }

  return(x)
}
