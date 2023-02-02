#'@include 0-utility-functions.R 3-stats-functions.R
NULL

#'Normalise Data
#'
#'@description This function is a dispatcher for the
#'selected normalisation method.
#'
#'@param x numerical matrix, features-by-samples
#'@param method string, the normalisation method.
#'Available options are:
#'\describe{
#'\item{quantile}{quantile normalisation}
#'}
#'@param ... further arguments to the normalisation method
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
#'}
#'@param ... further arguments to the transformation method
#'
#'@return A numerical matrix containing the transformed
#'values.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{stepFunctionTranformation}},
#'\code{\link{quantileNormalization}}
#'
#'@keywords internal
transformData <- function(
    x,
    f = c('none', 'stepFunction', 'quantile'),
    ...){

  #match arg
  f = match.arg(f)

  #Compute
  x = switch(
    f,
    'none'         = x,
    'stepFunction' = stepFunctionTranformation(x = x, ...),
    'quantile'     = quantileNormalization(x = x, ...)
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
#'}
#'
#'@return A data transformation function.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{stepFunctionTranformation}},
#'\code{\link{quantileNormalization}}
#'
#'@export
getDataTransformer <- function(
    f = c('stepFunction', 'quantile')
  ){
  #match arg
  f = match.arg(f)

  #Compute
  x = switch(
    f,
    'stepFunction' = stepFunctionTranformation,
    'quantile'     = quantileNormalization
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
