#'@include utility-functions.R
NULL

#'Compute Scores
#'
#'@description TThis function computes summary
#'score(s) of the signature \code{i} in input
#'considering each column vector in the input matrix
#'\code{x}.
#'
#'@param x numerical matrix, features-by-samples
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
#@param rm.dupli logical, whether to remove
#duplicated elements in \code{i}
#'@param scores character vector, indicating the
#'summary score(s) to compute
#'@param sample.id logical, whether to report the
#'sample ID as a column in the output data frame
#'@param args named list, where the names must match the
#'\code{scores}. Each element in the list is another list
#'containing the arguments to pass to the function used
#'for computing the named score. For example,
#'\code{args = list(trimmedMean = list(trim = 0.4))}
#'indicates to use \code{trim = 0.4} when computing the
#'trimmed mean scores
#'
#'@inheritParams forLoop
#'
#'@return A numerical vector containing the computed
#'score for each sample.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeMatrixScores <- function(
    x,
    i,
    # w,
    na.rm    = TRUE,
    # rm.dupli = TRUE,
    scores = c("mean", "trimmedMean", "weightedMean", "median",
               "mode", "midrange", "midhinge",
               "trimean", "bristow", "reviewedBristow",
               "IQM", "weightedSum",
               "ssGSEA", "gsva", "plage", "zscore"),
    sample.id = FALSE,
    args,
    cores     = 1L
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
  ##i
  if(isTRUE(missing(i) | is.null(i))){
    i = seq_len(length.out = nrow(x))
  }
  ##scores
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
  ##args
  if(isFALSE(missing(args) | is.null(args))){
    if(isTRUE(is.list(args) && !is.null(names(args)) && (any(names(args) %in% scores)))){
      #update
      args = args[(names(args) %in% scores)]
    } else{
      stop("Error: provided 'args' not valid. It must be a named list, with names matching the provided 'scores'.\n")
    }
  }

  #compute scores ------------------------------------------------
  ##loop over scores
  ###number of iteration
  n.iter = length(scores)
  ###compute scores
  out = forLoop(
    n.iter   = n.iter,
    cores    = cores,
    .inorder = T,
    fun = function(iloop, scores, args, ...) {
      #compute
      o = do.call(
        what = computeMatrixScore,
        args = c(
          list(score = scores[iloop]),
          list(...),
          args[[scores[iloop]]]
        )
      )
      #return
      return(o)
    },
    x      = x,
    na.rm  = na.rm,
    scores = scores,
    i      = i,
    # w      = w,
    args   = args
  )

  ##set names
  names(out) = scores

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

#'Compute Score
#'
#'@description This function computes one summary
#'score from each column vector in the input matrix.
#'
#'@param x numerical matrix, features-by-samples
#'@param i (optional) numerical vector giving the rows
#'in \code{x} or character vector matching the row
#'names in \code{x}.
#'If \code{missing} or \code{i = NULL}, all the rows
#'in \code{x} are considered for the computation of
#'the score
#'@param w numerical vector of weights the same length as
#'\code{i}
#'@param na.rm logical, whether to remove \code{NA}
#'values before computation
#@param rm.dupli logical, whether to remove
#duplicated elements in \code{i}
#'@param score character string indicating the summary score to compute
#'@param ... further arguments to the function used
#'for computing the scores
#'
#'@return A numerical vector containing the computed
#'score for each sample.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeMatrixScore <- function(
    x,
    i,
    w,
    na.rm,
    # rm.dupli = TRUE,
    score = c("mean", "trimmedMean", "weightedMean", "median",
              "mode", "midrange", "midhinge",
              "trimean", "bristow", "reviewedBristow",
              "IQM", "weightedSum",
              "ssGSEA", "gsva", "plage", "zscore"),
    ...
){

  #match
  score = match.arg(score)

  #compute
  out = switch(
    score,
    "mean"            =              meanScores(x = x, i = i, na.rm = na.rm),
    "trimmedMean"     =       trimmedMeanScores(x = x, i = i, na.rm = na.rm, ...),
    "weightedMean"    =      weightedMeanScores(x = x, i = i, na.rm = na.rm, w = w),
    "median"          =            medianScores(x = x, i = i, na.rm = na.rm),
    "mode"            =              modeScores(x = x, i = i, na.rm = na.rm),
    "midrange"        =          midrangeScores(x = x, i = i, na.rm = na.rm),
    "midhinge"        =          midhingeScores(x = x, i = i, na.rm = na.rm),
    "trimean"         =           trimeanScores(x = x, i = i, na.rm = na.rm),
    "bristow"         =           bristowScores(x = x, i = i, na.rm = na.rm),
    "reviewedBristow" =   reviewedBristowScores(x = x, i = i, na.rm = na.rm),
    "IQM"             = interquartileMeanScores(x = x, i = i, na.rm = na.rm),
    "weightedSum"     =       weightedSumScores(x = x, i = i, na.rm = na.rm, w = w, ...),
    "ssGSEA"          =            ssGseaScores(x = x, i = i, na.rm = na.rm, ...),
    "gsva"            =              gsvaScores(x = x, i = i, na.rm = na.rm, ...),
    "plage"           =             plageScores(x = x, i = i, na.rm = na.rm, ...),
    "zscore"          =                 zScores(x = x, i = i, na.rm = na.rm, ...)
  )

  #return
  return(out)
}



# computeMatrixScore <- function(
#     x,
#     i     = NULL,
#     na.rm = TRUE,
#     # rm.dupli = TRUE,
#     score = c("mean", "trimmedMean", "weightedMean", "median",
#               "mode", "midrange", "midhinge",
#               "trimean", "bristow", "reviewedBristow",
#               "IQM", "weightedSum",
#               "ssGSEA", "gsva", "plage", "zscore"),
#     ...
# ){
#
#   #match
#   score = match.arg(score)
#
#   #compute
#   out = switch(
#     score,
#     "mean"            =              meanScores,
#     "trimmedMean"     =       trimmedMeanScores,
#     "weightedMean"    =      weightedMeanScores,
#     "median"          =            medianScores,
#     "mode"            =              modeScores,
#     "midrange"        =          midrangeScores,
#     "midhinge"        =          midhingeScores,
#     "trimean"         =           trimeanScores,
#     "bristow"         =           bristowScores,
#     "reviewedBristow" =   reviewedBristowScores,
#     "IQM"             = interquartileMeanScores,
#     "weightedSum"     =       weightedSumScores,
#     "ssGSEA"          =            ssGseaScores,
#     "gsva"            =              gsvaScores,
#     "plage"           =             plageScores,
#     "zscore"          =                 zScores
#   )
#
#   #args
#   args = list(x = x, i = i, na.rm = na.rm)
#
#   #get args
#   args = getArgsInFunFormals(fun = f, args = list(x = x, i = i, na.rm = na.rm), ...)
#
#   #compute
#   out = do.call(what = f, args = args)
#
#   #return
#   return(out)
# }





#'Arithmetic Mean Scores
#'
#'@description This function computes the arithmetic mean
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'@param ... further arguments to \code{\link[matrixStats]{colMeans2}}
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link[matrixStats]{colMeans2}}
meanScores <- function(
    x,
    i        = NULL,
    na.rm    = TRUE,
    ...
){

  #get dim
  numCol = ncol(x)
  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)
  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = matrixStats::colMeans2(
    x     = x,
    na.rm = na.rm,
    ...
    # dims = 1L
    # rows  = i
  )

  #return
  return(out)
}


#'Trimmed Arithmetic Mean Scores
#'
#'@description This function computes the trimmed arithmetic mean
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'@inheritParams base::mean
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{trimmedMeanScore}},
#'\code{\link[base]{mean}}
trimmedMeanScores <- function(
    x,
    i,
    trim  = 0,
    na.rm = TRUE
){

  #get dim
  numCol = ncol(x)
  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)
  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimmedMeanScore, trim = trim, na.rm = na.rm, simplify = F)

  #update
  out = unlist(out)

  #return
  return(out)
}

#'Weighted Arithmetic Mean Scores
#'
#'@description This function computes the weighted arithmetic mean
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#@inheritParams stats::weighted.mean
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{weightedMeanScore}},
#'\code{\link[stats]{weighted.mean}}
weightedMeanScores <- function(
    x,
    i,
    # w = rep(x = 1, times = nrow(x)),
    w,
    na.rm = TRUE
){

  #check weights
  if(missing(i) || is.null(i)) {i = seq_len(length.out = nrow(x))}
  if(missing(w) || is.null(w)) {w = rep(x = 1, times = length(i))}

  #get dim
  numCol = ncol(x)

  #subset input
  i = updateSig(x = x, i = i, rm.dupli = T)
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(i))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #update w
  w = w[getIndex(i)]

  #-----------------------------------------------------------#
  #compute
  out = apply(X = x, MARGIN = 2, FUN = weightedMeanScore, w = w, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #-----------------------------------------------------------#
  #return
  return(out)
}


#'Median Scores
#'
#'@description This function computes the sample median
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{medianScore}}
medianScores <- function(
    x,
    i,
    na.rm = TRUE
){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = medianScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}

#'Mode Scores
#'
#'@description This function computes the mode
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{modeScore}}
modeScores <- function(
    x,
    i,
    na.rm    = TRUE
){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = modeScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}





#'Midrange Scores
#'
#'@description This function computes the midrange score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{midrangeScore}}
midrangeScores <- function(
    x,
    i,
    na.rm    = TRUE
  ){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = midrangeScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}



#'Midhinge Scores
#'
#'@description This function computes the midhinge score
#' from each column vector in the input matrix.
#' The midhinge of a set of values is the mean of the first
#' and third quartiles.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{midhingeScore}}
midhingeScores <- function(
    x,
    i,
    na.rm    = TRUE
  ){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = midhingeScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}


#'Trimean Scores
#'
#'@description This function computes the trimean score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{trimeanScore}}
trimeanScores <- function(
    x,
    i,
    na.rm    = TRUE
){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimeanScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}


#'Bristow Scores
#'
#'@description This function computes the Bristow score
#'from each column vector in the input matrix.
#'See \code{\link{bristowScore}} for further details.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{bristowScore}}
bristowScores <- function(
    x,
    i,
    na.rm    = TRUE
  ){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = bristowScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}


#'Reviewed Bristow Scores
#'
#'@description This function computes the reviewed Bristow score
#'from each column vector in the input matrix.
#'See \code{\link{reviewedBristowScore}} for further details.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{reviewedBristowScore}}
reviewedBristowScores <- function(
    x,
    i,
    na.rm    = TRUE
  ){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #compute
  out = apply(X = x, MARGIN = 2, FUN = reviewedBristowScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}

#
# #'PC1 Scores
# #'@inheritParams computeMatrixScore
# #'
# #'@inherit computeMatrixScore return
# #'
# #'@inherit computeMatrixScore author
# pc1Scores <- function(x, features = NULL, na.rm = T){
#
#   x = subset_features(object = x, which = features)
#
#   #rm NAs
#   if(na.rm){
#     x = stats::na.omit(object = x)
#   }
#
#   #compute PCA
#   if(!is.null(dim(x)) && all(dim(x)>0)){
#     PCA <- stats::prcomp(x = x, scale = F)
#     #select first component
#     out = PCA$x[,1]
#   } else {
#     out = NULL
#   }
#   #return
#   return(out)
# }




#'Interquartile Mean (IQM) Scores
#'
#'@description This function computes the IQM score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#@param ... further arguments to \code{\link{interquartileMeanScore}}
#'
#'@inherit computeMatrixScore return
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{interquartileMeanScore}}
interquartileMeanScores <- function(
    x,
    i,
    na.rm    = TRUE
){
  #get dim
  numCol = ncol(x)

  #subset input
  x = subsetMatrix(x = x, i = i, rm.dupli = T)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}


  #compute
  out = apply(X = x, MARGIN = 2, FUN = interquartileMeanScore, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  #return
  return(out)
}




#'Weighted Sum Scores
#'
#'@description This function computes the weighted sum
#'from each column vector in the normalised input matrix.
#'The score for the \eqn{j}-th sample is computed as:
#'
#'\deqn{wss_{j} = \sum_{i=1}^{n}  xnorm_{ij} * w_{i}}
#'
#'where \eqn{n} is the number of \code{i} elements in \code{x};
#'\eqn{xnorm_{ij}} is the \eqn{i}-th row for the \eqn{j}-th column
#'of the normalised matrix;
#'\eqn{w_{i}} is the \eqn{i}-th element of weights vector.
#'
#'@inheritParams computeMatrixScore
#@inheritParams stats::weighted.mean
#'@param normalisation the normalisation method to perform before
#'computing the weighted sum
#'@param ... further arguments to \code{normalisation}
#'
#'
#'@inherit computeMatrixScore return
#'
#'@details The input data \code{x} is firstly normalised using
#'the technique chosen via the \code{normalisation} parameter.
#'Then, for each column of \code{x}, each element of \code{i}
#'in \code{x} is weighted by the corresponding element in \code{w}.
#'Finally, the sum of the weighted elements is calculated.
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link{normaliseData}}
#'
#'@examples
#'\dontrun{
#'x = matrix(
#'data = sample(x = 1:100, size = 15),
#'nrow = 5,
#'ncol = 3,
#'dimnames = list(LETTERS[1:5], letters[1:3]))
#'
#'weightedSumScores(x = x)
#'
#'}
weightedSumScores <- function(
    x,
    i,
    # w = rep(x = 1, times = nrow(x)),
    w,
    na.rm    = TRUE,
    # rm.dupli = TRUE,
    normalisation = "quantile",
    ...
){

  #match arg
  normalisation = match.arg(normalisation)

  #check input
  if(missing(i) || is.null(i)) {i = seq_len(length.out = nrow(x))}
  if(missing(w) || is.null(w)) {w = rep(x = 1, times = length(i))}

  #get dim
  numCol = ncol(x)

  #subset input
  i = updateSig(x = x, i = i, rm.dupli = T)
  x = subsetMatrix(x = x, i = i, rm.dupli = TRUE)

  #check input
  if(isTRUE(is.null(i))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}

  #update w
  w = w[getIndex(i)]

  #Check NAs
  if(isTRUE(na.rm & any(is.na(x)))){
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
    if(isFALSE(length(w)==nrow(x))){
      stop("Error after removing NAs: length of 'w' not matching number of rows in 'x'.\n")
    }
  }

  #1) Normalise data
  x = normaliseData(x = x, method = normalisation, na.rm = na.rm, ...)

  if(isTRUE(!is.null(dim(x)))){
    if(nrow(x)>0){

      #compute
      out = apply(X = x, MARGIN = 2, FUN = function(xi, w, na.rm){
        xi = xi * w
        out = sum(xi, na.rm = na.rm)
        return(out)
      }, w = w, na.rm = na.rm)

      if(isTRUE(any(is.na(out)))){
        out[is.na(out)] = getDefaultNaValue()
      }
    } else {
      out = rep(x = getDefaultNaValue(), times = ncol(x))
    }
  } else {
    out = getDefaultNaValue()
  }

  #return
  return(out)
}




#'Single Sample Gene Set Enrichment Analysis (ssGSEA) Scores
#'
#'@description This function computes the ssGSEA score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@details This function uses internally \code{\link[GSVA]{gsva}}
#'from \code{GSVA} R package. To avoid runtime
#'errors, if \code{x} has less than 2 rows the
#'function returns \code{NA}.
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link[GSVA]{gsva}}
ssGseaScores <- function(
    x,
    i,
    na.rm   = TRUE,
    verbose = FALSE,
    ...
){

  #remove
  if(isTRUE(na.rm)){x = stats::na.omit(object = x)}

  #Remove genes with constant expression values
  x = x[!areRowValuesEqual(x),,drop=F]

  #get dim
  numRow = nrow(x)
  numCol = ncol(x)


  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(numRow < 2)){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(rownames(x))) & is.character(i)){stop("Error: 'x' must have row names.\n")}

  #check rownames
  rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")

  #check i
  if(isTRUE(is.numeric(i))){
    #check dim
    i = i[i<=numRow]
    #get names
    i = rownames(x)[i]
    #bool
    hasElements = length(i) > 0
  } else {
    hasElements = length(intersect(x = i, y = rownames(x))) > 0
  }

  #compute
  ##Check if has matching elements to avoid Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
  if(hasElements){
    out = unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = i), verbose=verbose, method = "ssgsea", ...))[,1])
  } else {
    out = getDefaultNaValue()
  }


  #return
  return(out)
}


#'Gene Set Variation Analysis (GSVA) Scores
#'
#'@description This function computes the GSVA score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'@inheritParams GSVA::gsva
#'
#'@inherit computeMatrixScore return
#'
#'@details This function uses internally \code{\link[GSVA]{gsva}}
#'from \code{GSVA} R package. To avoid runtime
#'errors, if \code{x} has less than 2 rows the
#'function returns \code{NA}.
#'
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link[GSVA]{gsva}}
gsvaScores <- function(
    x,
    i,
    na.rm   = TRUE,
    verbose = FALSE,
    kcdf    = c("Gaussian", "Poisson", "none"),
    # error.as.na = FALSE,
    ...
){

  #remove
  if(isTRUE(na.rm)){x = stats::na.omit(object = x)}

  #Remove genes with constant expression values
  x = x[!areRowValuesEqual(x),,drop=F]

  #get dim
  numRow = nrow(x)
  numCol = ncol(x)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(numRow < 2)){   return(rep(x = getDefaultNaValue(), times = numCol))}

  if(isTRUE(is.null(rownames(x))) & is.character(i)){
    # if(isFALSE(error.as.na)){
      stop("Error: 'x' must have row names.\n")
    # } else {
    #   return(rep(x = getDefaultNaValue(), times = numCol))
    # }

  }

  #check rownames
  rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")

  #check i
  if(isTRUE(is.numeric(i))){
    #check dim
    i = i[i<=numRow]
    #get names
    i = rownames(x)[i]
    #bool
    hasElements = length(i) > 0
  } else {
    hasElements = length(intersect(x = i, y = rownames(x))) > 0
  }

  #match
  kcdf = match.arg(kcdf)

  #compute
  ##Check if has matching elements to avoid Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
  if(hasElements){
    # out = tryCatch(
    #   expr = {
    #     unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = i), verbose=FALSE, method = "gsva", kcdf = kcdf, ...))[,1])
    #   }, error = function(e){
    #     if(isFALSE(error.as.na)){
    #       stop(e)
    #     } else {
    #       out = getDefaultNaValue()
    #       return(out)
    #     }
    #   }
    # )

    out = unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = i), verbose=verbose, method = "gsva", kcdf = kcdf, ...))[,1])

  } else {
    out = getDefaultNaValue()
  }

  #return
  return(out)
}


#'Pathway Level Analysis of Gene Expression (PLAGE) Scores
#'
#'@description This function computes the PLAGE score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@details This function uses internally \code{\link[GSVA]{gsva}}
#'from \code{GSVA} R package. To avoid runtime
#'errors, if \code{x} has less than 2 rows the
#'function returns \code{NA}.
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link[GSVA]{gsva}}
plageScores <- function(
    x,
    i,
    na.rm   = TRUE,
    verbose = FALSE,
    ...
){

  #remove
  if(isTRUE(na.rm)){x = stats::na.omit(object = x)}

  #Remove genes with constant expression values
  x = x[!areRowValuesEqual(x),,drop=F]

  #get dim
  numRow = nrow(x)
  numCol = ncol(x)


  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(numRow < 2)){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(rownames(x))) & is.character(i)){stop("Error: 'x' must have row names.\n")}

  #check rownames
  rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")

  #check i
  if(isTRUE(is.numeric(i))){
    #check dim
    i = i[i<=numRow]
    #get names
    i = rownames(x)[i]
    #bool
    hasElements = length(i) > 0
  } else {
    hasElements = length(intersect(x = i, y = rownames(x))) > 0
  }

  #compute
  ##Check if has matching elements to avoid Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
  if(hasElements){
    out = unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = i), verbose=verbose, method = "plage", ...))[,1])
  } else {
    out = getDefaultNaValue()
  }

  #return
  return(out)
}


#'Z-Scores
#'
#'@description This function computes the z-score
#'from each column vector in the input matrix.
#'
#'@inheritParams computeMatrixScore
#'
#'@inherit computeMatrixScore return
#'
#'@details This function uses internally \code{\link[GSVA]{gsva}}
#'from \code{GSVA} R package. To avoid runtime
#'errors, if \code{x} has less than 2 rows the
#'function returns \code{NA}.
#'
#'@inherit computeMatrixScore author
#'
#'@seealso
#'\code{\link[GSVA]{gsva}}
zScores <- function(
    x,
    i,
    na.rm   = TRUE,
    verbose = FALSE,
    ...
){

  #remove
  if(isTRUE(na.rm)){x = stats::na.omit(object = x)}

  #Remove genes with constant expression values
  x = x[!areRowValuesEqual(x),,drop=F]

  #get dim
  numRow = nrow(x)
  numCol = ncol(x)

  #check input
  if(isTRUE(is.null(x))){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(all(is.na(x)))){return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(numRow < 2)){   return(rep(x = getDefaultNaValue(), times = numCol))}
  if(isTRUE(is.null(rownames(x))) & is.character(i)){stop("Error: 'x' must have row names.\n")}

  #check rownames
  rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")

  #check i
  if(isTRUE(is.numeric(i))){
    #check dim
    i = i[i<=numRow]
    #get names
    i = rownames(x)[i]
    #bool
    hasElements = length(i) > 0
  } else {
    hasElements = length(intersect(x = i, y = rownames(x))) > 0
  }

  #compute
  ##Check if has matching elements to avoid Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr))
  if(hasElements){
    out = unlist(t(GSVA::gsva(expr = x, gset.idx.list = list(gs1 = i), verbose=verbose, method = "zscore", ...))[,1])
  } else {
    out = getDefaultNaValue()
  }

  #return
  return(out)
}
