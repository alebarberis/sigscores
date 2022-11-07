#'Quantile Normalization
#'
#'@inheritParams normaliseData
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


#'Normalise Data
#'
#'@description This function is a dispatcher for the
#'selected normalisation method.
#'
#'@param x numerical matrix, features-by-samples
#'@param na.rm logical, whether to remove \code{NA}
#'values before computation
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

#'Default NA Value
#'
#'@description This function returns a default choice
#'for \code{NA} values. It is used to have uniform
#'behavior.
#'
#'@return The default \code{NA} value.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
getDefaultNaValue <- function(){NA_real_}



#'Subset Data Vector
#'
#'@description This function subset \code{x} given a
#'vector of values \code{i}. It is used to have uniform
#'behavior.
#'
#'@param x (named) numerical vector
#'@param i numerical vector giving the position in \code{x}
#'or character vector matching the names in \code{x}
#'@param rm.dupli logical, whether to remove
#'duplicated elements in \code{i}
#'
#'@return The subset vector. It returns \code{NULL}
#'if \code{i} values are not present in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'\dontrun{
#'x = setNames(object = c(1,2,3,4,5), nm = LETTERS[1:5])
#'subsetVector(x, i = c(1,3))
#'subsetVector(x, i = LETTERS[c(1,5)])
#'subsetVector(x, i = 20)
#'}
#'
#'@keywords internal
subsetVector <- function(
    x, i, rm.dupli = TRUE
  ){
  if(!missing(i) && !is.null(i)){

    #subset i
    i = updateSig(x = x, i = i, rm.dupli = rm.dupli)

    #subset
    if(length(i)>0){
      x = x[i]
    } else {
      x = NULL
    }

  }
  return(x)
}


#'Subset Data Matrix
#'
#'@description This function subset \code{x} given a
#'vector of values \code{i}. It is used to have uniform
#'behavior.
#'
#'@param x numerical matrix
#'@param i numerical vector giving the rows
#'in \code{x} or character vector matching the row
#'names in \code{x}
#'@param rm.dupli logical, whether to remove
#'duplicated elements in \code{i}. Useful argument when
#'computing random signatures
#'
#'@return The subset matrix. It returns \code{NULL}
#'if the subset is empty.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'\dontrun{
#'x = matrix(
#'   data = sample(x = 1:100, size = 15),
#'   nrow = 5,
#'   ncol = 3,
#'   dimnames = list(LETTERS[1:5], letters[1:3])
#')
#'subsetMatrix(x, i = c(1,3))
#'subsetMatrix(x, i = LETTERS[c(1,5)])
#'subsetMatrix(x, i = 20)
#'}
#'
#'@keywords internal
subsetMatrix <- function(
    x, i, rm.dupli = TRUE
){
  if(!missing(i) && !is.null(i)){

    #subset i
    i = updateSig(x = x, i = i, rm.dupli = rm.dupli)

    #subset
    if(length(i)>0){
      x = x[i,,drop=F]
    } else {
      x = NULL
    }

  }
  return(x)
}

#'Subset Signature Vector
#'
#'@description This function subset \code{i} given a
#'matrix \code{x}. It is used to have uniform
#'behavior.
#'
#'@param x numerical vector or matrix
#'@param i numerical vector giving the rows (or positions)
#'in \code{x} or character vector matching the row
#'names (or names) in \code{x}
#'@param rm.dupli logical, whether to remove
#'duplicated elements in \code{i}
#'
#'@return The subset vector. If \code{updateSig}
#'removes cases, the indices of the cases form the
#'\code{"index"} attribute of the object.
#'
#'It returns \code{NULL} if the subset is empty.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'\dontrun{
#'
#'#x is a vector
#'x = setNames(
#'   object = c(1,2,3,4,5),
#'   nm = LETTERS[1:5]
#')
#'updateSig(x, i = c(1,3))
#'updateSig(x, i = LETTERS[c(1,5)])
#'updateSig(x, i = LETTERS[c(1,6,5)])
#'updateSig(x, i = 20)
#'
#'#x is a matrix
#'x = matrix(
#'   data = sample(x = 1:100, size = 15),
#'   nrow = 5,
#'   ncol = 3,
#'   dimnames = list(LETTERS[1:5], letters[1:3])
#')
#'updateSig(x, i = c(1,3))
#'updateSig(x, i = LETTERS[c(1,5)])
#'updateSig(x, i = LETTERS[c(1,6,5)])
#'updateSig(x, i = 20)
#'}
#'
#'@keywords internal
updateSig <- function(
    x, i, rm.dupli = TRUE
){
  if(!missing(i) && !is.null(i)){

    #is matrix
    isM = is.matrix(x)

    #check input
    if(isTRUE(is.numeric(i))){
      #n elements
      n = if(isTRUE(isM)){nrow(x)}else{length(x)}
      #which
      index = which(i<=n)
    } else {
      if(isFALSE(is.character(i))){
        stop("Input vector 'i' must contain 'character' values.\n")
      }

      #names
      xnames = if(isTRUE(isM)){rownames(x)}else{names(x)}

      if(isTRUE(isM && is.null(xnames))){
        stop("Input matrix 'x' must have row names.\n")
      } else if(isTRUE(!isM && is.null(xnames))){
        stop("Input vector 'x' must have names.\n")
      }

      #which
      index = which(i %in% xnames)
    }

    #subset
    i = i[index]

    #unique
    if(isTRUE(rm.dupli)){
      keep = !duplicated(x = i)
      #update
      i = i[keep]
      index = index[keep]
    }

    #set attribute
    attr(x=i, which = "index") = index

    #check
    if(isFALSE(length(i)>0)){i = NULL}

  } else {
    i = NULL
  }
  #return
  return(i)
}


#'Index of Retained Elements in Signature Vector
#'
#'@description This function returns the \code{index}
#'attribute of an object returned from \code{updateSig}.
#'
#'@param x numerical vector, returned from \code{updateSig}
#'
#'@return The \code{index} attribute of \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'\dontrun{
#'
#'#x is a vector
#'x = setNames(
#'   object = c(1,2,3,4,5),
#'   nm = LETTERS[1:5]
#')
#'i = updateSig(x, i = LETTERS[c(1,6,5)])
#'getIndex(i)
#'
#'#x is a matrix
#'x = matrix(
#'   data = sample(x = 1:100, size = 15),
#'   nrow = 5,
#'   ncol = 3,
#'   dimnames = list(LETTERS[1:5], letters[1:3])
#')
#'i = updateSig(x, i = LETTERS[c(1,6,5)])
#'getIndex(i)
#'}
#'
#'@keywords internal
getIndex <- function(x){
  #get
  out = attr(x=x, which = "index")
  #return
  return(out)
}


#'Check Row Values Are Equal
#'
#'@description This function checks each row of the
#'input matrix to see if there are different elements.
#'
#'@param x numerical matrix
#'@return A logical vector, indicating for each row
#'if there are different elements (\code{FALSE}) or
#'if the elements are all the same (\code{TRUE}).
#'
#'If there \code{x} has one column, the function
#'returns \code{TRUE} for each row.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
areRowValuesEqual <- function(
    x){

  #check all row values are equal
  out = rowSums(x == x[,1]) == ncol(x)

  #return
  return(out)
}



#'Available Summary Scores
#'
#'@description This function returns the currently available
#'summary scores.
#'
#'@return A data frame with two columns:
#'
#'\describe{
#'\item{id}{the id of the summary score, to be used in the function calls}
#'\item{name}{the name of the summary score}
#'}
#'
#'@author Alessandro Barberis
#'
#'@examples
#'
#'getAvailableScores()
getAvailableScores <- function(){
  out = data.frame(
    id = c(
      "mean",
      "trimmedMean",
      "weightedMean",
      "IQM",
      "median",
      "mode",
      "midrange",
      "midhinge",
      "trimean",
      "briskow",
      "reviewedBriskow",
      "weightedSum",
      "ssGSEA",
      "gsva",
      "plage",
      "zscore"
      ),
    name = c(
      "Arithmetic Mean",
      "Trimmed Mean",
      "Weighted Mean",
      "Interquartile Mean",
      "Median",
      "Mode",
      "Midrange",
      "Midhinge",
      "Trimean",
      "Briskow Score",
      "Reviewed Briskow Score",
      "Weighted Sum of Quantile Normalised Data",
      "Single Sample Gene Set Enrichment Analysis",
      "Gene Set Variation Analysis",
      "Pathway Level Analysis of Gene Expression",
      "Z-Score"

    ),
    stringsAsFactors = F
  )

  return(out)
}




#'Looper
#'
#'@description Common interface to run a for loop sequentially or in parallel.
#'The core of the loop is defined by the function \code{fun}.
#'
#'@param n.iter number of iterations of the loop
#@param logger a \linkS4class{Logger}
#'@param fun A function to be executed inside the loop. It should contain an \code{iloop} parameter.
#'@param cores number of cores to use for parallel execution.
#'@param ... further arguments to \code{fun}
#'@param verbose logical, indicating whether R should report extra information on progress
#'@inheritParams foreach::foreach
#'
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
forLoop <- function(
    n.iter,
    cores    = 1L,
    .inorder = TRUE,
    fun,
    ...,
    verbose  = TRUE
  ){

  #-----------------------------------------------------------------------------------#
  #logger
  # logger = get_logger(looper)
  # logger = open_con(logger)
  # verbose = get_verbose(logger)

  #-----------------------------------------------------------------------------------#
  #Check cores
  if(isTRUE(is.null(cores) | is.na(cores)) | cores < 0){cores = 1L}

  #-----------------------------------------------------------------------------------#
  parallel = ifelse(test = cores > 1, yes = TRUE, no = FALSE)

  #-----------------------------------------------------------------------------------#
  if(isTRUE(parallel)){
    #Setup parallel
    # log_info(object = logger, message = "Detecting the number of available cores: ", sep = "", add.level = TRUE, add.time = TRUE)

    numCores.max = parallel::detectCores();

    # log_info(object = logger, message = numCores.max, sep = "\n", add.level = F, add.time = F)


    if(isTRUE(!is.na(numCores.max) && numCores.max>1)){
      # log_info(object = logger, message = "Setting up parallel environment...", sep = "", add.level = TRUE, add.time = TRUE)

      if(isTRUE(cores >= numCores.max)) {cores = numCores.max - 1}

      cl <- parallel::makeCluster(spec = cores)
      doParallel::registerDoParallel(cl);

      # log_info(object = logger, message = paste("DONE:", cores, "cores selected"), sep = "\n", add.level = F, add.time = F)
    } else {
      parallel     = FALSE;
      warning("The number of detected cores is 1 or not available. Execution set to sequential.")
    }
  }

  #-----------------------------------------------------------------------------------#
  # create list to store results
  outlist = as.list(seq(n.iter));

  #-----------------------------------------------------------------------------------#
  #execute

  if(isTRUE(parallel)){


    `%dopar%` <- foreach::`%dopar%`;
    tryCatch(expr = {

      outlist = foreach::foreach(i = seq(n.iter),
                                 .packages = c(getPackageName()),
                                 .inorder = .inorder,
                                 .verbose = verbose) %dopar%
        {
          #Fit the model
          res = do.call(what = fun, args = c(list(iloop = i), list(...)))

          #----------------------------------------------------------------------#

          # #Open connection to log file
          # logger.task = open_con(logger);
          #
          # log_info(object = logger.task, message = paste0("[",i, "] : DONE"), sep = "\n", add.level = T, add.time = T)
          # log_info(object = logger.task, message = get_log_line("long.line.1"), sep = "\n", add.level = F, add.time = F)
          #
          # #Close connection to log file
          # close_con(logger.task);
          #----------------------------------------------------------------------#

          return(res);
        }

    }, finally = {
      #Stop cluster
      # log_info(object = logger, message = "Shut down workers and stop parallel socket cluster...", sep = "", add.level = T, add.time = T)

      parallel::stopCluster(cl = cl);

      # log_info(object = logger, message = "DONE", sep = "\n", add.level = F, add.time = F)
    })

  } else {
    for(i in 1L:n.iter){
      # log_info(object = logger, message = paste("Iteration", i), sep = "\n", add.level = T, add.time = T)
      # log_debug(object = logger, message = paste("Iteration", i), sep = "\n", add.level = T, add.time = T)

      #Fit the model
      outlist[[i]] = do.call(what = fun, args = c(list(iloop = i), list(...)))
    }#END FOR(time in 1L:repeats)
  }
  #-----------------------------------------------------------------------------------#

  # close_con(logger)

  #-----------------------------------------------------------------------------------#
  return(outlist)
}



#'Package Name
#'
#'@description This function returns the name of
#'the package.
#'
#'@return A string containing the name of the
#'package.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
#'
getPackageName <- function(){
  name = "sigscores"
  return(name)
}



#'Compute Correlation
#'
#'@description This function computes the correlation
#'between summary scores.
#'
#'It is a wrapper around \code{\link[psych]{corr.test}}
#'function.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@inheritParams psych::corr.test
#'@param ... further arguments to internal function call
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link[psych]{corr.test}}
#'
#'@keywords internal
computeSigScoresCorrelation <- function(
    data,
    adjust = "BH",
    method = c("pearson", "spearman", "kendall"),
    ...
  ){

  #check input ---------------------------------------------------
  method = match.arg(method)

  #prepare data   ------------------------------------------------
  ##convert to matrix
  data = getSigScoresAsMatrix(data)

  #compute corr   ------------------------------------------------
  out = psych::corr.test(
    x      = data,
    adjust = adjust,
    method = method,
    ...
  )

  #return         ------------------------------------------------
  return(out)
}



#'Random Sampling
#'
#'@description This function takes \code{n} samples of
#'size \code{nrow(x)} from the elements of \code{x}.
#'Two types of sampling are available: sampling with
#'or without replacement.
#'
#'@param x numerical matrix, features-by-samples
#'@param method string, the sampling method to use.
#'Available options are:
#'
#'Three options are available:
#' \describe{
#'   \item{\code{none}}{a sequence from 1 to \code{nrow(x)} is returned}
#'   \item{\code{permutation}}{random sampling without replacement}
#'   \item{\code{bootstrap}}{random sampling with replacement}
#' }
#'
#'@param n.repeat integer, number of repeated samples
#'to generate
#'
#'@return A list of length \code{n} where each element
#'is a vector of indices representing the sampled data.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link[base]{sample}}
#'
#'@keywords internal
sampleData <- function(
  x,
  method = c('none', 'permutation', 'bootstrap'),
  n      = 1L
  ){

  #check input ------------------------------------------------
  ##repeats
  n = as.integer(n)
  if(isTRUE(n<0)){
    stop("Error: provided 'n' not valid. Please, insert a positive integer number (e.g. 10).\n")
  }

  ##sampling
  method = match.arg(method)

  ##dim
  nrowX = dim(x)[1]

  #sample --------------------------------------------------------
  ##get sequence
  out = seq_len(length.out = nrowX)

  ##sample
  if(isTRUE(identical(method, "none"))){
    out = list(out)
  } else {
    ###convert to bool
    replace = switch(
      method,
      "permutation" = F,
      "bootstrap"   = T
    )

    out = replicate(
      n        = n,
      expr     = sample(x = out, size = nrowX, replace = replace),
      simplify = FALSE
    )
  }

  #return         ------------------------------------------------
  return(out)
}


#'Check Function Arguments
#'
#'@description This function check if the arguments \code{args}
#'are present in the formal arguments of \code{fun}, and returns
#'only the ones that can be used.
#'
#'@param fun function
#'@param args list of arguments
#'@param ... arguments passed as \code{name = vale}
#'@return List of arguments that can be used in \code{fun}.
#'
#'@author Alessandro Barberis
#'@keywords internal
getArgsInFunFormals <- function(
    fun,
    args = NULL,
    ...,
    verbose = F){

  #get formals
  args.def = names(formals(fun = fun))

  #get one list
  args = c(args, list(...))

  #get arguments present in formals
  keep = names(args) %in% args.def

  #check passed arguments
  if(isTRUE(any(!(keep)))){
    if(isTRUE(verbose)){
      # warning("Provided argument not found in formals\n")
      warning("Provided argument not found in formals:",
              paste(setdiff(x = names(args), y = args.def), collapse = ", "),
              "\n")
    }

    args = args[keep]
  }

  #return
  return(args)
}

#'Get Signature Scores as Matrix
#'
#'@description This function takes a data frame output
#'of \code{\link{computeSigScores}} and returns a
#'numerical matrix.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'
#'@return A numerical matrix, containing the signature scores.
#'
#'@author Alessandro Barberis
#'@keywords internal
getSigScoresAsMatrix <- function(
    data
  ){

  #get scores id
  scores = getAvailableScores()$id

  #match columns
  scores = colnames(data) %in% scores

  #check
  if(isTRUE(sum(scores)>0)){
    #update data frame
    data = data[,scores,drop=F]

    #convert to matrix
    data = as.matrix(x = data, ncol = ncol(data), nrow = nrow(data))
  } else {
    stop("Error: columns in 'data' not valid. 'data' must have column names matching the score ids. Please, check your input.\n")
  }

  return(data)
}
