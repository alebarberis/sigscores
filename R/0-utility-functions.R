# Default Values ----------------------------------------------------------

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
#'
#'@export
getAvailableScores <- function(){

  out = c(
    "sum"            = "Sum",
    "weightedSum"    = "Weighted Sum",
    "mean"           = "Arithmetic Mean",
    "trimmedMean"    = "Trimmed Mean",
    "weightedMean"   = "Weighted Mean",
    "median"         = "Median",
    "mode"           = "Mode",
    "midrange"       = "Midrange",
    "midhinge"       = "Midhinge",
    "trimean"        = "Trimean",
    "iqr"            = "Interquartile Range",
    "iqm"            = "Interquartile Mean",
    "mad"            = "Median Absolute Deviation",
    "aad"            = "Average Absolute Deviation",
    "ssgsea"         = "Single Sample Gene Set Enrichment Analysis",
    "gsva"           = "Gene Set Variation Analysis",
    "plage"          = "Pathway Level Analysis of Gene Expression",
    "zscore"         = "Z-Score"
  )

  out = data.frame(
    id = names(out),
    name = out,
    stringsAsFactors = F
  )

  return(out)
}


#'Available Data Transformation Methods
#'
#'@description This function returns the currently available
#'data transformation methods.
#'
#'@return A data frame with two columns:
#'
#'\describe{
#'\item{id}{the id of the data transformer, to be used in the function calls}
#'\item{name}{the name of the data transformation method}
#'}
#'
#'@author Alessandro Barberis
#'
#'@examples
#'
#'getAvailableDataTransformers()
#'
#'@export
getAvailableDataTransformers <- function(){

  out = c(
    "stepFunction" = "Step Function",
    "quantile"     = "Quantile Normalization",
    "z"            = "Z-Standardisation"
  )

  out = data.frame(
    id = names(out),
    name = out,
    stringsAsFactors = F
  )

  return(out)
}


#'Available Random Sampling Methods
#'
#'@description This function returns the currently available
#'random sampling methods.
#'
#'@return A data frame with two columns:
#'
#'\describe{
#'\item{id}{the id of the random sampling method, to be used in the function calls}
#'\item{name}{the name of the random sampling method}
#'}
#'
#'@author Alessandro Barberis
#'
#'@examples
#'
#'getAvailableRandomSamplers()
#'
#'@export
getAvailableRandomSamplers <- function(){

  out = c(
    "permutation" = "Permutation (Sampling without Replacement)",
    "bootstrap"   = "Bootstrap (Sampling with Replacement)",
    "rndsig"      = "Random Signatures",
    "rndsigsub"   = "Quasi-Random Signatures"
  )

  out = data.frame(
    id = names(out),
    name = out,
    stringsAsFactors = F
  )

  return(out)
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
getThisPackageName <- function(){
  name = "sigscores"
  return(name)
}

# Subset Data -------------------------------------------------------------

#'Subset Data
#'
#'@description This function subset \code{x} given a
#'vector of values \code{i}. It is used to have uniform
#'behavior.
#'
#'@param x (named) numerical vector or a matrix features-by-samples
#'@param i numerical vector giving the (row) position in \code{x}
#'or character vector matching the (row) names in \code{x}
#'@param rm.dupli logical, whether to remove
#'duplicated elements in \code{i}
#'
#'@return The subset data. It returns \code{NULL}
#'if \code{i} values are not present in \code{x}.
#'
#'@author Alessandro Barberis
#'
#'@examples
#'\dontrun{
#'
#'#vector
#'x = setNames(object = c(1,2,3,4,5), nm = LETTERS[1:5])
#'subsetData(x, i = c(1,3))
#'subsetData(x, i = LETTERS[c(1,5)])
#'subsetData(x, i = 20)
#'
#'#matrix
#'x = matrix(
#'   data = sample(x = 1:100, size = 15),
#'   nrow = 5,
#'   ncol = 3,
#'   dimnames = list(LETTERS[1:5], letters[1:3])
#')
#'subsetData(x, i = c(1,3))
#'subsetData(x, i = LETTERS[c(1,5)])
#'subsetData(x, i = 20)
#'}
#'
#'@keywords internal
subsetData <- function(
  x, i, rm.dupli = TRUE,
  ...){

  if(isTRUE(is.vector(x))){
    out = subsetVector(x=x,i=i,rm.dupli=rm.dupli,...)
  } else if(isTRUE(is.matrix(x))){
    out = subsetMatrix(x=x,i=i,rm.dupli=rm.dupli,...)
  } else {
    stop("Error: 'x' class is not supported.\n")
  }
  return(out)
}

# subsetData.vector <- function(x, ...){subsetVector(x=x,...)}
# subsetData.matrix  <- function(x, ...){subsetMatrix(x=x,...)}

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
      if(isFALSE(all(seq_len(length(x)) %in% getIndex(i)))){#to avoid making a copy of x when all elements are selected
        x = x[i]
      }
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
      if(isFALSE(all(seq_len(nrow(x)) %in% getIndex(i)))){#to avoid making a copy when all rows are selected
        x = x[i,,drop=F]
      }
    } else {
      x = NULL
    }

  }
  return(x)
}

#'Subset Signature Vector
#'
#'@description This function subset \code{i} given a
#'vector or matrix \code{x}. It is used to have uniform
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

# Other Functions ---------------------------------------------------------

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
  # scores = getAvailableScores()$id
  scores = setdiff(x = colnames(data), c("sampleID", "run"))

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

# Sampling Functions ------------------------------------------------------

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
#'@param n integer, number of repeated samples
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
#'@examples
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
#'#Take one sample via permutation
#'#(duplicated elements are not possible)
#'sampleData(
#'  x = x,
#'  method = 'permutation',
#'  n = 1
#')
#'
#'#Take one sample via bootstrap
#'#(duplicated elements are possible)
#'sampleData(
#'  x = x,
#'  method = 'bootstrap',
#'  n = 1
#')
#'
#'@export
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

#'Random Sampling
#'
#'@description This function takes \code{n} samples of
#'size \code{length(i)} from the elements of \code{x}.
#'Two types of sampling are available: using all
#'elements of \code{x} or using the elements of \code{x}
#'that are not part of \code{i}.
#'
#'@param x a numerical vector or matrix
#'@param i (optional) numerical vector giving the (row)
#'elements in \code{x} or character vector matching the
#'(row) names in \code{x}
#'@param n integer, number of repeated samples
#'to generate
#'@param exclude logical, whether to exclude
#'\code{i} from the possible values
#'
#'@return A list of length \code{n} where each element
#'is a vector of indices representing the sampled data.
#'
#'@author Alessandro Barberis
#'
#'@examples
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
#'#Set signature
#'i = rownames(x)[1:10]
#'
#'#Generate 1 random signatures
#'#(elements of i are possible)
#'randomSignatures(
#'  x = x,
#'  i = i,
#'  n = 1,
#'  exclude = FALSE
#')
#'
#'#Generate 1 random signatures
#'#(elements of i are excluded)
#'randomSignatures(
#'  x = x,
#'  i = i,
#'  n = 1,
#'  exclude = TRUE
#')
#'
#'@seealso
#'\code{\link[base]{sample}}
#'
#'@export
randomSignatures <- function(
  x,
  i,
  n       = 1L,
  exclude = TRUE){

  #check input ------------------------------------------------
  ##repeats
  n = as.integer(n)
  if(isTRUE(n<0)){
    stop("Error: provided 'n' not valid. Please, insert a positive integer number (e.g. 10).\n")
  }

  ##size
  size = length(i)

  ##Elements from which to choose
  universeSize = dim(x)
  universeSize = if(isTRUE(is.null(universeSize))){length(x)}else{universeSize[1]}


  #sample --------------------------------------------------------
  ##get sequence
  out = seq_len(length.out = universeSize)

  ##exclude
  if(isTRUE(exclude)){
    #get index
    out = setdiff(out, getIndex(updateSig(x=x,i=i)))
  }

  ##sample
  out = replicate(
    n        = n,
    expr     = sample(x = out, size = size, replace = F),
    simplify = FALSE
  )

  #return         ------------------------------------------------
  return(out)
}

# File and Connection -----------------------------------------------------

#'Connection Utility Functions
#'
#'@description Functions to check existence, create, open, and close
#'connections.
#'
#'@name connectionUtilityFunctions
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
NULL

#'Connection Utility Functions
#'@describeIn connectionUtilityFunctions Checks the existence of a file path
#'
#'@param path character string, the file path
#'@return A logical value, whether the path to the file exists.
#'@keywords internal
filePathExists <- function(path = NULL){
  if(isTRUE(!missing(path) && !is.null(path) && length(path)>0 && dir.exists(dirname(path)))){
    bool = TRUE;
  } else {
    bool = F;
  }
  return(bool);
}

#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Checks if \code{con}
#'is a connection.
#'@param con object to be tested
#'@return Returns \code{TRUE} or \code{FALSE} depending on whether
#'its argument is of \code{connection} class or not.
#'@keywords internal
is.connection <- function(con){
  return(isTRUE(sum(grepl(pattern = "connection", x = class(con)))==1))
}

#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Checks the existence of a
#'file path and returns an opened connection.
#'
#'@param path character string, the file path
#'@return An opened connection.
#'@keywords internal
checkFilePathAndOpenConnection <- function(path = NULL){
  if(isTRUE(filePathExists(path))){
    #Open connection to log and warning files
    con = file(description = path, open = "a");
  } else {
    con = NULL;
    # con = character();
  }
  return(con);
}


#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Returns \code{TRUE} or
#'\code{FALSE} depending on whether its argument is an opened
#'\code{connection}.
#'
#'@param con object to be tested
#'@return Returns \code{TRUE} or \code{FALSE} depending on whether
#'its argument is an opened \code{connection}.
#'@keywords internal
isOpenConnection <- function(con){

  out = FALSE;

  if(isTRUE(!missing(con) && !is.null(con))){
    if(isTRUE((sum(grepl(pattern = "connection", x = class(con)))==1) && isOpen(con, rw = ""))){
      out=TRUE;
    }
  }

  return(out);
}


#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Closes an
#'opened connection.
#'
#'@param con object to be tested
#'
#'@keywords internal
checkConAndClose <- function(con){
  if(isTRUE(isOpenConnection(con))){
    close(con);
  }
}

#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Closes an
#'opened connection.
#'
#'@param object a \code{\link{Logger}}
#'
#'@keywords internal
closeCon <- function(object){
  con = getCon(object = object)
  con = checkConAndClose(con = con)
}

#'Connection Utility Functions
#'
#'@describeIn connectionUtilityFunctions Opens a
#'connection.
#'
#'@param object a \code{\link{Logger}}
#'
#'@keywords internal
openCon <- function(object){
  fpath  = getPath(object = object)
  con    = checkFilePathAndOpenConnection(path = fpath)
  object = setCon(object = object, value = con)
  return(object)
}

# Directory ---------------------------------------------------------------

#'Check the existence of a directory
#'
#'@description This function checks the existence of a directory in input.
#'See the **Details** section below for further information.
#'
#'@param path character string, a directory path
#'
#'@return \code{TRUE} if \code{output} is a path to a directory,
#'\code{FALSE} otherwise.
#'
#'@details The function checks that the path exists.
#'If not, it creates the entire path by using
#'\code{\link[base]{dir.create}} with \code{recursive = TRUE}.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
dirExists <- function(path){
  if(!missing(path) && !is.null(path)){hasDir = T} else {hasDir = F}
  if(hasDir && !dir.exists(path)){dir.create(path, recursive = TRUE)}

  return(hasDir)
}
