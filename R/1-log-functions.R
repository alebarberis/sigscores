#'@include 0-utility-functions.R
NULL

# Class Constructor -------------------------------------------------------

#'Logger Structure
#'
#'@description Creates a Logger object.
#'
#' @param path a path to the log file to use
#' @param con (optional) connection to a file
#' @param verbose logical, whether R should report extra information on progress
#' @param level the log level of the logger. Five levels are supported:
#' \describe{
#' \item{\code{"OFF"}}{indicates the logger is inactive}
#' \item{\code{"INFO"}}{indicates informational messages}
#' \item{\code{"DEBUG"}}{indicates fine-grained informational messages}
#' \item{\code{"TRACE"}}{indicates finer-grained informational events than the \code{"DEBUG"}}
#' \item{\code{"ALL"}}{all levels are considered}
#' }
#' Levels are considered ordered: \code{OFF < INFO < DEBUG < TRACE < ALL}
#'
#' @return a Logger object
#'
#' @author Alessandro Barberis
#'
#' @keywords internal
Logger <- function(path, con, verbose, level){

  if(isFALSE(is.character(path))){stop("Error: 'path' must be a character string. Please, check your input.\n")}
  if(isFALSE(is.null(con) || is.connection(con))){stop("Error: 'con' must be 'NULL' or of class 'connection'. Please, check your input.\n")}
  if(isFALSE(is.logical(verbose))){stop("Error: 'verbose' must be a boolean. Please, check your input.\n")}
  if(isFALSE(is.character(level) && (level %in% c("OFF", "INFO", "DEBUG", "TRACE", "ALL")) )){stop("Error: 'level' must be a valid character string. Please, check your input.\n")}


  return(structure(
    .Data = list(
      path    = path,
      con     = con,
      verbose = verbose,
      level   = level),
    class = "Logger")
  )
}



#' Logger Class Constructor
#'
#' @description
#' Constructor for the object.
#'
#' @details A log request in a Logger of \code{level} **l**
#' is enabled if the level of the request is **<= l**.
#'
#'
#'
#' @param path a path to the log file to use
#' @param con (optional) connection to a file
#' @param verbose logical, whether R should report extra information on progress
#' @param level the log level of the logger. Five levels are supported:
#' \describe{
#' \item{\code{"OFF"}}{indicates the logger is inactive}
#' \item{\code{"INFO"}}{indicates informational messages}
#' \item{\code{"DEBUG"}}{indicates fine-grained informational messages}
#' \item{\code{"TRACE"}}{indicates finer-grained informational events than the \code{"DEBUG"}}
#' \item{\code{"ALL"}}{all levels are considered}
#' }
#' Levels are considered ordered: \code{OFF < INFO < DEBUG < TRACE < ALL}
#'
#' @return a Logger object
#'
#'@seealso
#'\code{\link{logAll}},
#'\code{\link{logTrace}},
#'\code{\link{logDebug}},
#'\code{\link{logInfo}}
#'
#' @author Alessandro Barberis
#'
#' @export
createLogger <- function(
    path    = character(),
    con     = NULL,
    verbose = TRUE,
    level   = c("INFO", "DEBUG", "TRACE", "ALL", "OFF")
){

  #----------------------------------------------------------------------#
  level = match.arg(level)

  #----------------------------------------------------------------------#
  use.path = FALSE

  #----------------------------------------------------------------------#
  if(isTRUE(!is.null(con))){
    #check con
    bool = isOpenConnection(con)
    if(!bool){
      use.path = TRUE
    }
  }

  #----------------------------------------------------------------------#
  if(isTRUE(use.path && !is.null(path) && length(path) > 0)){
    con = checkFilePathAndOpenConnection(path = path)
  }

  #----------------------------------------------------------------------#
  return(
    Logger(path = path, con = con, verbose = verbose, level = level)
  )
}

# Getters and Setters -----------------------------------------------------
#methods

#'Logger Getters and Setters
#'
#'@name loggerGettersAndSetters
#'
#'@param object a \code{\link{Logger}}
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
NULL

#'@describeIn loggerGettersAndSetters Returns the object element
getCon     <- function(object){return(object$con)}
#'@describeIn loggerGettersAndSetters Returns the object element
getVerbose <- function(object){return(object$verbose)}
#'@describeIn loggerGettersAndSetters Returns the object element
getPath    <- function(object){return(object$path)}
#'@describeIn loggerGettersAndSetters Returns the object element
getLevel   <- function(object){return(object$level)}
#'@describeIn loggerGettersAndSetters Returns the updated object
setCon     <- function(object, value){object$con = value; return(object)}
#'@describeIn loggerGettersAndSetters Returns the updated object
setVerbose <- function(object, value){object$con = value; return(object)}
#'@describeIn loggerGettersAndSetters Returns the updated object
setPath    <- function(object, value){object$con = value; return(object)}
#'@describeIn loggerGettersAndSetters Returns the updated object
setLevel   <- function(object, value){object$con = value; return(object)}

# Log Functions  ----------------------------------------------------------

#'Get log levels
#'@description Returns all the log levels included in the selected level
#'@param level character string, a log level
#'@return A vector containing all the levels corresponding to \code{level}
#'@author Alessandro Barberis
#'@keywords internal
getLogLevels = function(
    level = c("OFF", "INFO", "DEBUG", "TRACE", "ALL")
){

  level = match.arg(level)

  out = switch(
    level,
    ALL   = c("INFO", "DEBUG", "TRACE") ,
    TRACE = c("INFO", "DEBUG", "TRACE"),
    DEBUG = c("INFO", "DEBUG"),
    INFO  = c("INFO"),
    OFF   = ""
  )

}

#'Get log levels
#'@description Returns the log \code{level} as a string with fixed length.
#'@param level character string, a log level
#'@return A vector containing all the levels corresponding to \code{level}
#'@author Alessandro Barberis
#'@keywords internal
logLevelNameFixedLength <- function(level){
  out = switch(
    level,
    ALL   = "[ALL]  ",
    TRACE = "[TRACE]",
    DEBUG = "[DEBUG]",
    INFO  = "[INFO] ",
    OFF   = ""
  )
}


#'Log Method
#'
#'@description Generic function to print the log information.
#'
#'@param object an object of class \code{\link{Logger}}
#'@param log.level a character string, the log level of the information
#'@param message a character string, the message to print
#'@param sep a character vector of strings to append after \code{message}
#'@param add.level logical, whether to add the log level to \code{message}
#'@param add.time logical, whether to add the time to \code{message}
#'
#'@return None (invisible \code{NULL}).
#'
#'@seealso
#'\code{\link{logAll}},
#'\code{\link{logTrace}},
#'\code{\link{logDebug}},
#'\code{\link{logInfo}}
#'
#'@keywords internal
#'
#'@author Alessandro Barberis
logDefault = function(object, log.level, message, sep, add.level = FALSE, add.time = FALSE){
  level = getLevel(object)

  if(isTRUE(log.level %in% getLogLevels(level = level))){

    if(isTRUE(add.time)){
      message = paste(Sys.time(), message)
    }

    if(isTRUE(add.level)){
      # message = paste0("[",log.level,"] ", message)
      message = paste(logLevelNameFixedLength(log.level), message)
    }

    printLog(object = object, message = message, sep = sep)
  }
}

#'Log Method
#'
#'@name logAll
#'
#'@description This function prints the log information of level \code{ALL}.
#'
#'@details A log information of level \code{ALL} is printed if the provided logger
#'has a >= log level.
#'
#'@inheritParams logDefault
#'@inherit logDefault return
#'
#'@seealso
#'\code{\link{logTrace}},
#'\code{\link{logDebug}},
#'\code{\link{logInfo}}
#'
#'@inherit logDefault author
#'
#'@keywords internal
logAll <- function(object, message, sep, add.level = FALSE, add.time = FALSE){
    logDefault(object = object, log.level = "ALL", message = message, sep = sep, add.level = add.level, add.time = add.time)
}

#'Log Method
#'
#'@name logTrace
#'
#'@description This function prints the log information of level \code{TRACE}.
#'
#'@details A log information of level \code{TRACE} is printed if the provided logger
#'has a >= log level.
#'
#'@inheritParams logDefault
#'@inherit logDefault return
#'
#'@seealso
#'\code{\link{logAll}},
#'\code{\link{logDebug}},
#'\code{\link{logInfo}}
#'
#'@inherit logDefault author
#'
#'@keywords internal
logTrace <- function(object, message, sep, add.level = FALSE, add.time = FALSE){
  logDefault(object = object, log.level = "TRACE", message = message, sep = sep, add.level = add.level, add.time = add.time)
}

#'Log Method
#'
#'@name logDebug
#'
#'@description This function prints the log information of level \code{DEBUG}.
#'
#'@details A log information of level \code{DEBUG} is printed if the provided logger
#'has a >= log level.
#'
#'@inheritParams logDefault
#'@inherit logDefault return
#'
#'@seealso
#'\code{\link{logAll}},
#'\code{\link{logTrace}},
#'\code{\link{logInfo}}
#'
#'@inherit logDefault author
#'
#'@keywords internal
logDebug <- function(object, message, sep, add.level = FALSE, add.time = FALSE){
  logDefault(object = object, log.level = "DEBUG", message = message, sep = sep, add.level = add.level, add.time = add.time)
}

#'Log Method
#'
#'@name logInfo
#'
#'@description This function prints the log information of level \code{INFO}.
#'
#'@details A log information of level \code{INFO} is printed if the provided logger
#'has a >= log level.
#'
#'@inheritParams logDefault
#'@inherit logDefault return
#'
#'@seealso
#'\code{\link{logAll}},
#'\code{\link{logTrace}},
#'\code{\link{logDebug}}
#'
#'@inherit logDefault author
#'
#'@keywords internal
logInfo <- function(object, message, sep, add.level = FALSE, add.time = FALSE){
  logDefault(object = object, log.level = "INFO", message = message, sep = sep, add.level = add.level, add.time = add.time)
}


#'Log Line
#'
#'@description Returns a pre-defined character vector.
#'
#'@param line character vector, the type of line to return
#'
#'@return A character vector, a line break.
#'
#'@inherit logDefault author
#'
#'@keywords internal
getLogLine <- function(line = c("long.line.1", "long.line.2", "long.line.3",
                                  "short.line.1", "short.line.2", "short.line.3")){
  line = match.arg(line);

  message = switch(line,
                   long.line.1  = "----------------------------------------------------------------------------",
                   long.line.2  = "#--------------------------------------------------------------------------#",
                   long.line.3  = "############################################################################",
                   short.line.1 = "------------------------------------------",
                   short.line.2 = "#----------------------------------------#",
                   short.line.3 = "##########################################")

  return(message)
}

# Print Functions ---------------------------------------------------------

#'Logger Print Functions
#'
#'@description This function prints the provided \code{message}
#'to the output defined via \code{object}.
#'
#'@param object a \code{\link{Logger}}
#'@param message character string, the message to print
#'@inheritParams base::cat
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
printLog <- function(object, message, sep){
  con     = getCon(object = object)
  verbose = getVerbose(object = object)

  if(isTRUE(verbose)){
    cat(message, sep = sep);
  }

  if(isTRUE(isOpenConnection(con))){
    cat(message, file = con, sep = sep);
  }

}


#'Logger Print Functions
#'
#'@describeIn printLog Prints a log line.
#'
#'@inheritParams printLog
#'@inheritParams getLogLine
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
printLogLine <- function(
    object,
    line = c("long.line.1", "long.line.2", "long.line.3",
             "short.line.1", "short.line.2", "short.line.3"),
    sep) {
  line = match.arg(line);

  message = getLogLine(line)

  printLog(object = object, message = message, sep = sep)
}
