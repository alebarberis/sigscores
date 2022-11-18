#'@include 0-utility-functions.R 1-log-functions.R
NULL

#'Looper
#'
#'@description Common interface to run a for loop sequentially or in parallel.
#'The core of the loop is defined by the function \code{fun}.
#'
#'@param n.iter number of iterations of the loop
#'@param fun A function to be executed inside the loop. It should contain an \code{iloop} parameter.
#'@param cores number of cores to use for parallel execution.
#'@param ... further arguments to \code{fun}
#'@param logger a \code{\link{Logger}}
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
    logger   = NULL,
    verbose  = TRUE
){

  #-----------------------------------------------------------------------------------#
  #logger
  if(isTRUE(is.null(logger))){logger = createLogger(verbose = FALSE)}
  logger = openCon(logger)
  # verbose = get_verbose(logger)
  logDebug(object = logger, message = getLogLine("short.line.1"), sep = "\n", add.level = F, add.time = F)
  #-----------------------------------------------------------------------------------#
  #Check cores
  if(isTRUE(is.null(cores) | is.na(cores)) | cores < 0){cores = 1L}

  #-----------------------------------------------------------------------------------#
  parallel = ifelse(test = cores > 1, yes = TRUE, no = FALSE)

  #-----------------------------------------------------------------------------------#
  if(isTRUE(parallel)){
    #Setup parallel
    logDebug(object = logger, message = "Detecting the number of available cores: ", sep = "", add.level = TRUE, add.time = TRUE)

    numCores.max = parallel::detectCores();

    logDebug(object = logger, message = numCores.max, sep = "\n", add.level = F, add.time = F)


    if(isTRUE(!is.na(numCores.max) && numCores.max>1)){
      logDebug(object = logger, message = "Setting up parallel environment...", sep = "", add.level = TRUE, add.time = TRUE)

      if(isTRUE(cores >= numCores.max)) {cores = numCores.max - 1}

      cl <- parallel::makeCluster(spec = cores)
      doParallel::registerDoParallel(cl);

      logDebug(object = logger, message = paste("DONE:", cores, "cores selected"), sep = "\n", add.level = F, add.time = F)
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
                                 .packages = c(getThisPackageName()),
                                 .inorder = .inorder,
                                 .verbose = verbose) %dopar%
        {
          #Fit the model
          res = do.call(what = fun, args = c(list(iloop = i), list(...)))

          #----------------------------------------------------------------------#

          #Open connection to log file
          logger.task = openCon(logger);

          logDebug(object = logger.task, message = paste0("[",i, "] : DONE"), sep = "\n", add.level = T, add.time = T)

          #Close connection to log file
          closeCon(logger.task);
          #----------------------------------------------------------------------#

          return(res);
        }

    }, finally = {
      #Stop cluster
      logDebug(object = logger, message = "Shut down workers and stop parallel socket cluster...", sep = "", add.level = T, add.time = T)

      parallel::stopCluster(cl = cl);

      logDebug(object = logger, message = "DONE", sep = "\n", add.level = F, add.time = F)
    })

  } else {
    for(i in 1L:n.iter){
      logDebug(object = logger, message = paste("Iteration", i), sep = "\n", add.level = T, add.time = T)

      #Fit the model
      outlist[[i]] = do.call(what = fun, args = c(list(iloop = i), list(...)))
    }#END FOR(time in 1L:repeats)
  }
  #-----------------------------------------------------------------------------------#
  logDebug(object = logger, message = getLogLine("short.line.1"), sep = "\n", add.level = F, add.time = F)
  closeCon(logger)

  #-----------------------------------------------------------------------------------#
  return(outlist)
}
