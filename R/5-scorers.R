#'@include 0-utility-functions.R 1-log-functions.R 2-parallel-functions.R 3-stats-functions.R 4-data-preprocessing.R
NULL

# Generic Scorer ----------------------------------------------------------

#'Generic Scorer Function
#'
#'@description This function computes a summary score
#'from an input vector or from each column vector
#'of an input matrix.
#'
#'@param x (named) numerical vector or matrix
#'@param i (optional) numerical vector giving the (row) position
#'in \code{x} or character vector matching the (row) names in \code{x}.
#'If \code{missing} or \code{i = NULL}, the entire \code{x} is
#'considered for the computation of the score
#'@param w numerical vector of weights the same length as
#'\code{nrow(x)}. It is used for the computation of the weighted scores
#'@param na.rm logical, whether to remove \code{NA}
#'values from \code{x} before computation
#'@param score character string indicating the summary score to compute
#'@param ... further arguments to \code{score}
#@param transform.fun (optional) data transformation function.
#If provided, it must accept \code{x} in input.
#The data is transformed before the computation
#of the score
#'@param transform.fun function to transform the data.
#'If provided, \code{x} is transformed
#'(\code{x = transform.fun(x, transform.args)})
#'before the computation of the scores.
#'See \code{\link{getDataTransformer}} for further details
#'about built-in options
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
#'@seealso
#'\code{\link{getDataTransformer}}
#'
#'@keywords internal
#'
#'@name genericScorer
NULL

# Scorers -----------------------------------------------------------------

#'Sum Scorer
#'
#'@description This scorer returns the *sum* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colSummations
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colSummations details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colSums2}}
#'
#'@export
sumScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F
){

  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = matrixStats::colSums2(
    x     = x,
    na.rm = na.rm,
    rows  = i
  )
  return(out)

}

#'Weighted Sum Scorer
#'
#'@description This scorer returns the *weighted sum* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colWeightedSums
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colWeightedSums details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colSums2}}
#'
#'@export
weightedSumScorer <- function(
    x,
    i = NULL,
    w = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){

  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
    w = w[i]
  }

  #compute
  if(isTRUE(!is.null(w))){
    ##multiply
    out = w * x
  } else {
    out = x
  }

  ##sum
  out = matrixStats::colSums2(
    x     = out,
    na.rm = na.rm
  )

  #update
  out = unlist(out)

  return(out)
}

#'Arithmetic Mean Scorer
#'
#'@description This scorer returns the *arithmetic mean* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colArithmeticMeans
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colArithmeticMeans details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colMeans2}}
#'
#'@export
meanScorer <- function(
    x,
    i     = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = matrixStats::colMeans2(
    x     = x,
    na.rm = na.rm,
    rows  = i)
  return(out)
}

#'Trimmed Mean Scorer
#'
#'@description This scorer returns the *trimmed mean* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colTrimmedMeans
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colTrimmedMeans details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{trimmedMean}}
#'
#'@export
trimmedMeanScorer <- function(
    x,
    i = NULL,
    trim  = 0,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimmedMean, trim = trim, na.rm = na.rm,simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}


#'Weighted Mean Scorer
#'
#'@description This scorer returns the *weighted mean* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colWeightedArithmeticMeans
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colWeightedArithmeticMeans details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colWeightedMeans}}
#'
#'@export
weightedMeanScorer <- function(
    x,
    i = NULL,
    w = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = matrixStats::colWeightedMeans(
    x     = x,
    w     = w,
    na.rm = na.rm,
    rows  = i)
  return(out)
}

#'Median Scorer
#'
#'@description This scorer returns the *median* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colMidpoints
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colMidpoints details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colMedians}}
#'
#'@export
medianScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #compute
  out = matrixStats::colMedians(
    x     = x,
    na.rm = na.rm,
    rows  = i)

  return(out)
}

#'Mode Scorer
#'
#'@description This scorer returns the *mode* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colModes
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colModes details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{modalValue}}
#'
#'@export
modeScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = modalValue, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}

#'Midrange Scorer
#'
#'@description This scorer returns the *midrange* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colMidranges
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colMidranges details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colRanges}}
#'
#'@export
midrangeScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
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

#'Midhinge Scorer
#'
#'@description This scorer returns the *midhinge* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colMidhinges
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colMidhinges details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{midhinge}}
#'
#'@export
midhingeScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = midhinge, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}

#'Trimean Scorer
#'
#'@description This scorer returns the *trimean* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colTrimeans
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colTrimeans details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{trimean}}
#'
#'@export
trimeanScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = trimean, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)

  return(out)
}

#'Interquartile Range (IQR) Scorer
#'
#'@description This scorer returns the *interquartile range* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colIQRs
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colIQRs details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{iqr}}
#'
#'@export
iqrScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = iqr, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}

#'Interquartile Mean (IQM) Scorer
#'
#'@description This scorer returns the *interquartile mean* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colIQMs
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colIQMs details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link{IQM}}
#'
#'@export
iqmScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }

  #compute
  out = apply(X = x, MARGIN = 2, FUN = IQM, na.rm = na.rm, simplify = FALSE)

  #update
  out = unlist(out)
  return(out)
}

#'Median Absolute Deviation (MAD) Scorer
#'
#'@description This scorer returns the *median absolute deviation* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colMADs
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colMADs details
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colMads}}
#'
#'@export
madScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  #compute
  out = matrixStats::colMads(
    x     = x,
    na.rm = na.rm,
    rows  = i)

  return(out)
}


#'Average Absolute Deviation (AAD) Scorer
#'
#'@description This scorer returns the *average absolute deviation* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colAADs
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@details
#'Internally, it uses \code{\link[matrixStats]{colMeans2}} to
#'compute the measure(s).
#'
#'@inherit genericScorer author
#'
#'@seealso
#'\code{\link[matrixStats]{colMeans2}}
#'
#'@export
aadScorer <- function(
    x,
    i = NULL,
    center = NULL,
    na.rm  = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F){
  #update
  if(isTRUE(!is.null(i))){
    x = x[i,,drop=F]
  }
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
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
    na.rm = na.rm
  )
  return(out)
}

#'Single Sample Gene Set Enrichment Analysis (ssGSEA) Scorer
#'
#'@description This scorer returns the *single sample gene set enrichment analysis (ssGSEA)* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colSsgsea
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colSsgsea details
#'
#'@inherit genericScorer author
#'
#'@inherit colEnrichment seealso
#'
#'@export
ssgseaScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F,
    ...
){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = colEnrichment(method = "ssgsea", x = x, rows = i, na.rm = na.rm, ...)
  return(out)
}

#'Gene Set Variation Analysis (GSVA) Scorer
#'
#'@description This scorer returns the *gene set variation analysis (GSVA)* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colGsva
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colGsva details
#'
#'@inherit genericScorer author
#'
#'@inherit colEnrichment seealso
#'
#'@export
gsvaScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F,
    ...){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = colEnrichment(method = "gsva", x = x, rows = i, na.rm = na.rm, ...)
  return(out)
}

#'Pathway Level Analysis of Gene Expression (PLAGE) Scorer
#'
#'@description This scorer returns the *pathway level analysis of gene expression (PLAGE)* score(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colPlage
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colPlage details
#'
#'@inherit genericScorer author
#'
#'@inherit colEnrichment seealso
#'
#'@export
plageScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F,
    ...){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = colEnrichment(method = "plage", x = x, rows = i, na.rm = na.rm, ...)
  return(out)
}

#'Z-Score Scorer
#'
#'@description This scorer returns the *z-score*(s).
#'See the **Details** section below for further information.
#'
#'@inheritParams colZscore
#'@inheritParams genericScorer
#'
#'@inherit genericScorer return
#'
#'@inherit colZscore details
#'
#'@inherit genericScorer author
#'
#'@inherit colEnrichment seealso
#'
#'@export
zscoreScorer <- function(
    x,
    i = NULL,
    na.rm = TRUE,
    transform.fun = NULL,
    transform.args = list(),
    transform.sub = F,
    ...){
  if(isTRUE(!is.null(transform.fun) && is.function(transform.fun))){
    x = do.call(what = transform.fun, args = c(list(x = x), transform.args))
  }

  out = colEnrichment(method = "zscore", x = x, rows = i, na.rm = na.rm, ...)
  return(out)

}





# Get Scorer(s) -----------------------------------------------------------

#'Get Scorer
#'
#'@description This function is a dispatcher for the score
#'function selected in input.
#'
#'@param score character string, one of the supported measures
#'@return A scoring function:
#'
#'\describe{
#' \item{\code{"sum"         }}{\code{\link{sumScorer}}}
#' \item{\code{"weightedSum" }}{\code{\link{weightedSumScorer}}}
#' \item{\code{"mean"        }}{\code{\link{meanScorer}}}
#' \item{\code{"trimmedMean" }}{\code{\link{trimmedMeanScorer}}}
#' \item{\code{"weightedMean"}}{\code{\link{weightedMeanScorer}}}
#' \item{\code{"median"      }}{\code{\link{medianScorer}}}
#' \item{\code{"mode"        }}{\code{\link{modeScorer}}}
#' \item{\code{"midrange"    }}{\code{\link{midrangeScorer}}}
#' \item{\code{"midhinge"    }}{\code{\link{midhingeScorer}}}
#' \item{\code{"trimean"     }}{\code{\link{trimeanScorer}}}
#' \item{\code{"iqr"         }}{\code{\link{iqrScorer}}}
#' \item{\code{"iqm"         }}{\code{\link{iqmScorer}}}
#' \item{\code{"mad"         }}{\code{\link{madScorer}}}
#' \item{\code{"aad"         }}{\code{\link{aadScorer}}}
#' \item{\code{"ssgsea"      }}{\code{\link{ssgseaScorer}}}
#' \item{\code{"gsva"        }}{\code{\link{gsvaScorer}}}
#' \item{\code{"plage"       }}{\code{\link{plageScorer}}}
#' \item{\code{"zscore"      }}{\code{\link{zscoreScorer}}}
#'}
#'
#'@details
#'The scoring functions are created via the generic scorer
#'\code{\link{genericScorer}} which handles vector or matrix
#'input by calling \code{\link{computeColMeasures}}.
#'
#'Internally, \code{\link{computeColMeasures}} uses these
#'functions to compute the measures:
#'\code{\link{colSummations}},
#'\code{\link{colWeightedSums}},
#'\code{\link{colArithmeticMeans}},
#'\code{\link{colTrimmedMeans}},
#'\code{\link{colWeightedArithmeticMeans}},
#'\code{\link{colMidpoints}},
#'\code{\link{colModes}},
#'\code{\link{colMidranges}},
#'\code{\link{colMidhinges}},
#'\code{\link{colTrimeans}},
#'\code{\link{colIQRs}},
#'\code{\link{colIQMs}},
#'\code{\link{colMADs}},
#'\code{\link{colAADs}},
#'\code{\link{colSsgsea}},
#'\code{\link{colGsva}},
#'\code{\link{colPlage}},
#'\code{\link{colZscore}}
#'
#'@author Alessandro Barberis
#'
#'@export
getScorer <- function(
  score = c("sum", "weightedSum",
            "mean", "trimmedMean", "weightedMean",
            "median", "mode", "midrange", "midhinge",
            "trimean", "iqr", "iqm", "mad", "aad",
            "ssgsea", "gsva", "plage", "zscore")
  ){
  score = match.arg(score)

  #compute
  out = switch(
    score,
    "sum"             = sumScorer         ,
    "weightedSum"     = weightedSumScorer ,
    "mean"            = meanScorer       ,
    "trimmedMean"     = trimmedMeanScorer ,
    "weightedMean"    = weightedMeanScorer,
    "median"          = medianScorer      ,
    "mode"            = modeScorer        ,
    "midrange"        = midrangeScorer    ,
    "midhinge"        = midhingeScorer    ,
    "trimean"         = trimeanScorer     ,
    "iqr"             = iqrScorer         ,
    "iqm"             = iqmScorer         ,
    "mad"             = madScorer         ,
    "aad"             = aadScorer         ,
    "ssgsea"          = ssgseaScorer      ,
    "gsva"            = gsvaScorer        ,
    "plage"           = plageScorer       ,
    "zscore"          = zscoreScorer
  )

  return(out)
}


#'Get Scorers
#'
#'@description This function is a dispatcher for the score
#'functions selected in input.
#'
#'@param scores character vector, the selected supported measures
#'@return A list containing the scoring functions.
#'
#'@author Alessandro Barberis
#'
#'@inherit getScorer details
#'
#'@export
getScorers <- function(
    scores = c("sum", "weightedSum",
              "mean", "trimmedMean", "weightedMean",
              "median", "mode", "midrange", "midhinge",
              "trimean", "iqr", "iqm", "mad", "aad",
              "ssgsea", "gsva", "plage", "zscore")){

  #check
  isId = scores %in% getAvailableScores()$id
  scores = scores[isId]
  if(isFALSE(all(isId))){
    if(isTRUE(length(scores)>0)){
      warning("Some score ids are not valid. Please, check your input.\n")
    } else {
      stop("Error: 'scores' is not valid. Please, check your input.\n")
    }
  }

  out = lapply(X = scores, FUN = getScorer)

  names(out) = scores

  return(out)
}

# Compute Scores ----------------------------------------------------------

#'Compute Scores
#'
#'@description This function computes summary
#'score(s) of the signature \code{i} in input
#'considering each column vector in the input matrix
#'\code{x}.
#'
#'@param x features-by-samples matrix
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
#'@param scores (optional) character vector, indicating the
#'summary score(s) to compute
#'@param scorers named list of scoring functions.
#'If provided, \code{scores} is not considered.
#'Each function must accept some specific arguments,
#'i.e. \code{x}, \code{i}, \code{na.rm}, \code{...}
#'and is expected to compute a score for each column
#'in \code{x}
#'@param args named list, where the names must match the
#'\code{scores} or the names of \code{scorers}.
#'Each element in the list is another list
#'containing the arguments to pass to the function used
#'for computing the named score. For example,
#'\code{args = list(trimmedMean = list(trim = 0.4))}
#'indicates to use \code{trim = 0.4} when computing the
#'trimmed mean scores (\code{scores = "trimmedMean"} or
#'\code{scorers = list(trimmedMean = getScorer("trimmedMean"))})
#'@param sample.id logical, whether to report the
#'sample ID as a column in the output data frame
#'@inheritParams forLoop
#'@param logger (optional) a \code{\link{Logger}} object.
#'If provided, it will be used to report extra information
#'on progress. To create a Logger use \code{\link{createLogger}}
#'@param filename (optional) character string, a name without
#'extension for the output file
#'@param outdir (optional) character string, path to
#'the output directory. If provided the returned data
#'will be stored
#'
#'@return A data frame containing the computed
#'score(s) for each sample. Each row corresponds to
#'a different sample.
#'
#'@author Alessandro Barberis
#'
#'@details
#'\code{\link{computeScores}} uses the provided scorers
#'to calculate the measures.
#'
#'If \code{scorers = NULL} and \code{scores} are provided,
#'the scorers are retrieved by using the function
#'\code{\link{getScorers}}:
#'
#'\describe{
#' \item{\code{"sum"         }}{\code{\link{sumScorer}}}
#' \item{\code{"weightedSum" }}{\code{\link{weightedSumScorer}}}
#' \item{\code{"mean"        }}{\code{\link{meanScorer}}}
#' \item{\code{"trimmedMean" }}{\code{\link{trimmedMeanScorer}}}
#' \item{\code{"weightedMean"}}{\code{\link{weightedMeanScorer}}}
#' \item{\code{"median"      }}{\code{\link{medianScorer}}}
#' \item{\code{"mode"        }}{\code{\link{modeScorer}}}
#' \item{\code{"midrange"    }}{\code{\link{midrangeScorer}}}
#' \item{\code{"midhinge"    }}{\code{\link{midhingeScorer}}}
#' \item{\code{"trimean"     }}{\code{\link{trimeanScorer}}}
#' \item{\code{"iqr"         }}{\code{\link{iqrScorer}}}
#' \item{\code{"iqm"         }}{\code{\link{iqmScorer}}}
#' \item{\code{"mad"         }}{\code{\link{madScorer}}}
#' \item{\code{"aad"         }}{\code{\link{aadScorer}}}
#' \item{\code{"ssgsea"      }}{\code{\link{ssgseaScorer}}}
#' \item{\code{"gsva"        }}{\code{\link{gsvaScorer}}}
#' \item{\code{"plage"       }}{\code{\link{plageScorer}}}
#' \item{\code{"zscore"      }}{\code{\link{zscoreScorer}}}
#'}
#'
#'Look at the different functions to know which specific
#'arguments they accept (arguments can be passed via the
#'\code{args} parameter).
#'
#'Scorers also accepts a transformation function
#'via the \code{transform.fun} argument, which
#'is used to transform the data before the computation
#'of the scores so that:
#'\code{x = transform.fun(x = x, transform.args)},
#'where \code{transform.args} is a list of parameters passed
#'to the transformation function.
#'Look at the different functions for further details.
#'
#'@examples
#'\dontrun{
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
#'#Compute all scores
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10]
#')
#'
#'#Compute all scores using 'scorers' argument
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'     'score1' = getScorer('sum'),
#'     'score2' = getScorer('mean')
#'  )
#')
#'
#'#Pass different parameters to the same scorer
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'   'score1' = getScorer('trimmedMean'),
#'   'score2' = getScorer('trimmedMean')
#'  ),
#'  args = list(
#'   'score1' = list(trim = 0),
#'   'score2' = list(trim = 0.3)
#'  )
#')
#'
#'#Transform data and compute the scores
#'computeScores(
#'  x = x,
#'  i = rownames(x)[1:10],
#'  scorers = list(
#'   'score1' = getScorer('weightedSum'),
#'   'score2' = getScorer('mean')
#'  ),
#'  args = list(
#'   'score1' = list(transform.fun = getDataTransformer('quantile')),
#'   'score2' = list(
#'      transform.fun = getDataTransformer('stepFunction'),
#'      transform.args = list(
#'        method = 'median',
#'        by = 'rows'
#'      )
#'    )
#'  )
#')
#'}
#'
#'@keywords internal
computeScores <- function(
    x,
    i         = NULL,
    na.rm     = TRUE,
    scorers   = NULL,
    scores    = NULL,
    args      = NULL,
    sample.id = TRUE,
    cores     = 1L,
    logger    = NULL,
    filename  = "sigscores",
    outdir    = NULL
){

  #check input ------------------------------------------------
  ##x
  if(isFALSE(is.matrix(x))){
    x = matrix(data = x, ncol = 1, dimnames = list(names(x)))
  }
  ##sample names
  if(isTRUE(is.null(colnames(x)))){colnames(x) = colnames(x = x, do.NULL = FALSE, prefix = "S")}
  ##gene names
  if(isTRUE(is.null(rownames(x)))){rownames(x) = rownames(x = x, do.NULL = FALSE, prefix = "g")}

  ##scorers
  if(isTRUE(missing(scorers) | is.null(scorers))){
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
    scorers = getScorers(scores = scores)
    areScoresInternal = T
  } else {
    if(isTRUE(is.list(scorers))){
      if(isFALSE(all(sapply(scorers, FUN = is.function)))){
        stop("Error: provided 'scorers' not valid. Each element of the list must be a function computing a score.\n")
      }
    } else {
      stop("Error: provided 'scorers' not valid. It should be a named list, where each element is a function computing a score.\n")
    }
    areScoresInternal = F
  }
  ##args
  if(isFALSE(missing(args) | is.null(args))){
    if(isTRUE(is.list(args) && !is.null(names(args)) && (any(names(args) %in% names(scorers))))){
      #update
      args = args[(names(args) %in% names(scorers))]
    } else{
      stop("Error: provided 'args' not valid. It must be a named list, with names matching the provided 'scorers'.\n")
    }
  }

  ##signature
  if(isTRUE(missing(i) | is.null(i))){
    i = seq_len(length.out = nrow(x))
  } else {
    #check i
    if(isTRUE(is.character(i))){
      #get as integer
      i = match(table = rownames(x), x = i)
    }
  }
  ##get size
  n = length(i)

  ##logger
  logger = openCon(logger)

  #compute scores ------------------------------------------------
  logDebug(object = logger, message = "Starting the computation.", sep = "\n", add.level = TRUE, add.time = TRUE)

  ##loop over scores
  ###number of iteration
  n.iter = length(scorers)
  ###compute scores
  out = forLoop(
    n.iter   = n.iter,
    cores    = cores,
    .inorder = T,
    fun = function(iloop, scorers, args, ...) {
      #compute
      o = tryCatch(
        expr = {
          do.call(
            what = scorers[[iloop]],
            args = c(
              list(...),
              args[[names(scorers)[[iloop]]]]
            )
          )
        },
        error = function(err){
          o = rep(x = getDefaultNaValue(), times = ncol(list(...)[['x']]))
          return(o)
        }
      )

      #return
      return(o)
    },
    x      = x,
    na.rm  = na.rm,
    scorers = scorers,
    i      = i,
    logger = logger,
    # w      = w,
    args   = args
  )

  logDebug(object = logger, message = "Computation finished.", sep = "\n", add.level = TRUE, add.time = TRUE)

  #output --------------------------------------------------------

  logDebug(object = logger, message = "Creating output object...", sep = "", add.level = TRUE, add.time = TRUE)

  ##set names
  names(out) = names(scorers)

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

  logDebug(object = logger, message = "DONE", sep = "\n", add.level = F, add.time = F)


  #save         --------------------------------------------------
  saveRes = dirExists(outdir);
  if(isTRUE(saveRes)){
    if(isTRUE(!missing(filename) & is.character(filename))){
      #Save locally
      logDebug(object = logger, message = "Saving the results...", sep = "", add.level = TRUE, add.time = TRUE)
      saveRDS(object = out, file = file.path(outdir, paste0(filename, ".rds")));
      logDebug(object = logger, message = "DONE.", sep = "\n", add.level = F, add.time = F)

    } else {
      warning("Output directory was provided but 'filename' is missing without default. Data won't be saved.\n")
    }
  }

  #clean ---------------------------------------------------------
  ##logger
  closeCon(logger)

  #return         ------------------------------------------------
  return(out)
}
