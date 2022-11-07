#'@include utility-functions.R
NULL

#'Plot Summary Scores as Heatmap
#'
#'@description This function generates a
#'\code{ComplexHeatmap} object.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@param scores (optional) character vector, indicating
#'the summary score(s) to plot from \code{data}
#'@param ... further arguments to internal
#'function call
#'
#'@return A \code{ComplexHeatmap} object.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{heatmapCorSigScores}}
#'
#'@examples
#'#Set seed
#'set.seed(seed = 5381L)
#'
#'#Define row/col size
#'n = 10
#'
#'#Create input matrix
#'x = matrix(
#'  data = stats::runif(n = n*n, min = 0, max = 1000),
#'  nrow = n,
#'  ncol = n,
#'  dimnames = list(
#'     paste0("g",seq(n)),
#'     paste0("S",seq(n))
#'  )
#')
#'
#'#Compute Summary Scores
#'x = computeSigScores(
#'  x = x,
#'  i = rownames(x)
#')
#'
#'#Plot scores
#'heatmapSigScores(data = x)
heatmapSigScores <- function(
    data,
    scores = NULL,

    #ComplexHeatmap arguments
    name         = "summary\nscores",
    column_title = "Summary Scores",
    ...
  ){

  #check input    ------------------------------------------------
  ##scores
  if(isTRUE(!is.null(scores) & is.character(scores))){
    #inter
    scores = intersect(x = scores, colnames(data))
    #check
    if(isFALSE(length(scores)>0)){
      stop(paste("Error: 'scores' is not valid.",
      "It should be NULL or a vector containing elements matching column names in 'data'.",
      "Please, check the provided input.\n"))
    }
  } else {
    scores = getAvailableScores()$id
  }

  #prepare data   ------------------------------------------------
  #data
  ##update data
  data = data[,(colnames(data) %in% c(scores)), drop=F]
  ##convert to matrix
  data = as.matrix(x = data, ncol = ncol(data), nrow = nrow(data))

  #plot           ------------------------------------------------
  out = heatMap(
    data         = data,
    name         = name,
    column_title = column_title,
    ...
  )

  #return         ------------------------------------------------
  return(out)
}



#'Plot Correlation of Summary Scores as Heatmap
#'
#'@description This function generates a
#'\code{ComplexHeatmap} object.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@param cor (optional) matrix containing the data to plot.
#'If missing, a correlation matrix is internally computed using
#'the default settings of function \code{\link{computeSigScoresCorrelation}}
#'@param scores (optional) character vector, indicating
#'the summary score(s) to plot from \code{data}
#'@param ... further arguments to internal
#'function call
#'
#'@return A \code{ComplexHeatmap} object.
#'
#'@author Alessandro Barberis
#'@seealso
#'\code{\link{computeSigScoresCorrelation}},
#'\code{\link{heatmapSigScores}}
#'
#'@examples
#'#Set seed
#'set.seed(seed = 5381L)
#'
#'#Define row/col size
#'n = 10
#'
#'#Create input matrix
#'x = matrix(
#'  data = stats::runif(n = n*n, min = 0, max = 1000),
#'  nrow = n,
#'  ncol = n,
#'  dimnames = list(
#'     paste0("g",seq(n)),
#'     paste0("S",seq(n))
#'  )
#')
#'
#'#Compute 5 Summary Scores
#'x = computeSigScores(
#'  x = x,
#'  i = rownames(x),
#'  scores = c("mean", "median", "mode", "midrange", "midhinge")
#')
#'
#'#Plot correlation (correlation is internally computed)
#'heatmapCorSigScores(data = x)
#'
#'
#'
#'
#'#Use simulated correlation matrix
#'
#'#Define row/col size
#'n = 5
#'
#'#Create matrix
#'x = matrix(
#'  data = 1,
#'  nrow = n,
#'  ncol = n,
#'  dimnames = list(
#'    getAvailableScores()$id[1:n],
#'    getAvailableScores()$id[1:n]
#'  )
#')
#'
#'#Create data
#'data = stats::runif(n = (((n*n)-n) / 2), min = -1, max = 1)
#'
#'#Update upper triangular
#'x[upper.tri(x = x, diag = F)] = data
#'
#'#Update lower triangular to be symmetric
#'x[lower.tri(x = x, diag = F)] = t(x)[lower.tri(x)]
#'
#'#Plot
#'heatmapCorSigScores(cor = x)
heatmapCorSigScores <- function(
    data,
    cor    = NULL,
    scores = NULL,

    #ComplexHeatmap
    name            = "test\nstatistic",
    column_title    = "Correlation",
    col             = NULL,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    ...
){

  #check input    ------------------------------------------------
  ##cor
  if(isTRUE(missing(cor) | is.null(cor))){
    if(isFALSE(missing(data) | is.null(data))){
      #compute
      cor = computeSigScoresCorrelation(
        data   = data,
        method = "pearson"
      )
      #extract statistic
      cor = cor$r

      #update
      if(isTRUE(missing(column_title))){
        column_title = "Pearson Correlation"
      }
    } else {
      stop("Error: 'cor' or 'data' must be provided. Please, check your input.\n")
    }
  }

  ##scores
  if(isTRUE(!is.null(scores) & is.character(scores))){
    #inter
    scores = intersect(x = scores, colnames(cor))
    #check
    if(isFALSE(length(scores)>0)){
      stop(paste("Error: 'scores' is not valid.",
                 "It should be NULL or a vector containing elements matching column names in 'data' or 'cor'.",
                 "Please, check the provided input.\n"))
    }
  } else {
    scores = getAvailableScores()$id
  }

  #prepare data   ------------------------------------------------
  #cor
  ##keep
  keep = colnames(cor) %in% c(scores)
  ##update data
  cor = cor[keep, keep, drop=F]
  ##convert to matrix
  cor = as.matrix(x = cor, ncol = ncol(cor), nrow = nrow(cor))


  #check          ------------------------------------------------
  #color
  if(isTRUE(missing(col))){
    col = circlize::colorRamp2(
      breaks = c(-1, 0 , 1),
      colors = c("blue", "#EEEEEE", "red")
    )
  }

  #plot           ------------------------------------------------
  out = heatMap(
    data = cor,
    name            = name,
    col             = col,
    column_title    = column_title,
    cluster_rows    = cluster_rows,
    cluster_columns = cluster_columns,
    ...
  )

  #return         ------------------------------------------------
  return(out)
}



#'Heatmap
#'
#'@description This function generates a
#'\code{ComplexHeatmap} object.
#'
#'It is a wrapper around the \code{\link[ComplexHeatmap]{Heatmap}}
#'function.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@inheritParams ComplexHeatmap::Heatmap
#'@param ... further arguments to internal
#'function call
#'
#'@return A \code{ComplexHeatmap} object.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
heatMap <- function(
    #Data
    data,

    #ComplexHeatmap default
    na_col       = "grey",
    color_space  = "LAB",
    rect_gp      = grid::gpar(col = "white"),
    col          = NULL,
    column_title = character(0),

    #Clustering
    cluster_rows    = TRUE,
    cluster_columns = TRUE,

    #Further arguments
    ...
  ){

  #prepare data   ------------------------------------------------
  ##convert to matrix
  data = as.matrix(x = data, ncol = ncol(data), nrow = nrow(data))

  #check NA       ------------------------------------------------
  #NA data
  if(isTRUE(any(is.na(data)))){
    #no cluster
    cluster_rows    = FALSE
    cluster_columns = FALSE
  }

  #plot           ------------------------------------------------
  out = ComplexHeatmap::Heatmap(
    matrix          = data,
    na_col          = na_col,
    color_space     = color_space,
    rect_gp         = rect_gp,
    col             = col,
    column_title    = column_title,
    cluster_rows    = cluster_rows,
    cluster_columns = cluster_columns,
    ...
  )
  # ComplexHeatmap::draw(out)
  #return         ------------------------------------------------
  return(out)
}



#'Plot Summary Scores as Box Plots
#'
#'@description This function generates a
#'\code{ggplot} object.
#'
#'\code{boxplotSigScores} uses the internal function
#'\code{\link{ggPlot}} which contains different function
#'calls to \code{ggplot2} functions in order to create a
#'pre-defined plot.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@param scores (optional) character vector, indicating
#'the summary score(s) to plot from \code{data}
#'@inheritParams ggPlot
#'
#'@param ... further arguments to \code{\link{ggPlot}}
#'
#'@return A \code{ggplot} object.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
#'
#'@seealso
#'\code{\link{ggPlot}}
#'
#'@examples
#'#Set seed
#'set.seed(seed = 5381L)
#'
#'#Define row/col size
#'n = 10
#'
#'#Create input matrix
#'x = matrix(
#'  data = stats::runif(n = n*n, min = 0, max = 100),
#'  nrow = n,
#'  ncol = n,
#'  dimnames = list(
#'     paste0("g",seq(n)),
#'     paste0("S",seq(n))
#'  )
#')
#'
#'#Compute Summary Scores
#'x = computeSigScores(
#'  x = x,
#'  i = rownames(x)
#')
#'
#'#Plot scores
#'boxplotSigScores(data = x)
boxplotSigScores <- function(
    #data
    data,
    scores            = NULL,

    #violin
    add.violin        = F,

    #boxplot
    add.boxplot       = T,
    boxplot.width     = 0.1,

    #points
    add.points        = F,
    point.shape       = 16,
    point.jitter.w    = 0.2,
    point.jitter.h    = NULL,

    #plot labels
    labs.title        = "Summary Scores",
    labs.x            = "Summary Scores",
    labs.y            = "Scores",
    labs.col          = "Summary\nScores",

    #theme
    axis.text.x.angle = 90,
    axis.text.x.size  = 9,
    axis.text.y.angle = 0,
    axis.text.y.size  = 9,

    #further arguments
    ...
  ) {

  #prepare data   ------------------------------------------------
  #data
  data = prepareDataForPlot(
    data   = data,
    scores = scores
  )

  #plot           ------------------------------------------------
  ##create ggplot
  out = ggPlot(
    #data
    data              = data,
    x                 = "summaryScore",
    y                 = "score",
    color             = "summaryScore",

    #violin
    add.violin        = add.violin,

    #boxplot
    add.boxplot       = add.boxplot,
    boxplot.width     = boxplot.width,

    #points
    add.points        = add.points,
    point.shape       = point.shape,
    point.jitter.w    = point.jitter.w,
    point.jitter.h    = point.jitter.h,

    #plot labels
    labs.title        = labs.title,
    labs.x            = labs.x,
    labs.y            = labs.y,
    labs.col          = labs.col,

    #theme
    axis.text.x.angle = axis.text.x.angle,
    axis.text.x.size  = axis.text.x.size,
    axis.text.y.angle = axis.text.y.angle,
    axis.text.y.size  = axis.text.y.size,
    ...
  )

  #return         ------------------------------------------------
  return(out)

}


#'Plot Summary Scores as Scatter Plots
#'
#'@description This function generates a
#'\code{ggplot} object.
#'
#'\code{scatterplotSigScores} uses the internal function
#'\code{\link{ggPlot}} which contains different function
#'calls to \code{ggplot2} functions in order to create a
#'pre-defined plot.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'@inheritParams prepareDataForPlot
#'@inheritParams ggPlot
#'
#'@param ... further arguments to \code{\link{ggPlot}}
#'
#'@return A \code{ggplot} object.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
#'
#'@seealso
#'\code{\link{ggPlot}}
#'
#'@examples
#'#Set seed
#'set.seed(seed = 5381L)
#'
#'#Define row/col size
#'n = 10
#'
#'#Create input matrix
#'x = matrix(
#'  data = stats::runif(n = n*n, min = 0, max = 100),
#'  nrow = n,
#'  ncol = n,
#'  dimnames = list(
#'     paste0("g",seq(n)),
#'     paste0("S",seq(n))
#'  )
#')
#'
#'#Compute Summary Scores
#'x = computeSigScores(
#'  x = x,
#'  i = rownames(x)
#')
#'
#'#Plot scores
#'scatterplotSigScores(data = x)
#'
#'#Plot scatter plots per summary score
#'scatterplotSigScores(
#'  data       = x,
#'  scores     = c("mean", "median"),
#'  facet.rows = "summaryScore"
#')
scatterplotSigScores <- function(
  #data
  data,
  scores            = NULL,
  runs              = NULL,

  #points
  add.points        = T,
  point.shape      = 16,
  point.jitter.w   = 0.2,
  point.jitter.h   = NULL,

  #plot labels
  labs.title        = "Summary Scores",
  labs.x            = "Sample IDs",
  labs.y            = "Scores",
  labs.col          = "Summary\nScores",

  #theme
  axis.text.x.angle = 90,
  axis.text.x.size  = 9,
  axis.text.y.angle = 0,
  axis.text.y.size  = 9,

  #further arguments
  ...
) {

  #prepare data   ------------------------------------------------
  #data
  data = prepareDataForPlot(
    data   = data,
    scores = scores,
    runs   = runs
  )

  #plot           ------------------------------------------------
  ##create ggplot
  out = ggPlot(
    #data
    data              = data,
    x                 = "sampleID",
    y                 = "score",
    color             = "summaryScore",

    #violin
    add.violin        = F,

    #boxplot
    add.boxplot       = F,

    #scores
    add.points        = add.points,
    point.shape       = point.shape,
    point.jitter.w    = point.jitter.w,
    point.jitter.h    = point.jitter.h,

    #plot labels
    labs.title        = labs.title,
    labs.x            = labs.x,
    labs.y            = labs.y,
    labs.col          = labs.col,

    #theme
    axis.text.x.angle = axis.text.x.angle,
    axis.text.x.size  = axis.text.x.size,
    axis.text.y.angle = axis.text.y.angle,
    axis.text.y.size  = axis.text.y.size,
    ...
  )

  #return         ------------------------------------------------
  return(out)

}





#'Plot Function
#'
#'@description This function generates a
#'\code{ggplot} object.
#'
#'Inside \code{ggPlot} there are different function calls
#'to \code{ggplot2} functions in order to create a pre-defined
#'plot.
#'
#'@param data \code{data.frame} containing the data to plot
# The function expects specific columns:
#
# \describe{
#   \item{\code{run}}{contains run ID}
#   \item{\code{sampleID}}{contains the sample ID}
#   \item{\code{score}}{contains the summary score for each sample}
#   \item{\code{summaryScore}}{contains summary score name}
# }
#'
#'@param violin.data,boxplot.data,point.data \code{data.frame}
#'containing the data to be displayed in the respective layers.
#'If \code{NULL}, the data is inherited from \code{data}
#'
#'@param x,y,color character string, columns in \code{data}
#'describing which variables
#'should be mapped to the aesthetics \code{x}, \code{y}, and \code{color}
#'in \code{\link[ggplot2]{ggplot}}
#'@param add.violin logical, whether to include a violin plot in the figure
#'@param add.boxplot logical, whether to include a boxplot in the figure
#'@param boxplot.width numeric, the width of boxplot
#'
#'@param add.points logical, whether to add the computed scores
#'of the individual summary measures as points in the plot
#'@param point.shape the shape to use to plot the scores. It can take
#'five types of values:
#' \itemize{
#'   \item An integer in [0, 25]
#'   \item The name of the shape
#'   \item A single character, used as a plotting symbol
#'   \item A . to draw the smallest rectangle that is visible,
#'   usually 1 pixel
#'   \item An \code{NA}, to draw nothing
#' }
#'
#' See \code{vignette("ggplot2-specs")} for further details
#'
#'@param point.jitter.w,point.jitter.h Amount of vertical and
#'horizontal jitter. The jitter is added in both positive and negative
#'directions, so the total spread is twice the value specified here.
#'See \code{\link[ggplot2]{position_jitter}} for further details
#'
#'@param labs.title The text for the title
#'@param labs.x The title of the x axis
#'@param labs.y The title of the y axis
#'@param labs.col The title of the legend
#'
#'@param axis.text.x.angle,axis.text.y.angle Specify the x and y
#'axis tick labels angles (in [0, 360])
#'See \code{\link[ggplot2]{element_text}} for further details
#'@param axis.text.x.size,axis.text.y.size Specify the x and y
#'axis tick labels size in pts.
#'See \code{\link[ggplot2]{element_text}} for further details
#'
#'@param facet.rows,facet.cols character string, variables
#'defining the faceting groups on the rows or columns
#'dimension.
#'See \code{\link[ggplot2]{facet_grid}} for further details
#'
#'@param ... further arguments to internal
#'function call
#'
#'@return A \code{ggplot} object.
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
ggPlot <- function(
    #data
    data              = NULL,

    #axes
    x,
    y,

    #color
    color,

    #violin
    add.violin        = F,
    violin.data       = NULL,

    #boxplot
    add.boxplot       = T,
    boxplot.data      = NULL,
    boxplot.width     = 0.1,

    #points
    add.points        = T,
    point.data        = NULL,
    point.shape       = 16,
    point.jitter.w    = 0.2,
    point.jitter.h    = NULL,

    #plot labels
    labs.title        = "Plots",
    labs.x            = "x",
    labs.y            = "y",
    labs.col          = character(),

    #theme
    axis.text.x.angle = 90,
    axis.text.x.size  = 9,
    axis.text.y.angle = 0,
    axis.text.y.size  = 9,

    #facets
    facet.rows        = NULL,
    facet.cols        = NULL,

    #further arguments
    ...
  ){

  #check input    ------------------------------------------------
  #check data column names
  if(isFALSE(x %in% colnames(data))){
    stop("Error: 'x' should be a column of 'data'.\n")
  }
  if(isFALSE(y %in% colnames(data))){
    stop("Error: 'y' should be a column of 'data'.\n")
  }
  if(isFALSE(color %in% colnames(data))){
    stop("Error: 'color' should be a column of 'data'.\n")
  }


  #plot           ------------------------------------------------
  ##create ggplot
  out = ggplot2::ggplot(
    data    = data,
    mapping = ggplot2::aes(
      x     = .data[[x]],
      y     = .data[[y]],
      color = .data[[color]],
    )
  )

  ##create violin plot
  if(isTRUE(add.violin)){
    out = out +
      ggplot2::geom_violin(
        data = violin.data
        # color = summaryScore
      )
  }

  ##add box plot
  if(isTRUE(add.boxplot)){
    out = out +
      ggplot2::geom_boxplot(
        data  = boxplot.data,
        width = boxplot.width
      )
  }

  ##add scatter plot
  if(isTRUE(add.points)){
    out = out +
      ggplot2::geom_jitter(
        data  = point.data,
        shape = point.shape,
        position = ggplot2::position_jitter(
          width  = point.jitter.w,
          height = point.jitter.h
        )
      )
  }

  ##theme
  out = out +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        # face  = "bold",
        # color = "#993333",
        angle = axis.text.x.angle,
        size  = axis.text.x.size
      ),
      axis.text.y = ggplot2::element_text(
        angle = axis.text.y.angle,
        size  = axis.text.y.size
      )
    )

  ##labels
  out = out +
    ggplot2::labs(
      title = labs.title,
      # x     = labs.x,
      # y     = labs.y,
      col   = labs.col
    ) +
    ggplot2::xlab(label = labs.x) +
    ggplot2::ylab(label = labs.y)

  ##facets
  if(isTRUE(!is.null(facet.rows))){
    facet.rows = ggplot2::vars(!! ggplot2::sym(facet.rows))
  }
  if(isTRUE(!is.null(facet.cols))){
    # facet.cols = ggplot2::sym(facet.cols)#NOT WORKING
    facet.cols = ggplot2::vars(!! ggplot2::sym(facet.cols))
  }
  if(isTRUE(!is.null(facet.rows) | !is.null(facet.cols))){
    out = out +
      ggplot2::facet_grid(
        rows = facet.rows,
        cols = facet.cols
      )
  }

  #return         ------------------------------------------------
  return(out)

}


#'Prepare Data For Plot
#'
#'@description This function reshape the input data frame
#'to create a long format version usable by ggplot2.
#'
#'@param data data frame, output of \code{\link{computeSigScores}}
#'
#'@param scores (optional) character vector, indicating
#'the summary score(s) to plot from \code{data}
#'
#'@param runs (optional) numeric vector, indicating the
#'repeats to plot from \code{data}
#'
#'@return A data frame in long format with 4 columns:
#' \describe{
#'   \item{\code{run}}{contains an integer run IDs; \code{0}
#'   corresponds to the signature, IDs greater than \code{0}
#'   correspond to the random signatures}
#'   \item{\code{sampleID}}{contains the sample ID}
#'   \item{\code{score}}{contains the summary score for each sample}
#'   \item{\code{summaryScore}}{contains summary score name}
#' }
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
prepareDataForPlot <- function(
  data,
  scores = NULL,
  runs   = NULL
){


  #check input    ------------------------------------------------

  ##data
  if(isTRUE(missing(data) | is.null(data))){stop("Error: 'data' not valid. Please, check your input.\n")}
  if(isFALSE(length(intersect(x = colnames(data), y = getAvailableScores()$id))>0)){
    stop(
      paste(
        "Error: 'data' not valid.",
        "It should be the output data frame obtained from 'computeSigScores'.",
        "Please, check your input.\n")
    )
  }
  if(isFALSE("sampleID" %in% colnames(data))){data$sampleID = rownames(data)}

  ##scores
  if(isTRUE(!is.null(scores) & is.character(scores))){
    #inter
    scores = intersect(x = scores, colnames(data))
    #check
    if(isFALSE(length(scores)>0)){
      stop(paste("Error: 'scores' is not valid.",
                 "It should be NULL or a vector containing elements matching column names in 'data'.",
                 "Please, check the provided input.\n"))
    }
  } else {
    scores = getAvailableScores()$id
  }

  ##run
  if(isTRUE(!is.null(runs)) & is.numeric(runs)){
    #coerce
    runs = unique(as.integer(runs))
  } else {runs = NULL}

  #reshape        ------------------------------------------------
  #check if data frame has a 'run' column
  if(isFALSE('run' %in% colnames(data))){
    ##data is output from original data
    ###update data
    data = data[,(colnames(data) %in% c(scores)), drop=F]
    ###convert to matrix
    data = as.matrix(x = data, ncol = ncol(data), nrow = nrow(data))
    ###convert to long format
    data = as.data.frame.table(x = data)
    ###rename
    colnames(data) = c("sampleID", "summaryScore", "score")
    ###add 'run' column
    data$run = 0L
  } else {
    ##data is output from repeated sampling
    ###split by group
    data  = split(
      x = data,
      f = as.factor(data$run)
    )
    ###convert to long format
    data = lapply(X = data, FUN = function(y, scores){
      #row names
      rownames(y) = y$sampleID
      #check
      y = y[,(colnames(y) %in% c(scores)), drop=F]
      #matrix
      y = as.matrix(x = y, ncol = ncol(y), nrow = nrow(y))
      #convert
      out = as.data.frame.table(y)
      #return
      return(out)
    }, scores = scores)
    ###convert to data frame
    data = data.table::setDF(x = data.table::rbindlist(l = data, use.names = T, idcol = "run"))
    ###rename
    colnames(data) = c("run", "sampleID", "summaryScore", "score")
  }

  #order
  data = data[,c("run", "sampleID", "summaryScore", "score"),drop=F]

  #update
  if(isTRUE(!is.null(runs))){
    #keep
    keep = data$run %in% runs
    #check
    if(isTRUE(any(keep))){
      #subset
      data = data[keep,,drop=F]
    } else {
      warning("'runs' value(s) not valid for selection and not considered.\n")
    }
  }

  #return         ------------------------------------------------
  return(data)
}
