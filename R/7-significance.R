#'@include 0-utility-functions.R
NULL

#'Compute Significance of a Summary Statistic
#'
#'@param obs a numeric scalar representing the observed value of the summary statistic.
#'@param rnd numerical vector containing the permutation/bootstrap
#'statistics.
#'@param type character string indicating the significance to
#'compute. Three options are available:
#' \describe{
#'   \item{\code{asl}}{achieved significance level}
#'   \item{\code{ci}}{confidence interval}
#'   \item{\code{se}}{standard error}
#' }
#'@param ci character string indicating the type of
#'confidence interval. Two options are available:
#' \describe{
#'   \item{\code{standard}}{the standard confidence interval}
#'   \item{\code{percentile}}{the percentile confidence interval}
#' }
#'@param conf.level numeric scalar, the desired confidence level. Used if \code{type = "ci"}.
#' @param ... further arguments passed to \code{\link{computeASL}},
#' \code{\link{computeCI}}, or \code{\link{computeSE}}. For example,
#' \code{alternative} for \code{type = "asl"}, \code{distribution} or \code{n}
#' for \code{type = "ci"} with \code{ci = "standard"}.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{computeASL}},
#'\code{\link{computeCI}},
#'\code{\link{computeSE}}
#'
#'@export
computeSignificance <- function(
    obs,
    rnd,
    type       = c("asl", "ci", "se"),
    ci         = c("standard", "percentile"),
    conf.level = 0.95,
    ...
){

  #match
  ci = match.arg(ci)
  type = match.arg(type)

  #check
  if (type == "asl" || (type == "ci" & ci == "standard")) {
    if (length(obs) != 1) {
      stop("Argument 'obs' must be a scalar (length 1).")
    }
  }
  if (!is.numeric(rnd) || length(rnd) < 1) {
    stop("Argument 'rnd' must be a non-empty numeric vector.")
  }

  #compute
  out = switch(
    type,
    'asl' = computeASL(obs = obs, rnd = rnd, ...),
    'ci'  = computeCI(obs = obs, rnd = rnd, type = ci, conf.level = conf.level, ...),
    'se'  = computeSE(rnd = rnd, ...)
  )

  #return
  return(out)
}


#'Achieved Significance Level (ASL)
#'
#'@description This function computes the approximate
#'achieved significance level (ASL) for a permutation
#'test or bootstrap hypothesis test.
#'
#'@param obs a numeric scalar representing the observed value of the summary statistic.
#'@param rnd numerical vector containing the permutation/bootstrap
#'replication of the statistics.
#'@param alternative a character string specifying the
#'alternative hypothesis. Must be one of \code{"two.sided"} (default),
#'\code{"greater"} or \code{"less"}.
#'
#'@return A numerical value representing the approximate
#'achieved significance level (ASL) of the test.
#'
#'@details It is computed as:
#'
#'\deqn{\hat{ASL}_{sample}(x) = \frac{1}{n}\sum_{i=1}^{n}  {\hat{\theta}_{i}}^{*} \ge \hat{\theta}}
#'
#'where \eqn{\hat{\theta}} is a summary statistic; \eqn{{\hat{\theta}}^{*}}
#'is the permutation/bootstrap distribution of \eqn{\hat{\theta}};
#'\eqn{{\hat{\theta}_{i}}^{*}} is the computed statistic on the
#'\eqn{i}-th permutation/bootstrap vector.
#'
#'\eqn{\hat{\theta}} is observed, and the ASL of the test
#'represents the probability of observing at least that large
#'a value when the null hypothesis is true.
# Roughly:
#
# \describe{
#   \item{\code{ASL < 0.10}}{borderline evidence against the null hypothesis}
#   \item{\code{ASL < 0.05}}{reasonably strong evidence against the null hypothesis}
#   \item{\code{ASL < 0.025}}{strong evidence against the null hypothesis}
#   \item{\code{ASL < 0.01}}{very strong evidence against the null hypothesis}
# }
#'
#'@seealso
#'Efron, B. and Tibshirani, R.J., An Introduction to the Bootstrap, pp.208 (1994)
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeASL <- function(
    obs,
    rnd,
    alternative = c("two.sided", "less", "greater")
){

  #match
  alternative = match.arg(alternative)

  #number of replications
  n = length(rnd)

  #compute asl
  ##alternative
  asl = switch(
    alternative,
    'two.sided' = abs(rnd) >= abs(obs),
    'less'      = rnd < obs,
    'greater'   = rnd >= obs
  )
  ##compute
  asl = sum(asl) / n

  #return
  return(asl)
}


#'Standard Error Estimate
#'
#'@description This function computes an estimate
#'of the standard error of the summary statistic.
#'
#'@param rnd sampling replication of the summary statistic
#'@inheritParams base::mean
#'
#'@return A numerical value, the estimate of the
#'standard error of the summary statistic.
#'
#'@details It is computed as:
#'
#'\deqn{\hat{se}_{sample}(x) = \sqrt{ \frac{1}{n - 1}\sum_{i=1}^{n}  {\hat{\theta}_{i}}^{*} -  \hat{\theta}_{M}}^{*} }
#'
#'where \eqn{n} is the number of replications;
#'\eqn{{\hat{\theta}}^{*}} is the sampling replication of the summary statistic;
#'\eqn{{\hat{\theta}_{i}}^{*}} is the computed statistic on the
#'\eqn{i}-th sampled vector;
#'\eqn{{\hat{\theta}_{M}}^{*} = \sum_{i=1}^{n}{{\hat{\theta}_{i}}^{*} / n} } is the mean
#'of the sampled statistics.
#'
#'@seealso
#'Efron, B. and Tibshirani, R.J., An Introduction to the Bootstrap, pp.47-49 (1994)
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeSE <- function(rnd, na.rm = FALSE){
  #compute mean
  m = mean(x = rnd, na.rm = na.rm)

  #number of replications
  n = length(rnd)

  if (n <= 1) stop("Standard error cannot be computed with less than 2 values.")

  #compute se
  se = sum((rnd - m)^2)/(n - 1)
  se = sqrt(x = se)

  #return
  return(se)
}


#'Confidence interval of estimate
#'@description This function computes an approximate confidence interval of
#'a point estimate (i.e. the test statistic).
#'
#'@param obs a numeric scalar representing the observed value of the summary
#'statistic. Used only if \code{type = "standard"}
#'@param rnd numerical vector containing the permutation/bootstrap
#'replication of the statistics.
#'@param conf.level the desired confidence level.
#'@param type character string indicating which formula to use
#'in the computation of the confidence interval. Available options are:
#' \describe{
#'   \item{\code{"standard"}}{the standard confidence interval}
#'   \item{\code{"percentile"}}{the percentile confidence interval}
#' }
#' @param na.rm logical, whether NA values should be removed before computation. Defaults to \code{FALSE}.
#' @param ... further arguments to \code{\link{computeStandardCI}} or \code{\link{computePercentileCI}}
#'
#'
#'@return A named numeric vector with 2 elements, \code{lower} and \code{upper},
#'the lower and upper bounds of the confidence interval
#'for the point estimate.
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{computeStandardCI}},
#'\code{\link{computePercentileCI}}
#'
#'@keywords internal
computeCI <- function(
    obs,
    rnd,
    conf.level = 0.95,
    type = c("standard", "percentile"),
    na.rm = FALSE,
    ...
  ){

  type = match.arg(type)

  #compute
  out = switch(
    type,
    'standard'   = computeStandardCI(obs = obs, rnd = rnd, conf.level = conf.level, na.rm = na.rm, ...),
    'percentile' = computePercentileCI(rnd = rnd, conf.level = conf.level, na.rm = na.rm, ...)
  )

  #return
  return(out)
}


#' Percentile Confidence Interval
#'
#' @description
#' Computes a percentile-based confidence interval for a summary statistic,
#' based on a vector of bootstrap or permutation replicates.
#'
#' @param rnd a numeric vector containing the bootstrap or permutation
#' replications of the summary statistic.
#' @param conf.level a numeric scalar between 0 and 1 specifying the desired
#' confidence level. Defaults to 0.95.
#' @param na.rm logical. Should missing values be removed before computing quantiles?
#' Defaults to \code{FALSE}.
#'
#' @return A named numeric vector with two elements:
#' \code{lower} and \code{upper}, corresponding to the lower and upper
#' bounds of the percentile confidence interval.
#'
#' @details
#' The percentile confidence interval is defined by taking the
#' \eqn{\alpha/2} and \eqn{1 - \alpha/2} quantiles of the bootstrap
#' distribution, where \eqn{\alpha = 1 - \text{conf.level}}.
#'
#' This method does not rely on the symmetry or normality of the
#' sampling distribution and is often used in bootstrap analysis.
#'
#' @seealso \code{\link{computeStandardCI}} for standard normal/t-based intervals.
#'
#' @references
#' Efron, B. and Tibshirani, R.J. (1994). *An Introduction to the Bootstrap*. CRC Press.
#'
#' @author Alessandro Barberis
#'
#'@keywords internal
computePercentileCI <- function(rnd, conf.level = 0.95, na.rm = FALSE){
  #alpha
  alpha = 1 - conf.level

  #compute quantiles
  q1 = stats::quantile(x = rnd, probs = (alpha / 2), na.rm = na.rm)
  qr = stats::quantile(x = rnd, probs = (1 - (alpha / 2)), na.rm = na.rm)

  #out
  out = c(q1, qr)
  names(out) = c("lower", "upper")

  #return
  return(out)
}


# TODO
#'Bias-corrected and Accelerated Confidence Interval
#'@keywords internal
computeBcaCI <- function(...){
  stop("Bias-corrected and accelerated CI not yet implemented.")
}


#'Standard Confidence Interval
#'
#'@description This function computes an approximate confidence interval of
#'a point estimate (i.e. the summary statistic).
#'
#'@details
#'The confidence interval provides additional information about the
#'variability of a point estimate and it is generally defined as
#'
#'\deqn{CI = estimate \pm \textit{margin of error} = estimate \pm \textit{critical value} \times \textit{standard error}}
#'
#'where \eqn{\textit{estimate}} is the sample statistic estimating the population parameter of interest;
#'\eqn{\textit{critical value}} is a value based on the sampling distribution of the estimate
#'and the desired confidence level; \eqn{\textit{standard error}} is the standard deviation of the
#'point estimate.
#'
#'@param obs a numeric scalar representing the observed value of the summary statistic.
#'@param rnd numerical vector containing the permutation/bootstrap
#'replication of the statistics.
#'@param conf.level the desired confidence level.
#'@param distribution sampling distribution of the estimate.
#'Use \code{"normal"} if the population has unknown mean and known variance (or
#'if the sample size \code{n} is large); \code{"t"} if the population variance is
#'unknown and sample size is small.
#'@param n sample size, used to compute the degrees of freedom if \code{distribution = "t"}.
#'@param na.rm logical, whether NA values should be removed before computing the
#'standard error. Defaults to \code{FALSE}.
#'
#'@inherit computeCI return
#'
#'@references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5723800/}
#@family quantifying uncertainty
#'@author Alessandro Barberis
#'
#'@keywords internal
computeStandardCI <- function(
    obs,
    rnd,
    conf.level = 0.95,
    distribution = c("normal", "t"),
    n,
    na.rm = FALSE){

  #alpha
  alpha = (1 - conf.level) / 2;

  #distribution
  distribution = match.arg(distribution)

  #se
  se = computeSE(rnd = rnd, na.rm = na.rm)

  #critical value
  if(identical(distribution, "normal")){
    Q <- stats::qnorm(p = 1 - alpha, lower.tail = TRUE);
  } else if(identical(distribution, "t")){
    Q <- stats::qt(p = 1 - alpha, df = n - 1, lower.tail = TRUE);
  }

  #Compute the interval
  interval <- c(lower = obs - Q*se, upper = obs + Q*se )

  #return
  return(interval);
}
