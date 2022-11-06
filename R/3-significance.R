#'@include utility-functions.R
NULL

#'Compute the Significance of the Test
#'
#'@param obs numerical vector containing the test statistic
#'@param rnd numerical vector containing the permutation/bootstrap
#'statistics
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
#'
#'@author Alessandro Barberis
#'
#'@seealso
#'\code{\link{computeASL}},
#'\code{\link{computeCI}},
#'\code{\link{computeSE}}
#'
#'@keywords internal
computeSignificance <- function(
    obs,
    rnd,
    type       = c("asl", "ci", "se"),
    ci         = c("standard", "percentile"),
    conf.level = 0.95
){

  #match
  ci = match.arg(ci)

  #compute
  out = switch(
    type,
    'asl' = computeASL(obs = obs, rnd = rnd),
    'ci'  = computeCI(obs = obs, rnd = rnd, type = ci, conf.level = conf.level),
    'se'  = computeSE(obs = obs, rnd = rnd)
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
#'@param obs numerical vector containing the observed value of
#'the test statistic
#'@param rnd numerical vector containing the permutation/bootstrap
#'replication of the statistics
#'@param alternative a character string specifying the
#'alternative hypothesis, must be one of \code{"two.sided"} (default),
#'\code{"greater"} or \code{"less"}
#'
#'@return A numerical value representing the approximate
#'achieved significance level (ASL) of the test.
#'
#'@details It is computed as:
#'
#'\deqn{\hat{ASL}_{sample}(x) = \frac{1}{n}\sum_{i=1}^{n}  {\hat{\theta}_{i}}^{*} \ge \hat{\theta}}
#'
#'where \eqn{\hat{\theta}} is a test statistic; \eqn{{\hat{\theta}}^{*}}
#'is the permutation/bootstrap distribution of \eqn{\hat{\theta}};
#'\eqn{{\hat{\theta}_{i}}^{*}} is the computed statistic on the
#'\eqn{i}-th permutation/boostrap vector.
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
    'two.sided' = abs(rnd) > abs(obs),
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
#'of the standard error of the test statistic.
#'
#'@param x sampling replication of the test statistic
#'@inheritParams base::mean
#'
#'@return A numerical value, the estimate of the
#'standard error of the test statistic.
#'
#'@details It is computed as:
#'
#'\deqn{\hat{se}_{sample}(x) = \sqrt{ \frac{1}{n - 1}\sum_{i=1}^{n}  {\hat{\theta}_{i}}^{*} -  \hat{\theta}_{M}}^{*} }
#'
#'where \eqn{n} is the number of replications;
#'\eqn{{\hat{\theta}}^{*}} is the sampling replication of the test statistic;
#'\eqn{{\hat{\theta}_{i}}^{*}} is the computed statistic on the
#'\eqn{i}-th sampled vector;
#'\eqn{{\hat{\theta}_{M}}^{*} = \sum_{i=1}^{n}{{\hat{\theta}_{i}}^{*} / n} } is the mean
#'of the sampled statistics.
#'
#'@seealso
#'Efron, B. and Tibshirani, R.J., An Introduction to the Bootstrap, pp.47 (1994)
#'
#'@author Alessandro Barberis
#'
#'@keywords internal
computeSE <- function(x, na.rm = T){
  #compute mean
  m = mean(x = x, na.rm = na.rm)

  #number of replications
  n = length(x)

  #compute se
  se = ((x - m)^2)/(n - 1)
  se = sqrt(x = se)

  #return
  return(se)
}


computeCI <- function(
    obs,
    rnd,
    conf.level = 0.95,
    type = c("standard", "percentile")
  ){

  #compute
  out = switch(
    type,
    'standard'   = computeStandardCi(obs = obs, rnd = rnd, conf.level = conf.level),
    'percentile' = computePercentileCI(rnd = rnd, conf.level = conf.level)
  )

  #return
  return(out)
}


computePercentileCI <- function(rnd, conf.level = 0.95){
  #alpha
  alpha = 1 - conf.level

  #compute quantiles
  q1 = quantile(x = rnd, probs = (alpha / 2))
  qr = quantile(x = rnd, probs = (1 - (alpha / 2)))

  #out
  out = c(q1, qr)
  names(out) = c("lower", "upper")

  #return
  return(out)
}



# #'Bias-corrected and Accelerated Confidence Interval
# computeBCaCI <- function(x, conf.level){
#   #alpha
# }


#'Standard Confidence Interval
#Confidence interval of estimate
#'@description This function computes an approximate confidence interval of
#'a point estimate (i.e. the test statistic).
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
#'@param obs numerical vector containing the observed value of
#'the test statistic
#'@param rnd numerical vector containing the permutation/bootstrap
#'replication of the statistics
#'@param conf.level the desired confidence level
#'@param distribution sampling distribution of the estimate.
#'Use \code{normal} if the population has unknown mean and known variance (or if \code{n} is large),
#'\code{t} if population has unknown mean and variance
#'@param n sample size, used to compute the degrees of freedom if \code{distribution = "t"}
#'
#'@return A named numeric vector with 2 elements, \code{lower} and \code{upper},
#'the lower and upper bounds of the confidence interval
#'for the point estimate.
#'
#'@references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5723800/}
#@family quantifying uncertainty
#'@author Alessandro Barberis
#'
#'@keywords internal
computeStandardCi <- function(
    obs,
    rnd,
    conf.level = 0.95,
    distribution = c("normal", "t"),
    n){

  #alpha
  alpha = (1 - conf.level) / 2;

  #distribution
  distribution = match.arg(distribution)

  #se
  se = computeSE(x = rnd)

  #critical value
  if(identical(distribution, "normal")){
    Q <- stats::qnorm(p = 1 - alpha, lower.tail = TRUE);
  } else if(identical(distribution, "t")){
    Q <- stats::qt(p = 1 - alpha, df = n - 1, lower.tail = TRUE);
  }

  #Compute the interval
  interval <- c(lower = estimate - Q*se, upper = estimate + Q*se )

  #return
  return(interval);
}
