################### This file includes all the functions to compute nonparametric sample estimators of risk/performance measures

#' Compute sample Sharpe Ratio from data
#'
#' @param data Vector of data
#' @param rf Risk-free rate
#'
#' @return Value of Sample Sharpe Ratio
#' @export
#'
#' @examples
#' SR(rnorm(10))
SR=function(data,rf=0){
  mu.hat=mean(data)
  sigma.hat=sd(data)
  return((mu.hat-rf)/sigma.hat)
}


#' Compute sample value-at-risk
#'
#' @param data Vector of data
#' @param alpha Tail probability
#'
#' @return Sample VaR
#' @export
#'
#' @examples
#' VaR(rnorm(10))
VaR=function(data,alpha=0.1){
  return(-quantile(data,alpha))
}

#' Compute sample Expected Shortfall
#'
#' @param data Vector of data
#' @param alpha Tail Probability
#'
#' @return sample ES
#' @export
#'
#' @examples
#' ES(rnorm(10))
ES=function(data,alpha=0.1){
  return(-mean(data[data<=quantile(data,alpha)]))
}

#' Compute sample Robust Expected Shortfall
#'
#' @param data Vector of data
#' @param alpha Lower Tail Probability
#' @param beta Upper Tail Probability
#'
#' @return sample RES
#' @export
#'
#' @examples
#' RES(rnorm(1000))

RES=function(data,alpha=0.01,beta=0.1){
  return(-mean(data[(data<=quantile(data,beta))&(data>=quantile(data,alpha))]))
}

