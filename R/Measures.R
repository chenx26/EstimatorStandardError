################### This file includes all the functions to compute nonparametric sample estimators of risk/performance measures

#' Compute sample Sharpe Ratio from data
#'
#' @param data Vector of data
#' @param rf Risk-free rate
#'
#' @return Value of Sample Sharpe Ratio
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
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
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' VaR.hist(rnorm(10))
VaR.hist=function(data,alpha=0.1){
  return(-quantile(data,alpha))
}

#' Compute sample Expected Shortfall
#'
#' @param data Vector of data
#' @param alpha Tail Probability
#'
#' @return sample ES
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' ES(rnorm(10))
ES.hist=function(data,alpha=0.1){
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
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' RES(rnorm(1000))

RES=function(data,alpha=0.01,beta=0.1){
  return(-mean(data[(data<=quantile(data,beta))&(data>=quantile(data,alpha))]))
}

#' Compute sample Sortino Ratio with mean threshold
#'
#' @param data vector of data
#' @param rf risk-free rate
#'
#' @return Sample SoR
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' SoR(rnorm(10))
SoR = function(data, rf = 0){
  mu.hat = mean(data)
  sigma.minus.hat = sqrt(mean((data-mu.hat)^2*(data<=mu.hat)))
  SoR.hat = (mu.hat-rf)/sigma.minus.hat
  return(SoR.hat)
}

#' Compute sample Sortino Ratio with constant threshold
#'
#' @param data Vector of data
#' @param ... other parameters
#' @param MAR Constant Threshold
#'
#' @return SoR
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' SoR.const.IF(rnorm(10),MAR = 0.1)
SoR.const = function(data, ..., MAR = 0){
  mu.hat = mean(data)
  sigma.minus.hat = sqrt(mean((data-MAR)^2*(data<=MAR)))
  SoR.const.hat = (mu.hat-MAR)/sigma.minus.hat
  return(SoR.const.hat)
}


#' Compute sample STARR
#'
#' @param data vector of data
#' @param alpha tail probability
#' @param rf risk-free interest rate
#' @param ... other parameters
#'
#' @return sample STARR
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' STARR(rnorm(10))
STARR = function(data, ..., alpha = 0.1, rf = 0){
  mu.hat=mean(data)
  return((mu.hat - rf)/ES.hist(data, alpha = alpha))
}

#' Compute sample lower partial moment
#'
#' @param data vector of data
#' @param ... other parameters
#' @param const the constant threshold
#' @param k order of the LPM to be computed
#'
#' @return sample lower partial moment
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' LPM(rnorm(10),const = -0.1)
LPM = function(data, ..., const = 0, k = 1){
  N = length(data)
  return(1/N*sum((const-data[data<=const])^k)^(1/k))
}

#' Compute sample Omega Ratio
#'
#' @param data vector of data
#' @param ... other parameters
#' @param const the constant threshold
#'
#' @return sample Omega Ratio
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' OmegaRatio(rnorm(10),const = -0.1)
OmegaRatio = function(data, ..., const = 0){
  N = length(data)
  OmegaPlus = sum(data[data>=const]-const)/N
  OmegaMinus = sum(const-data[data<=const])/N
  return(OmegaPlus/OmegaMinus)
}
