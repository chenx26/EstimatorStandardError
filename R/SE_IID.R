
#' Function to compute the asymptotic standard error of the sample mean for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#'
#' @return The asymptotic standard error of the sample mean
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.mu.iid(x)

SE.mu.iid=function(data){
  mu=mean(data)
  N=length(data)
  return(sqrt(sum((data-mu)^2)))
}

#' Function to compute the asymptotic standard error of the sample standard deviation for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#'
#' @return The asymptotic standard error of the sample standard deviation
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.SD.iid(x)

SE.SD.iid=function(data){
  data=data^2
  # mu=mean(data)
  # N=length(data)
  # variance=mean((data-mu)^2)
  return(sqrt(var(data)))
}

#' Function to compute the asymptotic standard error of the sample value-at-risk for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#'
#' @return The asymptotic standard error of the sample value-at-risk
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.VaR.iid(x)


SE.VaR.iid=function(data,alpha=0.05){
  VaR.hat=-quantile(data,alpha)
  mu=mean(data)
  N=length(data)
  variance=mean((data-mu)^2)
  pdf.fit <- approxfun(density(ar1.data))
  V=(1-alpha)*alpha*variance/(pdf.fit(-VaR.hat))^2
  return(sqrt(V))
}

#' Function to compute the asymptotic standard error of the sample expected shortfall for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#' @param alpha The tail probability
#'
#' @return The asymptotic standard error of the sample expected shortfall
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.ES.iid(x,alpha=0.01)

SE.ES.iid=function(data,alpha=0.05){
  VaR.hat=-quantile(data,alpha)
  ES.hat=-mean(data[data<=-VaR.hat])
  V=mean(data^2*(data<=-VaR.hat))/alpha^2
  +(1/alpha-1)*VaR.hat^2
  +(2-2/alpha)*ES.hat*VaR.hat
  -ES.hat^2
  return(sqrt(V))
}

#' Function to compute the asymptotic standard error of the sample sharpe ratio for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#'
#' @return The asymptotic standard error of the sample sharpe ratio
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.SR.iid(x)

SE.SR.iid=function(data,rf=0){
  mu.hat=mean(data)
  sigma.hat=sqrt(mean((data-mu.hat)^2))
  mu3.hat=mean((data-mu.hat)^3)
  mu4.hat=mean((data-mu.hat)^4)
  K3.hat=mu3.hat/sigma.hat^3
  K4.hat=(mu4.hat-3*sigma.hat^4)/sigma.hat^4
  SR.hat=(mu.hat-rf)/sigma.hat
  V=1+(K4.hat+2)/4*SR.hat^2-K3.hat*SR.hat
  return(sqrt(V))
}

#' Function to compute the asymptotic standard error of the sample Sortino Ratio for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#' @param rf The risk free rate
#'
#' @return The asymptotic standard error of the sample Sortino Ratio
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.SoR.iid(x)

SE.SoR.iid=function(data,rf=0){
  mu.hat=mean(data)
  sigma.hat=sqrt(mean((data-mu.hat)^2))
  sigma.minus.hat=sqrt(mean((data-mu.hat)^2*(data<=mu.hat)))
  SoR.hat=(mu.hat-rf)/sigma.minus.hat
  mu1.minus.hat=mean((data-mu.hat)*(data<=mu.hat))
  mu3.minus.hat=mean((data-mu.hat)^3*(data<=mu.hat))
  mu4.minus.hat=mean((data-mu.hat)^4*(data<=mu.hat))
  V=sigma.hat^2/sigma.minus.hat^2*(1+SoR.hat*mu1.minus.hat/sigma.minus.hat)
  +(mu4.minus.hat/4/sigma.minus.hat^4-1/4)*SoR.hat^2
  -mu3.minus.hat/sigma.minus.hat^3*(1+SoR.hat*mu1.minus.hat/sigma.minus.hat)*SoR.hat
  return(sqrt(V))
}

#' Function to compute the asymptotic standard error of the sample STARR Ratio for the given IID data
#'
#' This function only works if the data is IID
#'
#' @param data The vector of data
#' @param rf The risk free rate
#' @param alpha The tail probability used to compute expected shortfall
#'
#' @return The asymptotic standard error of the sample STARR Ratio
#' @export
#'
#' @examples
#' x=rnorm(10)
#' SE.STARR.iid(x)

SE.STARR.iid=function(data,alpha=0.05,rf=0){
  mu.hat=mean(data)
  sigma.hat=mean((data-mu.hat)^2)
  VaR.hat=-quantile(data,alpha)
  ES.hat=-mean(data[data<=-VaR.hat])
  STARR.hat=(mu.hat-rf)/ES.hat
  A=-1/alpha*(mean((data-mu.hat)*(data+VaR.hat)*(data<=-VaR.hat)))
  B=1/alpha^2*(mean(data^2*(data<=-VaR.hat)))
  +(1/alpha-1)*VaR.hat^2
  +(2-2/alpha)*ES.hat*VaR.hat-ES.hat^2
  V=sigma.hat^2/ES.hat^2
  -2*STARR.hat/ES.hat^2*A
  +STARR.hat^2/ES.hat^2*B
  return(sqrt(V))
}

#' Asymptotic Standard Error of standard deviation for xts object
#'
#' @param x xts object
#' @param na.rm whether NAs should be omitted or not
#'
#' @return Asymptotic SE of SD for the xts object
#' @export
#'
#' @examples
#' SE.SD.iid.xts(rnorm(10))
SE.SD.iid.xts=function(x,na.rm=FALSE){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    if(na.rm) x <- na.omit(x)
    return(SE.SD.iid(x))
  }
  else {
    x <- coredata(x)
    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SD.iid))
  }
}



