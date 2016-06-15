#### Collection of IF function computations

#' Compute the Influence Function value of Sharpe Ratio
#'
#' @param data vector of data 
#'
#' @return vector of influence function 
#' @export
#'
#' @examples
#' SR.IF(rnorm(10))
SR.IF=function(data){
  mu.hat=mean(data)
  sd.hat=sd(data)
  IF=1/sd.hat*(data-mu.hat)-1/2*mu.hat/sd.hat^3*((data-mu.hat)^2-sd.hat^2)
  return(IF)
}

VaR.IF=function(data,alpha=0.1){
  pdf.fit <- approxfun(density(data))
  qa=quantile(data,alpha)
  tmp=((data<=qa)-alpha)/pdf.fit(qa)
  return(tmp)
}





