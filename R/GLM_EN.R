#' Compute the periodogram as defined in H&W 1981
#'
#' @param data Vector of data
#' @param max.freq Maximum frequency to be computed
#'
#' @return list of frequencies and corresponding periodograms
#' @export
#'
#' @examples
#' myperiodogram(rnorm(10))
myperiodogram=function(data,..., max.freq=0.5, twosided = FALSE, keep = 1){
  ## data.fft=myfft(data) This is very slow
  data.fft=fft(data)
  N=length(data)
  #   tmp=1:N
  #   inset=tmp[(1:N)<floor(N/2)]
  tmp = Mod(data.fft[2:floor(N/2)])^2/N
  tmp = sapply(tmp, function(x) max(0.00001,x))
  freq = ((1:(floor(N/2)-1))/N)
  tmp = tmp[1:floor(length(tmp) * keep)]
  freq = freq[1:floor(length(freq) * keep)]
  if (twosided){
    tmp = c(rev(tmp), tmp)
    freq = c(-rev(freq), freq)
  }
  return(list(spec=tmp,freq=freq))
}

###### function to compute the standard error of sample mean of the data, to get asymptotic SE, multiply the SE by sqrt(length(data))
###### Using GLM with L1 regularization to fit a polynomial from the periodogram
###### then the intercept of the fitted polynomial is the spectral density p(0) at frequency 0
###### the variance of of the sample mean is p(0)/length(data)
###### and the asymptotic variance of the sample mean is p(0)
###### K is the proportion of frequencies to keep, 1 means keep all

#' Compute the variance of sample mean of the data
#'
#' Using GLM with L1 regularization to fit a polynomial from the periodogram then the intercept of the fitted polynomial is the spectral density p(0) at frequency 0
#'
#' @param data Vector of data
#' @param d Maximum order of the polynomial
#' @param alpha Weight of regularization. alpha=1 LASSO, alpha=0 Ridge
#' @param keep Percentage of periodograms to use for regression
#'
#' @return variance of the sample mean. To get the asymptotic variance, multiply by length(data)
#' @export
#'
#' @examples
#' require(h2o)
#' h2o.init()
#' SE.GLM.LASSO(rnorm(10))
#' # h2o.shutdown(prompt=FALSE)
SE.GLM.LASSO=function(data,d=5,alpha=1,keep=1){

  N=length(data)
  # Step 1: compute the periodograms
  my.periodogram=myperiodogram(data, ...)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec

  # remove values of frequency 0 as it does not contain information about the variance
  # my.freq=my.freq[-1]
  # my.periodogram=my.periodogram[-1]

  # implement cut-off
  nfreq=length(my.freq)
  my.freq=my.freq[1:floor(nfreq*keep)]
  my.periodogram=my.periodogram[1:floor(nfreq*keep)]

  # Step 2: use GLM with L1 regularization to fit polynomial

  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }
  x.df=as.data.frame(x.mat)
  colnames(x.df)=paste0("V",0:d)
  x.df=cbind(my.periodogram,x.df)
  x.h2o.df=as.h2o(x.df)

  # fit GLM with L1 regularization

  my.glm.lasso=h2o.glm(x=colnames(x.h2o.df)[-1],y=colnames(x.h2o.df)[1],
                       training_frame = x.h2o.df, family = "gamma",
                       link="log",lambda_search = TRUE,alpha=alpha,
                       ignore_const_cols = FALSE)
  my.glm.lasso@model$coefficients

  # predict with new data V0=1 and all other variables=0
  newx.mat=matrix(rep(0, ncol(x.h2o.df)),nrow=1)
  newx.df=as.data.frame(newx.mat)
  newx.h2o.df=as.h2o(newx.df)
  h2o::colnames(newx.h2o.df)=h2o::colnames(x.h2o.df)
  p0.hat=h2o.predict(my.glm.lasso,newx.h2o.df)[1,1]


  # Step 4: return the estimated variance
  return(p0.hat/N)
}

#' Compute the standard error of the mean of the data using glmnet_exp
#'
#' @param data vector of data
#' @param d degree of polynomial
#' @param alpha weight for Elastic Net regularizer
#' @param keep what portion of sample spectral density to use for model fitting
#'
#' @return variance of the mean of the data
#' @export
#'

SE.glmnet_exp=function(data, ...,
                       d=7, alpha.EN=0.5,
                       keep=1,
                       standardize = FALSE,
                       return.coeffs = FALSE,
                       prewhiten = FALSE){
  ##### perform prewhitening
  if(prewhiten){
    res.ar = ar(data, ...)
#    res.ar = ar(data)
    ar.coeffs = res.ar$ar
    data = na.omit(res.ar$resid)
  }

  N=length(data)
  # Step 1: compute the periodograms
  my.periodogram=myperiodogram(data, ...,  keep)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec

  # remove values of frequency 0 as it does not contain information about the variance
  # my.freq=my.freq[-1]
  # my.periodogram=my.periodogram[-1]

  # implement cut-off
  nfreq=length(my.freq)
  # my.freq=my.freq[1:floor(nfreq*keep)]
  # my.periodogram=my.periodogram[1:floor(nfreq*keep)]

  # standardize

  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,my.freq^col.iter)
  }

  # b0 = rnorm(d + 1)

  # standardize x.mat
  mean_vec = apply(x.mat[,-1], 2, mean)
  sd_vec = apply(x.mat[,-1], 2, sd)
  if (standardize){
  for(i in 2:ncol(x.mat)){
    tmp = x.mat[,i]
    x.mat[,i] = (tmp - mean(tmp)) / sd(tmp)
  }
  }

  # fit the glmnet_exp model
  res = glmnet_exp(x.mat, my.periodogram, ..., alpha.EN = alpha.EN)
#  res = glmnet_exp(x.mat, my.periodogram, alpha = alpha)

  # Step 3: return the estimated variance, and coeffs if return.coeffs = TRUE
  if(return.coeffs){
    if(standardize){
      variance = exp(sum(res * c(1, -mean_vec / sd_vec)))/N
      if(prewhiten)
        variance = variance / ( 1 - sum(ar.coeffs))^2
      coeffs = res
      return(list(variance, coeffs))
    }
    variance = exp(res[1])/N
    if(prewhiten)
      variance = variance / ( 1 - sum(ar.coeffs))^2
    coeffs = res
    return(list(variance, coeffs))
  }


  if (standardize){
    variance = exp(sum(res * c(1, -mean_vec / sd_vec)))/N
    if(prewhiten)
      variance = variance / ( 1 - sum(ar.coeffs))^2
    return(variance)
  }
  variance = exp(res[1])/N
  if(prewhiten)
    variance = variance / ( 1 - sum(ar.coeffs))^2
  return(variance)
}




