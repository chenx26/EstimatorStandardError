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
  my.periodogram=myperiodogram(data, ...,  keep = keep)
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




