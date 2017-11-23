#' Compute long run variance by using pspectrum to compute the tapered peridograms first then use glmnet_exp
#'
#' @param data vector of data
#' @param ...
#' @param d max order of polynomial
#' @param alpha.EN
#' @param keep
#' @param standardize
#' @param return.coeffs
#' @param prewhiten
#' @param twosided
#'
#' @return
#' @export
#'
#' @examples
SE.glmnet_exp.pspectrum = function(data, ...,
                                   d=7, alpha.EN=0.5,
                                   keep=1,
                                   standardize = FALSE,
                                   return.coeffs = FALSE,
                                   prewhiten = FALSE,
                                   twosided = FALSE){
  ##### perform prewhitening
  if(prewhiten){
    res.ar = ar(data, ...)
    #    res.ar = ar(data)
    ar.coeffs = res.ar$ar
    data = na.omit(res.ar$resid)
  }

  N = length(data)
  tmp = pspectrum(data, ...)
  my.freq = tmp$freq[-1]
  my.periodogram = tmp$spec[-1]/2

  # implement cut-off
  nfreq=length(my.freq)
  my.freq=my.freq[1:floor(nfreq*keep)]
  my.periodogram=my.periodogram[1:floor(nfreq*keep)]

  if (twosided){
    my.periodogram = c(rev(my.periodogram), my.periodogram)
    my.freq = c(-rev(my.freq), my.freq)
  }

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
