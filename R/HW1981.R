#' Exact implementation of H&W 1981, compute the variance of the sample mean of the data
#'
#' @param data Vector of data
#' @param K Number of periodograms to keep
#' @param d Order of the polynomial fit
#'
#' @return Variance of the sample mean. To get asymptotic variance, multiply by length(data)
#' @export
#'
#' @examples
#' SE.HW1981(rnorm(10))
SE.HW1981=function(data,K=25,d=2){
  N=length(data)

  # Step 1: compute I(n/N) for n=1..2K and J(fn) for n=1..K
  my.periodogram=myperiodogram(data)
  my.freq=my.periodogram$freq
  my.periodogram=my.periodogram$spec
  Is=my.periodogram[1:(2*K)]
  ns=1:K
  fn=(ns*4-1)/2/N
  Js=rep(0,K)
  for(n in 1:K){
    Js[n]=log((Is[2*n-1]+Is[2*n])/2)
  }

  # Step 2: use OLS to fit polynomial

  # create 1, x, x^2, ..., x^d
  x.mat=rep(1,K)
  for(col.iter in 1:d){
    x.mat=cbind(x.mat,fn^col.iter)
  }
  x.df=as.data.frame(x.mat)
  colnames(x.df)=paste0("V",0:d)
  Js=Js+0.27
  x.df=cbind(Js,x.df)
  lm.spec=lm(Js~.-1,data=x.df)

  a0.hat=as.numeric(coef(lm.spec)[1])
  # Step 3: adjust for bias

  s=solve(t(x.mat)%*%x.mat)
  c1=exp(-0.645*s[1,1]/2)

  # Step 4: return the estimated variance
  return(c1*exp(a0.hat)/N)
}
