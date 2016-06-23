#' Wrapper function to compute the standard error for VaR estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of VaR
#' @export
#'
#' @examples
#' data(edhec)
#' SE.VaR(edhec,se.method="IFiid")
SE.VaR = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.VaR.iid.xts(data,...),
               IFcor = SE.VaR.cor.xts(data,...),
               BOOTiid = SE.VaR.boot.iid.xts(data,...),
               BOOTcor = SE.VaR.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the VaR for an xts object
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#' @param alpha tail probability
#'
#' @return standard error(s) of the xts object
#' @export
#'
#' @examples
#' data(edhec)
#' SE.VaR.iid.xts(edhec)
SE.VaR.iid.xts = function(x,alpha=0.05){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.VaR.iid(x,alpha=alpha))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.VaR.iid,alpha=alpha))
  }
}


#' Standard Error of VaR for iid data
#'
#' @param data vector of data
#' @param alpha tail probability
#'
#' @return SE of VaR
#' @export
#'
#' @examples
#' SE.VaR.iid(rnorm(10))
SE.VaR.iid = function(data,alpha=0.05){
  VaR.hat=-quantile(data,alpha)
  mu=mean(data)
  N=length(data)
  variance=var(data)
  pdf.fit <- approxfun(density(data))
  V=(1-alpha)*alpha*variance/(pdf.fit(-VaR.hat))^2
  return(sqrt(V/N))
}

#' Compute standard error of VaR via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#' @param alpha tail probability
#'
#' @return SE of VaR
#' @export
#'
#' @examples
#' SE.VaR.boot.iid(rnorm(10))
SE.VaR.boot.iid = function(data, ..., alpha = 0.05){
  args=list(...)
  arg.names=names(args)
  boot.sim=100
  boot.run=100
  # check if the parameters are specified by the function call
  if(("boot.sim" %in% arg.names) & (!is.null(args["boot.sim"]))) boot.sim=args["boot.sim"]
  if(("boot.run" %in% arg.names) & (!is.null(args["boot.run"]))) boot.run=args["boot.run"]
  N=length(data)
  res=rep(0,boot.run)
  # compute boot.run standard errors
  for(run.iter in 1:boot.run){
    res.sd=rep(0,boot.sim)
    # each standard error is computed from boot.sim VaRs
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=-quantile(x,alpha)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the VaR of the boot.run standard errors
  return(mean(res))
}

