#' Wrapper function to compute the standard error for Sharpe Ratio estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated Sharpe Ratio. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of SR
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SR(edhec,se.method="IFiid")
SE.SR = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.SR.iid.xts(data,...),
               IFcor = SE.SR.cor.xts(data,...),
               BOOTiid = SE.SR.boot.iid.xts(data,...),
               BOOTcor = SE.SR.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the Sharpe Raio for an xts object
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#' @param ... other parameters
#' @param rf risk free interest rate
#'
#' @return standard error(s) of the xts object
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SR.iid.xts(edhec)
SE.SR.iid.xts = function(x,...,rf=0){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SR.iid(x,rf=rf))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SR.iid,rf=rf))
  }
}


#' Standard Error of SR for iid data
#'
#' @param data vector of data
#' @param rf risk free interest rate
#'
#' @return SE of SR
#' @export
#'
#' @examples
#' SE.SR.iid(rnorm(10))
SE.SR.iid = function(data,rf=0){
  N=length(data)
  mu.hat=mean(data)
  sigma.hat=sqrt(mean((data-mu.hat)^2))
  mu3.hat=mean((data-mu.hat)^3)
  mu4.hat=mean((data-mu.hat)^4)
  K3.hat=mu3.hat/sigma.hat^3
  K4.hat=(mu4.hat-3*sigma.hat^4)/sigma.hat^4
  SR.hat=(mu.hat-rf)/sigma.hat
  V=1+(K4.hat+2)/4*SR.hat^2-K3.hat*SR.hat
  return(sqrt(V/N))
}

#' Compute standard error of SR via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#' @param rf risk free interest rate
#'
#' @return SE of SR
#' @export
#'
#' @examples
#' SE.SR.boot.iid(rnorm(10))
SE.SR.boot.iid = function(data, ...,rf=0){
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
    # each standard error is computed from boot.sim SRs
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=SR(x,rf=rf)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the mean of the boot.run standard errors
  return(mean(res))
}

