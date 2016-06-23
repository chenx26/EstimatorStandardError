#' Wrapper function to compute the standard error for STARR estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of STARR
#' @export
#'
#' @examples
#' data(edhec)
#' SE.STARR(edhec,se.method="IFiid")
SE.STARR = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.STARR.iid.xts(data,...),
               IFcor = SE.STARR.cor.xts(data,...),
               BOOTiid = SE.STARR.boot.iid.xts(data,...),
               BOOTcor = SE.STARR.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the STARR for an xts object
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#'
#' @return standard error(s) of the xts object
#' @export
#'
#' @examples
#' data(edhec)
#' SE.STARR.iid.xts(edhec)
SE.STARR.iid.xts = function(x,alpha=0.05, rf = 0){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.STARR.iid(x,alpha=alpha))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.STARR.iid,alpha=alpha))
  }
}


#' Standard Error of STARR for iid data
#'
#' @param data vector of data
#'
#' @return SE of STARR
#' @export
#'
#' @examples
#' SE.STARR.iid(rnorm(10))
SE.STARR.iid = function(data,alpha=0.05, rf = 0){
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
  return(sqrt(V/length(data)))
}

#' Compute standard error of STARR via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#'
#' @return SE of STARR
#' @export
#'
#' @examples
#' SE.STARR.boot.iid(rnorm(10))
SE.STARR.boot.iid = function(data, ..., alpha = 0.05, rf = 0){
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
    # each standard error is computed from boot.sim STARRs
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=STARR(x,alpha)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the STARR of the boot.run standard errors
  return(mean(res))
}

