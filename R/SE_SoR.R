#' Wrapper function to compute the standard error for Sortino Ratio estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated Sharpe Ratio. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of SoR
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SoR(edhec,se.method="IFiid")
SE.SoR = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.SoR.iid.xts(data,...),
               IFcor = SE.SoR.cor.xts(data,...),
               BOOTiid = SE.SoR.boot.iid.xts(data,...),
               BOOTcor = SE.SoR.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the Sharpe Raio for an xts object
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
#' SE.SoR.iid.xts(edhec)
SE.SoR.iid.xts = function(x,...,rf=0){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SoR.iid(x,rf=rf))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SoR.iid,rf=rf))
  }
}


#' Standard Error of SoR for iid data
#'
#' @param data vector of data
#'
#' @return SE of SoR
#' @export
#'
#' @examples
#' SE.SoR.iid(rnorm(10))
SE.SoR.iid = function(data,rf=0){
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
  return(sqrt(V/length(data)))
}

#' Compute standard error of SoR via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#'
#' @return SE of SoR
#' @export
#'
#' @examples
#' SE.SoR.boot.iid(rnorm(10))
SE.SoR.boot.iid = function(data, ..., rf = 0){
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
    # each standard error is computed from boot.sim SoRs
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=SoR(x,rf=rf)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the mean of the boot.run standard errors
  return(mean(res))
}

