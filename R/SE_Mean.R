#' Wrapper function to compute the standard error for mean estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of mean
#' @export
#'
#' @examples
#' data(edhec)
#' SE.Mean(edhec,se.method="IFiid")
SE.Mean = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.Mean.iid.xts(data),
               IFcor = SE.Mean.cor.xts(data,...),
               BOOTiid = SE.Mean.boot.iid.xts(data,...),
               BOOTcor = SE.Mean.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the mean for an xts object
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
#' SE.Mean.iid.xts(edhec)
SE.Mean.iid.xts = function(x){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
#    if(na.rm) x <- na.omit(x)
    return(SE.Mean.iid(x))
  }
  else {
    x <- coredata(x)
#    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.Mean.iid))
  }
}


#' Standard Error of Mean for iid data
#'
#' @param data vector of data
#'
#' @return SE of Mean
#' @export
#'
#' @examples
#' SE.Mean.iid(rnorm(10))
SE.Mean.iid = function(data){
  N=length(data)
  return(sqrt(var(data)/N))
}

#' Compute standard error of mean via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#'
#' @return SE of mean
#' @export
#'
#' @examples
#' SE.Mean.boot.iid(rnorm(10))
SE.Mean.boot.iid = function(data, ...){
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
    # each standard error is computed from boot.sim means
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=mean(x)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the mean of the boot.run standard errors
  return(mean(res))
}

