#' Wrapper function to compute the standard error for standard deviation estimator
#'
#' @param data xts object
#' @param ... extra arugements to be passed to lower level functions, see details
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}. Currely, only \code{"IFiid"}
#' is implemented.
#'
#' @return SE of SD
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SD(edhec,se.method="IFiid")
SE.SD = function(data,...,se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){
  se.method=se.method[1]
  res = switch(se.method,
               none = NULL,
               IFiid = SE.SD.iid.xts(data),
               IFcor = SE.SD.cor.xts(data,...),
               BOOTiid = SE.SD.boot.iid.xts(data,...),
               BOOTcor = SE.SD.boot.cor.xts(data,...)
  )
  return(res)
}

#' Wrapper function to compute the standard error(s) of the standard deviation for an xts object
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
#' SE.SD.iid.xts(edhec)
SE.SD.iid.xts = function(x){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SD.iid(x))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SD.iid))
  }
}


#' Standard Error of SD for iid data
#'
#' @param data vector of inputs
#'
#' @return SE of SD
#' @export
#'
#' @examples
#' SE.SD.iid(rnorm(10))
SE.SD.iid = function(data){
  N=length(data)
  return(sqrt(var(data)/2/sd(data)/N))
}

#' Compute standard error of SD via bootstrapping for iid data
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#'
#' @return SE of SD
#' @export
#'
#' @examples
#' SE.SD.boot.iid(rnorm(10))
SE.SD.boot.iid = function(data, ...){
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
    # each standard error is computed from boot.sim SDs
    for(sim.iter in 1:boot.sim){
      x=sample(data,N,replace = TRUE)
      res.sd[sim.iter]=sd(x)
    }
    res[run.iter]=sd(res.sd)
  }
  # return the mean of the boot.run standard errors
  return(mean(res))
}

