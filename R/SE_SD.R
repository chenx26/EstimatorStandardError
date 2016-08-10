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

#' Wrapper function to compute the standard error(s) of the SD for an xts object
#' using GLM-EN method
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#'
#' @return standard error(s) of the xts object
#' @export

SE.SD.cor.xts = function(x){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SD.cor(x))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SD.cor))
  }
}

#' Compute the standard error of SD estimator for correlated data vector using
#' GLM-EN method
#'
#' @param x the vector for the data
#' @param d the maximum order of the polynomial
#' @param alpha.lasso the weight for lasso vs elastic net
#' @param keep the portion of peridograms to use to fit the polynomial
#'
#' @return the SE of the SD estimtor for the data
#' @export

SE.SD.cor = function(x, d = 5, alpha.lasso = 0.5, keep = 1){
  N=length(x)
  data.IF = SD.IF(x)
  tmp = SE.GLM.LASSO(data.IF, d = d, alpha = alpha.lasso, keep = keep)
  return(sqrt(tmp))
}

#' Wrapper function to compute the standard error(s) of the SD for an xts object
#' using iid bootstrapping
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#' @param ... other parameters
#'
#' @return standard error(s) of the xts object
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SD.boot.iid.xts(edhec)
SE.SD.boot.iid.xts = function(x,...){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SD.boot.iid(x,...))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SD.boot.iid,...=...))
  }
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

#' Wrapper function to compute the standard error(s) of the SD for an xts object
#' using tsboot
#'
#' The operation is performed column wise.
#'
#' @param x the xts object
#' @param ... other parameters
#'
#' @return standard error(s) of the xts object
#' @export
#'
#' @examples
#' data(edhec)
#' SE.SD.boot.cor.xts(edhec)
SE.SD.boot.cor.xts = function(x,...){
  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    #    if(na.rm) x <- na.omit(x)
    return(SE.SD.boot.cor(x,...))
  }
  else {
    x <- coredata(x)
    #    if(na.rm) x <- na.omit(x)
    return(apply(x, 2, SE.SD.boot.cor,...=...))
  }
}

#' Compute standard error of SD via bootstrapping for data (possibily serially correlated)
#'
#' @param data vector of data
#' @param ... parameters passed from upper calls
#' @param segment.length the length of the fixed block for bootstrapping
#'
#' @return SE of SD
#' @export
#'
#' @examples
#' SE.SD.boot.cor(rnorm(100))
SE.SD.boot.cor = function(data, ...,segment.length=length(data)/5){
  args=list(...)
  arg.names=names(args)
  boot.sim=100
  boot.run=100
  N=length(data)
  # check if the parameters are specified by the function call
  if(("boot.sim" %in% arg.names) & (!is.null(args["boot.sim"]))) boot.sim=args["boot.sim"]
  if(("boot.run" %in% arg.names) & (!is.null(args["boot.run"]))) boot.run=args["boot.run"]
  res=rep(0,boot.run)
  # compute boot.run standard errors
  for(run.iter in 1:boot.run){
    boot.res=tsboot(data, sd, R = boot.sim, sim = "fixed", l = segment.length)
    res[run.iter]=sd(boot.res$t)
  }
  # return the ES of the boot.run standard errors
  return(mean(res))
}
