#' Wrapper function that computes LPM and the standard error of the estimate
#'
#' @param data data
#' @param \dots any other passthru parameters. This include two types of parameters.
#' The first type is parameters associated with the risk/performance measure, such as tail
#' probability for VaR and ES. The second type is the parameters associated with the metohd
#' used to compute the standard error. See \code{\link{SE.IF.iid}}, \code{\link{SE.IF.cor}},
#' \code{\link{SE.BOOT.iid}}, \code{\link{SE.BOOT.cor}} for details.
#' @param const constant threshold
#' @param k order of the LPM, can only be 1 or 2 if se.method is not "none"
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}.
#'
#' @return a vector or a list depending on se.method
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}

LPM.SE = function(data, ..., const = 0, k = 1, se.method = "none"){
  data = checkData(data)
  myLPM = apply(data, 2, LPM, const = const, k = k, ...)
  names(myLPM) = colnames(data)
  if(se.method[1] == "none" & length(se.method)==1){
    return(myLPM)
  } else {
    res=list(LPM=myLPM)
    # for each of the method specified in se.method, compute the standard error
    for(mymethod in se.method){
      res[[mymethod]]=EstimatorSE(data, estimator.fun = "LPM",
                                  se.method = mymethod, const = const, k = k, ...)
    }
    return(res)
  }
}
