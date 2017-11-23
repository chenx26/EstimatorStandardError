#' Use pspectrum function to compute long run variance
#'
#' @param x data
#' @param ... other parameters
#' @param verbose verbose or not
#' @param plot plot or not
#'
#' @return long run variance of x
#' @export
#'
#' @examples
SE.pspectrum = function(x, ..., verbose = FALSE, plot = FALSE){
  n = length(x)
  res.pspectrum = pspectrum(x, verbose = verbose, plot = plot)
  return(res.pspectrum$spec[1]/n/2)
}
