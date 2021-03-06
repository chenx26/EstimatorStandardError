#' Newey-West NSE estimators.
#' @description Calculate the variance of the mean with the Newey West (1987, 1994) HAC estimator.
#' @description This is a wrapper around \link[sandwich]{lrvar} from the sandwich package.
#' @param x      A numeric vector or matrix.
#' @param prewhite  A bool indicating if the time-serie will be prewhitened before analysis.
#' @return The variance estimator in the univariate case or the variance-covariance matrix estimator in the multivariate case.
#' @references Andrews, Donald WK. "Heteroskedasticity and autocorrelation consistent covariance matrix estimation." Econometrica: Journal of the Econometric Society 59.03 (1991): 817-858.
#' @references Newey, Whitney K., and Kenneth D. West. "A simple, positive semi-definite, heteroskedasticity and autocorrelationconsistent covariance matrix.", Econometrica: Journal of the Econometric Society 55.03 (1987) : 703-708.
#' @references Newey, Whitney K., and Kenneth D. West. "Automatic lag selection in covariance matrix estimation." The Review of Economic Studies 61.4 (1994): 631-653.
#' @references Zeileis, Achim. "Econometric computing with HC and HAC covariance matrix estimators." (2004).
#' @import sandwich
#' @examples
#'n = 1000
#'ar = c(0.9,0.6)
#'mean = c(1,5)
#'sd = c(10,2)
#'
#'Ts1 = as.vector(arima.sim(n = n, list(ar = ar[1]), sd = sd[1]) + mean[1])
#'Ts2 = as.vector(arima.sim(n = n, list(ar = ar[2]), sd = sd[2]) + mean[2])
#'Ts = cbind(Ts1,Ts2)
#'
#'nse.nw(x = Ts1)
#'nse.nw(x = Ts)
#'nse.nw(x = Ts1, prewhite = TRUE)
#'nse.nw(x = Ts, prewhite = TRUE)
#'@export
nse.nw <- function(x,prewhite = FALSE) {
  out = sandwich::lrvar(x = x, type = "Newey-West", prewhite = prewhite, adjust = TRUE)
  out = unname(out)
  return(out)
}
