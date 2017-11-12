#' Compute long run variance by fitting AR model with AIC or BIC model selection
#'
#' @param x vector of input data time series
#' @param max.order maximum order for AR model, default 10
#' @param criteria model selection criteria, "AIC" or "BIC"
#' @param all.results if TRUE, return all fitted model along with longrun variance, otherwise return lrvar only
#' @param AR.method method used to fit the AR model, default "ML"
#'
#' @return lrvar or a list with first element lrvar and second element a list of all fitted model
#' @export
#'
#' @examples
#' SE.AR(rnorm(100))
SE.AR = function(x, max.order = 10,
                 criteria = c("AIC", "BIC"),
                 all.results = FALSE,
                 AR.method = c("ML", "CSS-ML", "CSS")){
  orders = seq(0, max.order, by = 1)
  criteria = criteria[1]
  AR.method = AR.method[1]
  res = list()
  for (i in 1:length(orders)){
    res[[i]] = arima(x, order = c(orders[i], 0, 0), method = AR.method)
  }
  if(criteria == "AIC"){
    IC_vector = sapply(res, AIC)
  }
  if(criteria == "BIC"){
    IC_vector = sapply(res, BIC)
  }
  best_AR = res[[which.min(IC_vector)]]
  ar.coeffs = best_AR$coef
  epsilon_variance = best_AR$sigma2
  variance = epsilon_variance/(1 - sum(ar.coeffs))^2
  if(all.results){
    return(list(variance = variance/length(x), results = res))
  }
  return(variance/length(x))
}
