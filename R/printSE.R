#' Formatted output of the results from xxx.SE() functions
#'
#' @param res the results from the xx.SE() functions
#' @param round.digit the number of digits to round to
#'
#' @return No return values
#' @export
#'
#' @examples
#' data(edhec)
#' res = ES.SE(edhec, p=.95, method="historical",se.method = "IFiid")
#' printSE(res, round.digit = 4)
printSE = function(res , round.digit = 3){
  N = length(res)
  if(N != 2) {
    cat("the results do not contain standard errors!\n")
    return()
  }
  list.names = names(res)
  res.df = data.frame(as.vector(res[[1]]), as.vector(res[[2]]))
  colnames(res.df) = list.names
  rownames(res.df) = colnames(res[[1]])
  res.df = round(res.df, digits = round.digit)
  res.df[2] = paste("(",res.df[,2],")")
  print(res.df)
}
