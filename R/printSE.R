#' Formatted output of the results from xxx.SE() functions
#'
#' @param res the results from the xx.SE() functions
#' @param valonly if true, return only the data frame, otherwise print and return data frame
#' @param round.digit the number of digits to round to
#'
#' @return data frame
#' @export
#' @author Xin Chen, \email{chenx26@uw.edu}
#'
#' @examples
#' data(edhec)
#' res = ES.SE(edhec, p=.95, method="historical",se.method = "IFiid")
#' printSE(res, round.digit = 4)
printSE = function(res , round.digit = 3, valonly = TRUE){
  N = length(res)
  # if(N != 2) {
  #   cat("the results do not contain standard errors!\n")
  #   return()
  # }
  list.names = names(res)
  res.df = data.frame(t(res[[1]]))
  for(i in 2:length(list.names)){
    res.df = cbind(res.df, res[[i]])
  }
  colnames(res.df) = list.names
  rownames(res.df) = colnames(res[[1]])
#  res.df = round(res.df, digits = round.digit)
  if(valonly){
    return(res.df)
  }
  print(res.df, digits = round.digit)
  # res.df[2] = paste("(",res.df[,2],")",sep="")
  # res.df[,-1] = apply(as.data.frame(res.df[,-1]),2,function(x) paste("(",x,")",sep=""))

}
