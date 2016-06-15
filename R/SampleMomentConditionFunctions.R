#### Collection of Sample Moment Condition Function computations


#' Compute the sample moment condition functions of SR from data
#'
#' @param data Vector of data
#'
#' @return 2 by N matrix of Sample Moment Condition Functions
#' @export
#'
#' @examples
#' SR.SMCF(rnorm(10))
SR.SMCF=function(data){
  mu.hat=mean(data)
  sd.hat=sd(data)
  ndata=length(data)
  psi.mat=matrix(0,2,ndata)
  for(iter in 1:ndata){
    psi.mat[1,iter]=data[iter]-mu.hat
    psi.mat[2,iter]=(data[iter]-mu.hat)^2-sd.hat^2
  }
  return(psi.mat)
}