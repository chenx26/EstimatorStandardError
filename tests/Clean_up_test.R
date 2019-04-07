# library(devtools)
# install_github("chenx26/glmExpENRcpp")
# install_github("chenx26/EstimatorStandardError")

rm(list=ls())
library(EstimatorStandardError)
#--------- Load data
data(edhec)
colnames(edhec)=c("CA", "CTAG", "DIS", "EM",
                  "EMN", "ED", "FIA", "GM", "L/S",
                  "MA", "RV", "SS", "FoF")
nboot=500
return.coeffs = TRUE
d.GLM.EN = 12
alpha.ES = 0.05
alpah.EN = 0.5
#--------- Set random seed
seed = 12345
k_fold_iter = 1000
set.seed(seed)
#--------- Expected Shortfall estimates and their S.E.'s
se.ES=ES.SE(edhec, p = 1 - alpha.ES, method="historical",nsim = nboot,
            se.method = c("IFiid","IFcor","BOOTiid","BOOTcor"),
            standardize = FALSE, return.coeffs = return.coeffs, 
            k_fold_iter = k_fold_iter,
            d.GLM.EN = d.GLM.EN,
            alpha = alpha.EN)


tmp = se.ES$IFcor
SEs = sapply(tmp, function(x) x[[1]])
se.ES$IFcor = SEs

se.ES.df = printSE(se.ES, round.digit  = 3, valonly = TRUE)

ar1 = apply(edhec, 2,
            function(x) {
              res = ar(x, order.max = 1)
              if(length(res$ar) == 0)
                return(0)
              return(res$ar[1])
            }
)

IF.ar1 = apply(edhec, 2,
               function(x) {
                 res = ar(ES.IF(x), order.max = 1)
                 if(length(res$ar) == 0)
                   return(0)
                 return(res$ar[1])
               }
) 

PI.IFcor = (se.ES.df$IFcor / se.ES.df$IFiid - 1) * 100

# use newey-west to compute SE's

NW = NULL
for(i in 1:ncol(edhec)){
  data = coredata(edhec[,i])
  data.IF = ES.IF(data, alpha.ES = alpha.ES)
  NW = c(NW, sqrt(nse.nw(data.IF)))
}
PI.NW = (NW / se.ES.df$IFiid - 1) * 100

se.ES.df = cbind(se.ES.df, PI.IFcor)[c(1,2,3,6,4,5)]
se.ES.df = cbind(se.ES.df, ar1, IF.ar1)

colnames(se.ES.df) = c("ES", "seIidIF", "seCorIF", "PI.seCorIF", "BOOTiid", "BOOTcor", "ar1", "IF.ar1")
knitr::kable(se.ES.df, digits = c(3,3,3,0,3,3,2,2), align = 'c')

