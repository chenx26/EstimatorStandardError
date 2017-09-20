rm(list=ls())
library(EstimatorStandardError)
#--------- Load data
data(edhec)
colnames(edhec)=c("CA", "CTAG", "DIS", "EM",
                  "EMN", "ED", "FIA", "GM", "L/S",
                  "MA", "RV", "SS", "FoF")
nboot=500
#--------- Set random seed
set.seed(123)
#--------- Expected Shortfall estimates and their S.E.'s

res.ES=ES.SE(edhec, p=.95, method="historical",nsim = nboot,
             se.method = c("IFiid","IFcor","BOOTiid","BOOTcor"),
             standardize = FALSE)

res.ES.df = printSE(res.ES, round.digit = 3, valonly = FALSE)

###### simulation test
tmp = xts(matrix(rnorm(1500, 0.05, 0.2), nrow = 100), Sys.Date() - seq(100, 1))
standardize = function(x){
  return(x.standardized = (x - mean(x))/sd(x), sigma = sd(x))
}

tmp.standardized = apply(tmp, 2, function(x) (x - mean(x))/sd(x))
tmp.sigmas = apply(tmp, 2, sd)

res.ES = ES.SE(R = tmp, p=.95, method="historical",nsim = nboot,
      se.method = c("IFiid","IFcor","BOOTiid","BOOTcor"))
printSE(res.ES)

res.ES.standardized = ES.SE(R = tmp.standardized, p=.95, method="historical",nsim = nboot,
                            se.method = c("IFiid","IFcor","BOOTiid","BOOTcor"))
printSE(res.ES.standardized)

apply(edhec, 2, mean)
apply(edhec, 2, sd)

(dd <- data.frame(x = 1:8, f = gl(2,4), ch = I(letters[1:8])))
# print() with defaults
print(dd, quote = TRUE, row.names = FALSE)
# suppresses row.names and quotes all entries

d=7
alpha=0.5
keep=1
data = ES.IF(coredata(edhec[,1]))

res = myperiodogram(tmp.IF)
plot(my.freq, my.periodogram, type = "l")


