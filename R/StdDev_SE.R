###############################################################################
# $Id$
###############################################################################





#' calculates Standard Deviation for univariate and multivariate series, also
#' calculates component contribution to standard deviation of a portfolio
#'
#' calculates Standard Deviation for univariate and multivariate series, also
#' calculates component contribution to standard deviation of a portfolio
#'
#' TODO add more details
#'
#' This wrapper function provides fast matrix calculations for univariate,
#' multivariate, and component contributions to Standard Deviation.
#'
#' It is likely that the only one that requires much description is the
#' component decomposition.  This provides a weighted decomposition of the
#' contribution each portfolio element makes to the univariate standard
#' deviation of the whole portfolio.
#'
#' Formally, this is the partial derivative of each univariate standard
#' deviation with respect to the weights.
#'
#' As with \code{\link{VaR}}, this contribution is presented in two forms, both
#' a scalar form that adds up to the univariate standard deviation of the
#' portfolio, and a percentage contribution, which adds up to 100%.  Note that
#' as with any contribution calculation, contribution can be negative.  This
#' indicates that the asset in question is a diversified to the overall
#' standard deviation of the portfolio, and increasing its weight in relation
#' to the rest of the portfolio would decrease the overall portfolio standard
#' deviation.
#'
#' @param R a vector, matrix, data frame, timeSeries or zoo object of asset
#' returns
#' @param \dots any other passthru parameters
#' @param clean method for data cleaning through \code{\link{Return.clean}}.
#' Current options are "none", "boudt", or "geltner".
#' @param portfolio_method one of "single","component" defining whether to do
#' univariate/multivariate or component calc, see Details.
#' @param weights portfolio weighting vector, default NULL, see Details
#' @param mu If univariate, mu is the mean of the series. Otherwise mu is the
#' vector of means of the return series , default NULL, , see Details
#' @param sigma If univariate, sigma is the variance of the series. Otherwise
#' sigma is the covariance matrix of the return series , default NULL, see
#' Details
#' @param use an optional character string giving a method for computing
#' covariances in the presence of missing values.  This must be (an
#' abbreviation of) one of the strings \code{"everything"}, \code{"all.obs"},
#' \code{"complete.obs"}, \code{"na.or.complete"}, or
#' \code{"pairwise.complete.obs"}.
#' @param method a character string indicating which correlation coefficient
#' (or covariance) is to be computed.  One of \code{"pearson"} (default),
#' \code{"kendall"}, or \code{"spearman"}, can be abbreviated.
#' @param se.method a character string indicating which method should be used to compute
#' the standard error of the estimated standard deviation. One of \code{"none"} (default),
#' \code{"IFiid"}, \code{"IFcor"}, \code{"BOOTiid"}, \code{"BOOTcor"}
#' @author Brian G. Peterson and Kris Boudt
#' @seealso \code{\link{Return.clean}} \code{sd}
###keywords ts multivariate distribution models

StdDev.se=function(R , ..., clean=c("none","boudt","geltner"),  portfolio_method=c("single","component"), weights=NULL, mu=NULL, sigma=NULL, use="everything", method=c("pearson", "kendall", "spearman"),
                    se.method=c("none","IFiid","IFcor","BOOTiid","BOOTcor")){

}
