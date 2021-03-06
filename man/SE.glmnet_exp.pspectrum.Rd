% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SE_pspectrum_glmnet_exp.R
\name{SE.glmnet_exp.pspectrum}
\alias{SE.glmnet_exp.pspectrum}
\title{Compute long run variance by using pspectrum to compute the tapered peridograms first then use glmnet_exp}
\usage{
SE.glmnet_exp.pspectrum(data, ..., d = 7, alpha.EN = 0.5, keep = 1,
  standardize = FALSE, return.coeffs = FALSE, prewhiten = FALSE,
  twosided = FALSE)
}
\arguments{
\item{data}{vector of data}

\item{d}{max order of polynomial}

\item{twosided}{}
}
\description{
Compute long run variance by using pspectrum to compute the tapered peridograms first then use glmnet_exp
}
