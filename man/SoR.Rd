% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Measures.R
\name{SoR}
\alias{SoR}
\title{Compute sample Sortino Ratio with mean threshold}
\usage{
SoR(data, ..., rf = 0)
}
\arguments{
\item{data}{vector of data}

\item{rf}{risk-free rate}
}
\value{
Sample SoR
}
\description{
Compute sample Sortino Ratio with mean threshold
}
\examples{
SoR(rnorm(10))
}
\author{
Xin Chen, \email{chenx26@uw.edu}
}
