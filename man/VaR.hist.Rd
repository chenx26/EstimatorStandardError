% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Measures.R
\name{VaR.hist}
\alias{VaR.hist}
\title{Compute sample value-at-risk}
\usage{
VaR.hist(data, ..., alpha = 0.05)
}
\arguments{
\item{data}{Vector of data}

\item{alpha}{Tail probability}
}
\value{
Sample VaR
}
\description{
Compute sample value-at-risk
}
\examples{
VaR.hist(rnorm(10))
}
\author{
Xin Chen, \email{chenx26@uw.edu}
}
