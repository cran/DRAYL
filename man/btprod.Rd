\name{btprod}
\alias{btprod}
\title{Auxilliary function computing intermediate products.}
\usage{
btprod(t,a,Jstar)
}
\arguments{
  \item{t}{Index number.}
  \item{a}{The respective Alpha matrix value.}
  \item{Jstar}{Matrix of the j-star indeces of the series expansion.}
}
\value{Returns the row-wise multiplication of the coefficients based on the indeces \code{j}.}
\description{
Auxilliary function. Based on the results of the \code{btcol} the row wise results are computed.
}
