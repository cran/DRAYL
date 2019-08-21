\name{zerooneoutput}
\alias{zerooneoutput}
\title{Non-zero value determination}
\usage{
zerooneoutput(j,rho,A)
}
\arguments{
  \item{j}{Vector of j indeces.}
  \item{rho}{Vector of the rho index.}
  \item{A}{Alpha matrix.}
}
\description{
Determines the contribution of sum terms, based on the index \code{j}, \code{rho} and the matrix \code{A}.
}
\value{Either 0 or 1, computes the integral contribution based on the alphamatrix \code{A}.}
\examples{
A = alphamatrix(3)
zerooneoutput(c(0,0,0),c(-1,-1,-1),A)
}
