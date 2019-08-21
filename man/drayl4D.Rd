\name{drayl4D}
\alias{drayl4D}
\title{Four dimensional Rayleigh density by series expansion}
\usage{
drayl4D(dK,Ccomp,lim)
}
\arguments{
  \item{dK}{Determinant of the covariance matrix.}
  \item{Ccomp}{"Compressed" cofactor matrix, leaving out zero value entries.}
  \item{lim}{Number of series terms.}
}
\description{
Returns a 4D Rayleigh density for arbitrary covariance values. The resulting function can then be evaluated at arbitrary points.
}
\value{The 4D Rayleigh density for the compressed cofactor matrix \code{Ccomp} of the covariance matrix. 
The function can then be evaluated for 4-dimensional vectors \code{r}.}
\examples{
library("RConics")

K4 = matrix(0,nrow = 8,ncol = 8)
sigma4 = sqrt(c(0.5,1,1.5,1))
rho4<-c(0.7,0.75,0.8,0.7,0.75,0.7)

K4[1,1]=K4[2,2]=sigma4[1]^2
K4[3,3]=K4[4,4]=sigma4[2]^2
K4[5,5]=K4[6,6]=sigma4[3]^2
K4[7,7]=K4[8,8]=sigma4[4]^2

K4[1,3]=K4[3,1]=K4[2,4]=K4[4,2]=sigma4[1]*sigma4[2]*rho4[1]
K4[1,5]=K4[5,1]=K4[2,6]=K4[6,2]=sigma4[1]*sigma4[3]*rho4[2]
K4[1,7]=K4[7,1]=K4[2,8]=K4[8,2]=sigma4[1]*sigma4[4]*rho4[3]
K4[3,5]=K4[5,3]=K4[4,6]=K4[6,4]=sigma4[2]*sigma4[3]*rho4[4]
K4[3,7]=K4[7,3]=K4[4,8]=K4[8,4]=sigma4[2]*sigma4[4]*rho4[5]
K4[5,7]=K4[7,5]=K4[8,6]=K4[6,8]=sigma4[3]*sigma4[4]*rho4[6]

C4=adjoint(K4)
n = nrow(K4)/2
Ccomp4<-C4[seq(1,(2*n-1),2),][,seq(1,(2*n-1),2)]
dK4<-det(K4)

pdf4D<-drayl4D(dK = dK4, Ccomp = Ccomp4, lim = 3)
pdf4D(rep(1,4))
}
