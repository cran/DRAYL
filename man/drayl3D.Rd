\name{drayl3D}
\alias{drayl3D}
\title{Three dimensional Rayleigh density by series expansion}
\usage{
drayl3D(dK,Ccomp,lim)
}
\arguments{
  \item{dK}{Determinant of the covariance matrix.}
  \item{Ccomp}{"Compressed" cofactor matrix, leaving out zero value entries.}
  \item{lim}{Number of series terms.}
}
\description{
Returns a 3D Rayleigh density for arbitrary covariance values. The resulting function can then be evaluated at arbitrary points.
}
\value{The 3D Rayleigh density for the compressed cofactor matrix \code{Ccomp} of the covariance matrix. 
The function can then be evaluated for 3-dimensional vectors \code{r}.}
\examples{
library("RConics")

# Matrix
K3 = matrix(0,nrow = 6,ncol = 6)
sigma3 = sqrt(c(0.5,1,1.5))
diag(K3) = c(0.5,0.5,1,1,1.5,1.5)

# rho_12 rho_13 rho_23
rho3<-c(0.9,0.8,0.7)

K3[1,3]=K3[3,1]=K3[2,4]=K3[4,2]=sigma3[1]*sigma3[2]*rho3[1]
K3[1,5]=K3[5,1]=K3[2,6]=K3[6,2]=sigma3[1]*sigma3[3]*rho3[2]
K3[3,5]=K3[5,3]=K3[4,6]=K3[6,4]=sigma3[2]*sigma3[3]*rho3[3]


C3=adjoint(K3)
n = nrow(K3)/2
Ccomp3<-C3[seq(1,(2*n-1),2),][,seq(1,(2*n-1),2)]
dK3<-det(K3)


pdf3D<-drayl3D(dK = dK3, Ccomp = Ccomp3, lim = 3)

pdf3D(rep(1,3))
}
