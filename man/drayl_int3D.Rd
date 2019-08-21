\name{drayl_int3D}
\alias{drayl_int3D}
\title{Three Dimensional Rayleigh Density by Integration}
\usage{
drayl_int3D(r,omega,sigma,cor,method)
}
\arguments{
  \item{r}{Evaluation point.}
  \item{omega}{Omega construct necessary for the Integration method.}
  \item{sigma}{Variances of the signals.}
  \item{cor}{Correlation structure.}
  \item{method}{Integration methods, either "Kronrod","Clenshaw","Simpson","Romberg","TOMS614" or "mixed".}
}
\description{
A three dimensional Rayleigh density by integration.
}
\value{Evaluates the 3D Rayleigh density at the point \code{r}, for the values 
\code{omega},\code{sigma} and \code{cor} as specified by Bealieu's method.}
\examples{
# Matrix
K3 = matrix(0,nrow = 6,ncol = 6)
sigma3 = sqrt(c(0.5,1,1.5))
diag(K3) = c(0.5,0.5,1,1,1.5,1.5)

# rho_12 rho_13 rho_23
rho3<-c(0.9,0.8,0.7)

K3[1,3]=K3[3,1]=K3[2,4]=K3[4,2]=sigma3[1]*sigma3[2]*rho3[1]
K3[1,5]=K3[5,1]=K3[2,6]=K3[6,2]=sigma3[1]*sigma3[3]*rho3[2]
K3[3,5]=K3[5,3]=K3[4,6]=K3[6,4]=sigma3[2]*sigma3[3]*rho3[3]

cor3 = rho3

mat<-diag(3)
mat[1,2]=mat[2,1]=cor3[1]
mat[1,3]=mat[3,1]=cor3[2]
mat[2,3]=mat[3,2]=cor3[3]

omega3=mat

drayl_int3D(c(1,1,1),omega = omega3,sigma = sigma3,cor = cor3, method = "Romberg")

}
