\name{drayl_int4D}
\alias{drayl_int4D}
\title{Four Dimensional Rayleigh Density by Integration}
\usage{
drayl_int4D(r,omega,sigma,cor,method)
}
\arguments{
  \item{r}{Evaluation point.}
  \item{omega}{Omega construct necessary for the Integration method.}
  \item{sigma}{Variances of the signals.}
  \item{cor}{Correlation structure.}
  \item{method}{Integration methods, either "Romberg","Cubature" or "Quadrature".}
}
\description{
A four dimensional Rayleigh density by integration.
}
\value{Evaluates the 4D Rayleigh density at the point \code{r}, for the values 
\code{omega},\code{sigma} and \code{cor} as specified by Bealieu's method.}
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

sigma4 = c(sqrt(c(K4[1,1],K4[3,3],K4[5,5],K4[7,7])))

cor4 = c(K4[1,3]/(sigma4[1]*sigma4[2]),
         K4[1,5]/(sigma4[1]*sigma4[3]),
         K4[1,7]/(sigma4[1]*sigma4[4]),
         K4[3,5]/(sigma4[2]*sigma4[3]),
         K4[3,7]/(sigma4[2]*sigma4[4]),
         K4[5,7]/(sigma4[3]*sigma4[4]))

omega4=omega4<-matrix(data = c(1,cor4[1],cor4[2],cor4[3],cor4[1],1,cor4[4],
                      cor4[5],cor4[2],cor4[4],1,cor4[6],cor4[3],cor4[5],cor4[6],1),nrow = 4)

drayl_int4D(c(1,1,1,1),omega = omega4,sigma = sigma4,cor = cor4, method = "Cubature")

}
