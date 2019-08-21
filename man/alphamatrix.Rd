\name{alphamatrix}
\alias{alphamatrix}
\title{Computation of Alpha coefficient matrix}
\usage{
alphamatrix(n)
}
\arguments{
  \item{n}{Distribution dimension.}
}
\value{Returns a \code{n-1} dimensional matrix that contains the permutations of
all indeces.}
\description{
The alpha matrix is a necessary intermediate step in the series expansion approach.
It lists the different parameter combinations necessary for the series expansion.
}
\examples{
alphamatrix(3)
}
