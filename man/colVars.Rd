\name{colVars}
\alias{colVars}
\alias{colStdevs}
\title{
Column variances and standard deviations for matrices.
}
\description{
Quick and dirty function for column variances and standard deviations.
}
\usage{
colVars(x, na.rm = FALSE)
colStdevs(x, \dots)
}

\arguments{
  \item{x}{
    data frame, matrix, or vector. These versions do not support
  higher-dimensional arrays.
  }
  \item{na.rm}{
    logical. Should missing values (including \code{NaN})
    be omitted from the calculations?
  }
  \item{\dots}{
    other arguments passed to \code{colVars}.
  }
}
\value{
  A numeric or complex array of suitable size, or a vector if the result is
  one-dimensional.  The \code{dimnames} (or \code{names} for a vector
  result) are taken from the original array.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
  \code{\link[base]{colSums}}, \code{\link[stats]{var}}, \code{\link[stats]{sd}}.
}
\note{
  There are better versions of these functions in the aggregate package\cr
  \url{https://www.timhesterberg.net/r-packages}.
  }
\examples{
x <- matrix(rnorm(12), 4)
colVars(x)
colStdevs(x)
}
\keyword{algebra}
\keyword{arith}
