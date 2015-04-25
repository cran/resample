\name{Quantile}
\alias{Quantile}
\title{
Compute quantiles using type = 6
}
\description{
Front end to quantile, using type = 6 (appropriate for resampling)
}
\usage{
Quantile(x, ..., type = 6)
}
\arguments{
  \item{x}{
    \code{resample} object,
    numerical object, or other object with a method for
    \code{\link[stats]{quantile}}.
}
  \item{\dots}{
    Other arguments passed to \code{\link[stats]{quantile}}.
}
  \item{type}{
    With \code{type=6} and 99 observations, the k\% quantile
    is the k'th smallest observation; this corresponds to equal
    probability above the largest observation, below the smallest
    observation, and between each pair of adjacent observations.
}
}
\details{
  This is a front end to \code{\link[stats]{quantile}}.
}
\value{
  A vector or matrix of quantiles.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\seealso{
\code{\link[stats]{quantile}}
}
\examples{
quantile(1:9, .2)
Quantile(1:9, .2)
}
\keyword{univar}
