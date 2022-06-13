\name{cat0}
\alias{cat0}
\alias{catn}
\alias{cat0n}
\title{
Front end to cat
}
\description{
Call cat, with \code{sep=""} and/or newline at end.
}
\usage{
cat0(...)
cat0n(...)
catn(...)
}
\arguments{
  \item{\dots}{R objects, like for \code{\link[base]{cat}}
}
}
\details{
  \code{cat0} and \code{cat0n} call \code{cat} with \code{sep = ""}.
  \code{catn} and \code{cat0n} print a final newline).
}
\value{
None (invisible \code{NULL}).
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
\code{\link[base]{cat}}, \code{\link[base]{paste0}}.
}
\examples{
cat("Print this")
# That printed without a final newline.
catn("Print this")
cat0n("10,", "000")
}
\keyword{utilities}
