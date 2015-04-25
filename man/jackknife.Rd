\name{jackknife}
\alias{jackknife}
\title{
One sample jackknife
}
\description{
 Basic resampling. Supply the data and statistic to resample.
}
\usage{
jackknife(data, statistic, args.stat = NULL,
          label = NULL, statisticNames = NULL, trace = FALSE)
}
\arguments{
  \item{data}{
    vector, matrix, or data frame.
}
  \item{statistic}{
    a function, or expression (e.g. \code{mean(myData, trim = .2)}.
}
  \item{args.stat}{
    a list of additional arguments to pass to \code{statistic}, if
    it is a function.
}
  \item{label}{
    used for labeling plots (in a future version).
}
  \item{statisticNames}{
    a character vector the same length as the vector returned by
    \code{statistic}.
}
  \item{trace}{
    logical, if \code{TRUE} an indication of progress is printed.
}
}
\value{
  a list with class \code{"jackknife"}
  that inherits from \code{"resample"},
  with components:
  \item{observed}{the value of the statistic for the original data.}
  \item{replicates}{a matrix with \code{R} rows and \code{p} columns.}
  \item{n}{number of observations in the original data, or vector of
  length 2 in two-sample problems.}
  \item{p}{\code{length(observed)}.}
  \item{R}{number of replications.}
  \item{seed}{the value of the seed at the start of sampling.}
  \item{call}{the matched call.}
  \item{statistics}{a data frame with \code{p} rows, with columns
    \code{"observed"}, \code{"mean"} (the mean of the replicates),
    and other columns appropriate to resampling; e.g. the bootstrap
    objects have columns \code{"SE"} and \code{"Bias"},
    while the permutation test objects have \code{"Alternative"}
    and \code{"PValue"}.}

  There are functions for printing and plotting these objects,
  in particular \code{print}, \code{plot}, \code{hist}, \code{qqnorm},
  \code{quantile}.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\note{
The current version only handles a single sample.
}
\seealso{
\code{\link{resample-package}}.
}
\examples{
\dontrun{
# See full set of examples in resample-package
data(Verizon)
CLEC <- with(Verizon, Time[Group == "CLEC"])
jackknife(CLEC, mean)
}
}
\keyword{nonparametric}
