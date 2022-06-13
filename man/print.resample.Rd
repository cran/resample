% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{print.resample}
\alias{plot.resample}
\alias{hist.resample}
\alias{qqnorm.resample}
\alias{print.resample}
\alias{quantile.resample}
\title{
  Methods for common generic functions for resample objects
}
\description{
  Methods for common generic functions. The methods
  operate primarily on the replicates
  (resampled statistics).
}
\usage{
\method{print}{resample}(x, \dots)
\method{hist}{resample}(x, \dots, resampleColumns = 1:x$p, xlim = NULL,
                xlab = NULL, main = "", col = "blue", border = 0,
                breaks = "FD", showObserved = TRUE,
                legend = TRUE, args.legend = NULL)
\method{plot}{resample}(x, \dots)
\method{qqnorm}{resample}(y, \dots, resampleColumns = 1:y$p, ylab = NULL,
                pch = if(y$R < 100) 1 else ".")
\method{quantile}{resample}(x, \dots)
}
\arguments{
  \item{x,y}{
    a \code{"resample"} object, usually produced by one of
    \code{\link{bootstrap}},
    \code{\link{bootstrap2}},
    \code{\link{permutationTest}}, or
    \code{\link{permutationTest2}}.
}
  \item{\dots}{
    additional arguments passed to the corresponding generic function.\cr
    For \code{plot.resample}, these are passed to \code{hist.resample}.
}
  \item{resampleColumns}{
    integer subscripts, or names of statistics.
    When a statistic is a vector, resampleColumns may be used to select
    which resampling distributions to plot.
}
  \item{xlim}{
    limits for the x axis.
}
  \item{xlab, ylab}{
    x and y axis labels.
}
  \item{main}{
    main title
}
  \item{col}{
    color used to fill bars, see \code{\link[graphics]{hist}}.
}
  \item{border}{
    color of the order around the bars, see \code{\link[graphics]{hist}}.
}
  \item{breaks}{
    method for computing breaks, see \code{\link[graphics]{hist}}.
}
  \item{showObserved}{
    logical, if \code{TRUE} then vertical lines are shown at the
    observed statistic and mean of the bootstrap replicates.
}
  \item{legend}{
    logical, if \code{TRUE} a legend is added. Not used if
    \code{showObserved = FALSE}.
}
  \item{args.legend}{
    \code{NULL} or a list of arguments to pass to
    \code{\link[graphics]{legend}}.
}
  \item{pch}{
    plotting character, see \code{\link[graphics]{par}}.
}
}
\value{
  For \code{quantile.resample}, a matrix with one row for each
  statistic and one column for each value in \code{probs}.
  This uses \code{type=6} when calling
  \code{\link[stats]{quantile}}, for wider (more accurate) quantiles than
  the usual default.

The other functions are not called for their return values.
}
\details{
  \code{hist.resample} displays a histogram overlaid with a density
  plot, with the observed value of the statistic indicated.

  \code{plot.resample} currently just calls \code{hist.resample}.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
  \code{\link{resample-package}},
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{jackknife}},
  \code{\link{permutationTest}},\cr
  \code{\link{permutationTest2}},
  \code{\link[stats]{quantile}}.
}
\examples{
\donttest{
# See full set of examples in resample-package
data(Verizon)
CLEC <- with(Verizon, Time[Group == "CLEC"])
bootC <- bootstrap(CLEC, mean, seed = 0)
print(bootC)
hist(bootC)
qqnorm(bootC)
quantile(bootC, probs = c(.25, .975))
# That is the percentile interval with expand = FALSE
CI.percentile(bootC)
}
}
\keyword{nonparametric}
\keyword{htest}
