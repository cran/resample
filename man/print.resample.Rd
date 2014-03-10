% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{plot.resample}
\alias{plot.resample}
\alias{hist.resample}
\alias{qqnorm.resample}
\alias{print.resample}
\alias{quantile.resample}
\title{
  Methods for common generic functions for resample objects
}
\description{
  Methods for common generic functions. The methods for
  \code{\link[graphics]{plot}},
  \code{\link[graphics]{hist}},
  \code{\link[stats]{qqnorm}}, and
  \code{\link[stats]{quantile}} operate primarily on the replicates.
}
\usage{
\method{print}{resample}(x, ...)
\method{plot}{resample}(x, ..., resampleColumns = 1:x$p, xlim = NULL, xlab = NULL, main = "")
\method{hist}{resample}(x, ..., resampleColumns = 1:x$p, xlim = NULL, xlab = NULL, main = "")
\method{qqnorm}{resample}(y, ..., resampleColumns = 1:y$p, ylab = NULL)
\method{quantile}{resample}(x, ..., type = 6)
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
    additional arguments passed to the corresponding generic function.
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
  \item{type}{
    the \code{type} argument for \code{\link[stats]{quantile}}.
    The value used here gives somewhat wider quantiles than the usual
    default, which are appropriate for resampling.
}
}
\value{
For \code{quantile.resample}, a matrix with one column for each
statistic.

The other functions are not called for their return values.
}
\details{
  \code{hist.resample} and \code{plot.resample} are identical; they
  display a histogram overlaid with a density plot, with the observed value
  of the statistic indicated.
}
\author{
  Tim Hesterberg <timhesterberg@gmail.com>
}
\seealso{
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{permutationTest}},
  \code{\link{permutationTest2}}.
}
\examples{
mydata <- data.frame(a = runif(40), b = rnorm(40))
boot1 <- bootstrap(mydata, colMeans)
boot1
par(mfrow = c(2,1))
plot(boot1)
hist(boot1) # same as plot
qqnorm(boot1)
par(mfrow = c(1,1))
quantile(boot1, probs = c(.25, .75))
}
\keyword{nonparametric}
\keyword{htest}
