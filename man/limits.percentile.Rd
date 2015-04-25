% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{CI.percentile}
\alias{CI.percentile}
\alias{CI.t}
\alias{CI.bca}
\alias{CI.bootstrapT}
\title{
Bootstrap confidence intervals
}
\description{
Bootstrap confidence intervals - percentile method or t interval.
}
\usage{
CI.percentile(x, probs = c(0.025, 0.975), expand = TRUE, ...)
CI.t(x, probs = c(0.025, 0.975))
CI.bca(x, probs = c(0.025, 0.975), expand = TRUE, L = NULL, ...)
CI.bootstrapT(x, probs = c(0.025, 0.975))
}
\arguments{
  \item{x}{
    a \code{\link{bootstrap}} or \code{\link{bootstrap}} object.
}
  \item{probs}{
    probability values, between 0 and 1. The default vector
    \code{c(0.025, 0.975)} gives a 95\% two-sided interval.
}
  \item{expand}{
    logical, if \code{TRUE} then use modified percentiles for better
    small-sample accuracy.
}
  \item{\dots}{
    additional arguments to pass to \code{\link{quantile.resample}}
    and \code{\link[stats]{quantile}}.
}
  \item{L}{vector of length \code{n}, empirical influence function
    values. If not supplied this is computed using \code{\link{jackknife}}.
}
}
\value{
  a matrix with one column for each value in \code{probs} and one row
  for each statistic.
}
\details{
  \code{CI.bootstrapT} assumes the first dimension of the statistic
  is an estimate, and the second is proportional to a SE for the
  estimate. E.g. for bootstrapping the mean, they could be the mean and s.
  This is subject to change.

  \code{CI.bca} and \code{CI.bootstrapT} currently only support
  a single sample.
}
\references{
  This discusses the expanded percentile interval:
  Hesterberg, Tim (2014), What Teachers Should Know about the
  Bootstrap: Resampling in the Undergraduate Statistics Curriculum,
  \url{http://arxiv.org/abs/1411.5279}.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\seealso{
\code{\link{bootstrap}},
\code{\link{bootstrap2}},
\code{\link{ExpandProbs}} (for the expanded intervals).
}
\examples{
\dontrun{
# See full set of examples in resample-package, including different
# ways to call all four functions depending on the structure of the data.
data(Verizon)
CLEC <- with(Verizon, Time[Group == "CLEC"])
bootC <- bootstrap(CLEC, mean, seed = 0)
bootC2 <- bootstrap(CLEC, c(mean = mean(CLEC), sd = sd(CLEC)), seed = 0)
CI.percentile(bootC)
CI.t(bootC)
CI.bca(bootC)
CI.bootstrapT(bootC2)
}
}
\keyword{nonparametric}
\keyword{htest}
