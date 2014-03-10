% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{limits.percentile}
\alias{limits.percentile}
\alias{limits.t}
\title{
Bootstrap confidence intervals
}
\description{
Bootstrap confidence intervals - percentile method or t interval.
}
\usage{
limits.percentile(x, probs = c(0.025, 0.975), ...)
limits.t(x, probs = c(0.025, 0.975))
}
\arguments{
  \item{x}{
    a \code{\link{bootstrap}} or \code{\link{bootstrap}} object.
}
  \item{probs}{
    probability values, between 0 and 1. The default vector
    \code{c(0.025, 0.975)} gives a 95\% two-sided interval.
}
  \item{\dots}{
    additional arguments to pass to \code{\link{quantile.resample}}
    and \code{\link[stats]{quantile}}.
}
}
\value{
  a matrix with one row for each value in \code{probs} and one column
  for each statistic.
}
\author{
  Tim Hesterberg <timhesterberg@gmail.com>
}
\seealso{
\code{\link{bootstrap}},
\code{\link{bootstrap2}}.
}
\examples{
mydata <- data.frame(a = runif(40), b = rnorm(40))
boot1 <- bootstrap(mydata, colMeans)
limits.t(boot1)
limits.percentile(boot1)
}
\keyword{nonparametric}
\keyword{htest}
