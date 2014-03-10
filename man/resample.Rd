% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{resample}
\alias{resample}
\title{
  Nonparametric resampling
}
\description{
  This function is called by \code{\link{bootstrap}} and other
  resampling functions to actually perform resampling, but may also be
  called directly.
}
\usage{
resample(data, resampleFun, sampler, B = 1000, seed = NULL,
         statisticNames = NULL,
         block.size = 100, trace = FALSE, ..., observedIndices = 1:n)
}
\arguments{
  \item{data}{
    vector, matrix, or data frame.
}
  \item{resampleFun}{
    a function with argument \code{data} and \code{ii}, that calculates
    a statistic of interest for \code{data[ii]} or
    \code{data[ii, , drop=FALSE]}, for a vector or matrix, respectively.
}
  \item{sampler}{
    a function like \code{\link{samp.bootstrap}} or
    \code{\link{samp.permute}}.
}
  \item{B}{
    number of resamples.
}
  \item{seed}{
    old value of .Random.seed, or argument to set.seed.
}
  \item{statisticNames}{
    a character vector the same length as the vector returned by
    \code{statistic}, giving names.
}
  \item{block.size}{
    integer. The \code{B} replicates are done this many at a time.
}
  \item{trace}{
    logical, if \code{TRUE} an indication of progress is printed.
}
  \item{\dots}{
    addition arguments passed to \code{sampler}.
}
  \item{observedIndices}{
    integer vector of indices, used for calculating the observed value.
    When this is called by \code{\link{bootstrap2}} or
    \code{\link{permutationTest2}}, those should be indices
    corresponding to one sample in a merged data set.
}
}
\value{
  an object of class \code{"resample"}; this is
  a list with components:
  \item{observed}{the observed statistic, length \code{p}.}
  \item{replicates}{a matrix with \code{B} rows and \code{p} columns.}
  \item{n}{number of observations}
  \item{p}{the length of the statistic returned by \code{resampleFun}.}
  \item{B}{number of resamples.}
  \item{seed}{the value of \code{seed} when this function is called.}
}
\details{
  This is called by
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{permutationTest}}, and
  \code{\link{permutationTest2}} to actually perform resampling.
  The results are passed back to the calling function, which may
  add additional components and a class, which inherits from
  \code{"resample"}.

  This may also be called directly. In contrast to the other functions,
  where you have flexibility in how you specify the statistic, here
  \code{resampleFun} must be a function.
}
\author{
  Tim Hesterberg <timhesterberg@gmail.com>
}

\seealso{
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{permutationTest}},
  \code{\link{permutationTest2}},
  \code{\link{samp.bootstrap}},
  \code{\link{samp.permute}}.

  For an overview of all functions in the package, see
  \code{\link{resample-package}}.
}
\examples{
mydata <- data.frame(a = runif(40), b = rnorm(40))
boot <- bootstrap(mydata, function(x, i) colMeans(x[i, , drop = FALSE]))
boot
par(mfrow = c(2,1))
plot(boot)
qqnorm(boot)
par(mfrow = c(1,1))
}
\keyword{nonparametric}
\keyword{htest}
