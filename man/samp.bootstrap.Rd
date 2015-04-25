% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{samp.bootstrap}
\alias{samp.bootstrap}
\alias{samp.permute}
\title{
  Generate indices for resampling
}
\description{
  Generate indices for resampling.
}
\usage{
samp.bootstrap(n, R, size = n - reduceSize, reduceSize = 0)
samp.permute(n, R, size = n - reduceSize, reduceSize = 0,
             groupSizes = NULL, returnGroup = NULL)
}
\arguments{
  \item{n}{
    sample size. For two-sample permutation tests, this is the sum of
    the two sample sizes.
}
  \item{R}{
    number of vectors of indices to produce.
}
  \item{size}{
    size of samples to produce. For example, to do "what-if" analyses,
    to estimate the variability of a statistic had the data been a
    different size, you may specify the size.
}
  \item{reduceSize}{
    integer; if specified, then \code{size = n - reduceSize}
    (for each sample or stratum).
    This is an alternate way to specify size. Typically bootstrap standard
    errors are too small; they correspond to using \code{n} in the
    divisor of the sample variance, rather than \code{n-1}.
    By specifying \code{reduceSize = 1}, you can correct for that bias.
    This is particularly convenient in two-sample problems where the
    sample sizes differ.
}
  \item{groupSizes}{
    \code{NULL}, or vector of positive integers that add to \code{n}.
}
  \item{returnGroup}{
    \code{NULL}, or integer from 1 to \code{length(groupSizes)}.
    \code{groupSizes} and \code{returnGroup} must be supplied together;
    then full permutations are created, but only subsets of size
    \code{groupSizes[returnGroup]} is returned.
}
}
\value{
  matrix with \code{size} rows and \code{R} columns
  (or \code{groupSizes(returnGroup)} rows).
  Each column contains indices for one bootstrap sample, or one permutation.
}
\details{
  To obtain disjoint samples without replacement,
  call this function multiple times, after setting the same random
  number seed, with the same \code{groupSizes} but different values of
  \code{returnGroup}. This is used for two-sample permutation tests.

  If \code{groupSizes} is supplied then \code{size} is ignored.
}
\note{
  The value passed as \code{R} to this function is typically the
  \code{block.size} argument to \code{\link{bootstrap}} and other
  resampling functions.
}
\references{
  This discusses reduced sample size:
  Hesterberg, Tim C. (2004), Unbiasing the Bootstrap-Bootknife Sampling
  vs. Smoothing, Proceedings of the Section on Statistics and the
  Environment, American Statistical Association, 2924-2930,
  \url{http://www.timhesterberg.net/articles/JSM04-bootknife.pdf}.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\seealso{
  \code{\link{resample-package}}.
}
\examples{
samp.bootstrap(7, 8)
samp.bootstrap(7, 8, size = 6)
samp.bootstrap(7, 8, reduceSize = 1)

# Full permutations
set.seed(0)
samp.permute(7, 8)

# Disjoint samples without replacement = subsets of permutations
set.seed(0)
samp.permute(7, 8, groupSizes = c(2, 5), returnGroup = 1)
set.seed(0)
samp.permute(7, 8, groupSizes = c(2, 5), returnGroup = 2)
}
\keyword{nonparametric}
\keyword{htest}
