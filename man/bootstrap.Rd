% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap2}
\alias{permutationTest}
\alias{permutationTest2}
\title{
  One and two sample bootstrap sampling and permutation tests.
}
\description{
  Basic resampling. Supply the data and statistic to resample.
}
\usage{
bootstrap(data, statistic, R = 10000,
          args.stat = NULL, seed = NULL, sampler = samp.bootstrap,
          label = NULL, statisticNames = NULL, block.size = 100,
          trace = FALSE)
bootstrap2(data, statistic, treatment, data2 = NULL, R = 10000,
          ratio = FALSE,
          args.stat = NULL, seed = NULL, sampler = samp.bootstrap,
          label = NULL, statisticNames = NULL, block.size = 100,
          trace = FALSE)
permutationTest(data, statistic, R = 9999,
          alternative = "two.sided", resampleColumns = NULL,
          args.stat = NULL, seed = NULL, sampler = samp.permute,
          label = NULL, statisticNames = NULL, block.size = 100,
          trace = FALSE, tolerance = .Machine$double.eps ^ 0.5)
permutationTest2(data, statistic, treatment, data2 = NULL, R = 9999,
          alternative = "two.sided", ratio = FALSE, paired = FALSE,
          args.stat = NULL, seed = NULL, sampler = samp.permute,
          label = NULL, statisticNames = NULL, block.size = 100,
          trace = FALSE, tolerance = .Machine$double.eps ^ 0.5)
}

\arguments{
  \item{data}{
    vector, matrix, or data frame.
}
  \item{statistic}{
    a function, or expression (e.g. \code{mean(myData, trim = .2)}.
}
  \item{R}{
    number of replicates (bootstrap samples or permutation resamples).
}
  \item{treatment}{
    a vector with two unique values.
    For two-sample applications, suppy either \code{treatment} or
    \code{data2}.
}
  \item{data2}{
    an object like \code{data}; the second sample.
}
  \item{alternative}{
    one of \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
    If \code{statistic} returns a vector, this may be a vector
    of the same length.
}
  \item{ratio}{
    logical, if \code{FALSE} then statistics for two samples are combined using
    statistic1 - statistic2 (the statistics from the two samples).
    If \code{TRUE}, it uses statistic1 / statistic2.
}
  \item{resampleColumns}{
    integer, or character (a subset of the column names of \code{data});
    if supplied then only these columns of the data are permuted.
    For example, for a permutation test of the correlation of x and y,
    only one of the variables should be permuted.
}
  \item{args.stat}{
    a list of additional arguments to pass to \code{statistic}, if
    it is a function.
}
  \item{paired}{
    logical, if \code{TRUE} then observations in \code{data} and
    \code{data2} are paired, and permutations are done within each pair.
    Not yet implemented.
}
  \item{seed}{
    old value of .Random.seed, or argument to set.seed.
}
  \item{sampler}{
    a function for resampling, see \code{help(\link{samp.bootstrap})}.
}
  \item{label}{
    used for labeling plots (in a future version).
}
  \item{statisticNames}{
    a character vector the same length as the vector returned by
    \code{statistic}.
}
  \item{block.size}{
    integer. The \code{R} replicates are done this many at a time.
}
  \item{trace}{
    logical, if \code{TRUE} an indication of progress is printed.
}
  \item{tolerance}{
    when computing P-values, differences smaller than \code{tolerance}
    (absolute or relative)
    between the observed value and the replicates are considered equal.
}
}
\details{
  There is considerable flexibility in how you specify the data and statistic.

  For the \code{statistic}, you may supply a function, or an expression.
  For example, if \code{data = x}, you may specify any of
  \itemize{
    \item{\code{statistic = mean}}
    \item{\code{statistic = mean(x)}}
    \item{\code{statistic = mean(data)}}
  }

  If \code{data} is a data frame, the expression may refer to columns in
  the data frame, e.g.
  \itemize{
    \item{\code{statistic = mean(x)}}
    \item{\code{statistic = mean(myData$x)}}
    \item{\code{statistic = mean(myData[, "x"])}}
  }

  If \code{data} is not just the name of an object, e.g.
  \code{data = subset(myData, age > 17)}, or if \code{data2}
  is supplied, then use the name 'data', e.g.
  \itemize{
    \item{\code{statistic = colMeans(data)}}
  }
}
\value{
  a list with class
  \code{"bootstrap"}, \code{"bootstrap2"},
  \code{"permutationTest"}, \cr
  or \code{"permutationTest2"},
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

  The two-sample versions have an additional component:
  \item{resultsBoth}{containing resampling results from each data set.
  containing two components,
  the results from resampling each of the two samples. These are
  \code{bootstrap} objects; in the \code{permutationTest2} case
  they are the result of sampling without replacement.}

  There are functions for printing and plotting these objects,
  in particular \code{print}, \code{hist}, \code{qqnorm},
  \code{plot} (currently the same as \code{hist}),
  \code{quantile}.
}

\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
\code{\link{resample-package}},
\code{\link{samp.bootstrap}},
\code{\link{CI.percentile}},
\code{\link{CI.t}}.
}
\examples{
\donttest{
# See full set of examples in resample-package, including different
# ways to call the functions depending on the structure of the data.
data(Verizon)
CLEC <- with(Verizon, Time[Group == "CLEC"])
bootC <- bootstrap(CLEC, mean)
bootC
hist(bootC)
qqnorm(bootC)
}
}
\keyword{nonparametric}
\keyword{htest}
