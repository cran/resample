% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{resample-package}
\alias{resample-package}
\docType{package}
\title{Resampling Functions}
\description{
  Resampling functions, including one- and two-sample bootstrap and
  permutation tests, with an easy-to-use syntax.
}
\details{
  See \code{library(help = resample)} for version number, date, etc.
}
\section{Main resampling functions}{
  The main resampling functions are:
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{permutationTest}},
  \code{\link{permutationTest2}}.
}
\section{Methods}{
  Methods for generic functions include:
  \code{\link{print.resample}},
  \code{\link{plot.resample}},
  \code{\link{hist.resample}},
  \code{\link{qqnorm.resample}},
  \code{\link{quantile.resample}}
}
\section{Confidence Intervals}{
  Functions that calculate confidence intervals for \code{\link{bootstrap}}
  and \code{\link{bootstrap2}} objects:
  \code{\link{limits.percentile}},
  \code{\link{limits.t}}.
}
\section{Samplers}{
  Functions that generate indices for random samples:
  \code{\link{samp.bootstrap}},
  \code{\link{samp.permute}}.
}
\section{Low-level Resampling Function}{
  This is called by the main resampling functions, but can also be
  called directly:
  \code{\link{resample}}.
}
\author{
  Tim Hesterberg <timhesterberg@gmail.com>
}
\keyword{nonparametric}
\keyword{htest}
