% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{deprecated.resample}
\alias{limits.percentile}
\alias{limits.t}
\alias{limits.bootstrapT}
\title{
Deprecated functions.
}
\description{
Deprecated functions
}
% \usage{
% }
\arguments{
  \item{\dots}{
    arguments to pass to the replacement functions.
}
}
\value{
  See the replacement functions.
}
\details{
  \code{limits.percentile}, \code{limits.t} and \code{limits.bootstrapT}
  have been renamed "CI.*".
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
\code{\link{CI.percentile}},
\code{\link{CI.t}},
\code{\link{CI.bootstrapT}}.
}
% \examples{
% \donttest{
% # See full set of examples in resample-package, including different
% # ways to call all four functions depending on the structure of the data.
% data(Verizon)
% CLEC <- with(Verizon, Time[Group == "CLEC"])
% bootC <- bootstrap(CLEC, mean, seed = 0)
% bootC2 <- bootstrap(CLEC, c(mean = mean(CLEC), sd = sd(CLEC)), seed = 0)
% CI.percentile(bootC)
% CI.t(bootC)
% CI.bootstrapT(bootC2)
% }
% }
\keyword{nonparametric}
\keyword{htest}
