% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{resample-data}
\alias{resample-data}
\alias{TV}
\alias{Verizon}
\title{
Data sets for resampling examples
}
\description{
Data sets for use in examples.
}
\section{DATA SETS}{
TV            10 observations: Time,Cable
Verizon     1687 observations: Time,Group
}
\source{
  The TV and Verizon datasets are used in What Teachers Should Know
  about the Bootstrap: Resampling in the Undergraduate Statistics
  Curriculum
}
\details{
\code{TV} has measurements of minutes of commercials per half-hour, for
"Basic" and "Extended" (extra-cost) cable TV stations.

\code{Verizon} has repair times, with two groups, CLEC and ILEC,
customers of the "Competitive" and "Incumbent" local exchange carrior.
  }
\references{
  Hesterberg, Tim (2014), What Teachers Should Know about the
  Bootstrap: Resampling in the Undergraduate Statistics Curriculum,
  \url{http://arxiv.org/abs/1411.5279}.
}
\seealso{
See \code{\link{resample-package}} for an overview of resampling functions.
}
\examples{
\dontrun{
data(TV); summary(TV)
Basic <- with(TV, Time[Cable == "Basic"])
Extended <- with(TV, Time[Cable == "Extended"])

data(Verizon); summary(Verizon)
ILEC <- with(Verizon, Time[Group == "ILEC"])
CLEC <- with(Verizon, Time[Group == "CLEC"])
}
}
\keyword{datasets}
