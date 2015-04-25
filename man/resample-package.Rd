% Copyright 2014 Google Inc. All rights reserved.
%
% Use of this source code is governed by a BSD-style
% license that can be found in the LICENSE file or at
% http://opensource.org/licenses/BSD-3-Clause

\name{resample-package}
\alias{resample-package}
\docType{package}
\title{Overview of the resample package}
\description{
  Resampling functions, including one- and two-sample bootstrap and
  permutation tests, with an easy-to-use syntax.
}
\details{
  See \code{library(help = resample)} for version number, date, etc.
}
\section{Data Sets}{
  A list of datasets is at
  \code{\link{resample-data}},
  }
\section{Main resampling functions}{
  The main resampling functions are:
  \code{\link{bootstrap}},
  \code{\link{bootstrap2}},
  \code{\link{permutationTest}}, \cr
  \code{\link{permutationTest2}}.
}
\section{Methods}{
  Methods for generic functions include:
  \code{\link{print.resample}},
  \code{\link{plot.resample}}, \cr
  \code{\link{hist.resample}},
  \code{\link{qqnorm.resample}}, and
  \code{\link{quantile.resample}}.
}
\section{Confidence Intervals}{
  Functions that calculate confidence intervals for \code{\link{bootstrap}}
  and \code{\link{bootstrap2}} objects: \cr
  \code{\link{CI.bca}},
  \code{\link{CI.bootstrapT}},
  \code{\link{CI.percentile}},
  \code{\link{CI.t}}.
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
\section{New Versions}{
  I will post the newest versions to
  \url{http://www.timhesterberg.net/r-packages}.
  See that page to join a list for announcements of new versions.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\examples{
data(Verizon)
ILEC <- with(Verizon, Time[Group == "ILEC"])
CLEC <- with(Verizon, Time[Group == "CLEC"])

#### Sections in this set of examples
### Different ways to specify the data and statistic
### Example with plots and confidence intervals.


### Different ways to specify the data and statistic
# This code is flexible; there are different ways to call it,
# depending on how the data are stored and on the statistic.

## One-sample Bootstrap

\dontrun{
# Ordinary vector, give statistic as a function
bootstrap(CLEC, mean)

# Vector by name, give statistic as an expression
bootstrap(CLEC, mean(CLEC))

# Vector created by an expression, use the name 'data'
bootstrap(with(Verizon, Time[Group == "CLEC"]), mean(data))

# A column in a data frame; use the name of the column
temp <- data.frame(foo = CLEC)
bootstrap(temp, mean(foo))

# Put function arguments into an expression
bootstrap(CLEC, mean(CLEC, trim = .25))

# Put function arguments into a separate list
bootstrap(CLEC, mean, args.stat = list(trim = .25))
}

\dontshow{
# Ordinary vector, give statistic as a function
bootstrap(CLEC, mean, R = 100)

# Vector by name, give statistic as an expression
bootstrap(CLEC, mean(CLEC), R = 100)

# Vector created by an expression, use the name 'data'
bootstrap(with(Verizon, Time[Group == "CLEC"]), mean(data), R = 100)

# A column in a data frame; use the name of the column
temp <- data.frame(foo = CLEC)
bootstrap(temp, mean(foo), R = 100)

# Put function arguments into an expression
bootstrap(CLEC, mean(CLEC, trim = .25), R = 100)

# Put function arguments into a separate list
bootstrap(CLEC, mean, args.stat = list(trim = .25), R = 100)
}


## One-sample jackknife

# Like bootstrap. E.g.
jackknife(CLEC, mean)


## One-sample permutation test

# To test H0: two variables are independent, exactly
# one of them just be permuted. For the CLEC data,
# we'll create an artificial variable.
CLEC2 <- data.frame(Time = CLEC, index = 1:length(CLEC))

\dontrun{
permutationTest(CLEC2, cor(Time, index),
                resampleColumns = "index")
# Could permute "Time" instead.

# resampleColumns not needed for variables outside 'data'
permutationTest(CLEC, cor(CLEC, 1:length(CLEC)))
}

\dontshow{
permutationTest(CLEC2, cor(Time, index),
                resampleColumns = "index", R = 99)
# Could permute "Time" instead.

# resampleColumns not needed for variables outside 'data'
permutationTest(CLEC, cor(CLEC, 1:length(CLEC)), R = 99)
}


### Two-sample problems
## Different ways to specify data and statistic

## Two-sample bootstrap

# Two data objects (one for each group)
\dontshow{bootstrap2(CLEC, data2 = ILEC, mean, R=100)}
\dontrun{bootstrap2(CLEC, data2 = ILEC, mean)}

# data frame containing y variable(s) and a treatment variable
\dontshow{bootstrap2(Verizon, mean(Time), treatment = Group, R=100)}
\dontrun{bootstrap2(Verizon, mean(Time), treatment = Group)}

# treatment variable as a separate object
temp <- Verizon$Group
\dontshow{bootstrap2(Verizon$Time, mean, treatment = temp, R=100)}
\dontrun{bootstrap2(Verizon$Time, mean, treatment = temp)}


## Two-sample permutation test

# Like bootstrap2. E.g.
\dontshow{permutationTest2(CLEC, data2 = ILEC, mean, R=99)}
\dontrun{permutationTest2(CLEC, data2 = ILEC, mean}


### Example with plots and confidence intervals.
\dontrun{
boot <- bootstrap2(CLEC, data2 = ILEC, mean)
perm <- permutationTest2(CLEC, data2 = ILEC, mean,
                         alternative = "greater")
}
\dontshow{
boot <- bootstrap2(CLEC, data2 = ILEC, mean, R=100)
perm <- permutationTest2(CLEC, data2 = ILEC, mean,
                         alternative = "greater", R=99)
}
\dontrun{
par(mfrow = c(2,2))
hist(boot)
qqnorm(boot)
qqline(boot$replicates)
hist(perm)
}

# P-value
perm

# Standard error, and bias estimate
boot

# Confidence intervals
CI.percentile(boot) # Percentile interval
CI.t(boot)  # t interval using bootstrap SE
# CI.bootstrapT and CI.bca do't currently support two-sample problems.

# Statistic can be multivariate.
# For the bootstrap2, it must have the estimate first, and a standard
# error second (don't need to divide by sqrt(n), that cancels out).
bootC <- bootstrap(CLEC, mean, seed = 0)
bootC2 <- bootstrap(CLEC, c(mean = mean(CLEC), sd = sd(CLEC)), seed = 0)
identical(bootC$replicates[, 1], bootC2$replicates[, 1])

CI.percentile(bootC)
CI.t(bootC)
CI.bca(bootC)
CI.bootstrapT(bootC2)
# The bootstrapT is the most accurate for skewed data, especially
# for small samples.

# By default the percentile interval is "expanded", for better coverage
# in small samples. To turn this off:
CI.percentile(bootC, expand = FALSE)
}
\keyword{nonparametric}
\keyword{htest}
