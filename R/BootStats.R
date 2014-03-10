# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

.resampleBootStats <- function(x) {
  # return a data frame with basic summary statistics
  meanReplicates <- colMeans(x$replicates)
  data.frame(row.names = names(x$observed),
             Observed = x$observed,
             SE = colStdevs(x$replicates),
             Mean = meanReplicates,
             Bias = meanReplicates - x$observed)
}

.resamplePermutationStats <- function(x, alternative) {
  # return a data frame with basic summary statistics
  #
  # Args:
  #   x: a permutationTest or permutationTest2 object
  #   alternative: character vector: "greater", "less", or "two.sided"
  #      This is replicated to length = number of statistics

  if(any(!is.element(alternative, c("greater", "less", "two.sided"))))
    warning("alternative is not a recognized value")
  alternative <- rep(alternative, length = x$p)

  B <- x$B
  o <- x$observed
  r <- x$replicates
  pValueLess <- (1 + colSums(r <= rep(o, each=B))) / (B+1)
  pValuePlus <- (1 + colSums(r >= rep(o, each=B))) / (B+1)
  pValue <- rep(NA, x$p)
  pValue[alternative == "less"] <- pValueLess[alternative == "less"]
  pValue[alternative == "greater"] <- pValuePlus[alternative == "greater"]
  pValue[alternative == "two.sided"] <- pmin(1,
           2 * pmin(pValueLess, pValuePlus))[alternative == "two.sided"]

  data.frame(row.names = names(o),
             Observed = o,
             Mean = colMeans(x$replicates),
             Alternative = alternative,
             PValue = pValue)
}
