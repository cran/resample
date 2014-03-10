# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

limits.percentile <- function(x, probs = c(.025, .975), ...) {
  quantile(x, probs = probs, ...)
}

limits.t <- function(x, probs = c(.025, .975)) {
  # T interval with bootstrap SE
  obs <- x$observed
  SE <- x$stats$SE
  df <- min(x$n) - 1
  result <- rep(obs, each = length(probs)) + qt(probs, df = df) %o% SE
  dimnames(result) <- list(paste0(100*probs, "%"), names(obs))
  result
}


if(FALSE) {
  r <- bootstrap2(treatment = t9, xDF, colMeans)
  limits.t(r)
  limits.percentile(r)
}
