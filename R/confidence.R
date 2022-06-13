# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

CI.percentile <- function(x, confidence = 0.95, expand = TRUE, ...,
                          probs = sort(1 + c(-1, 1) * confidence) / 2) {
  # Bootstrap percentile confidence interval
  # Args:
  #   x: a bootstrap or bootstrap2 object
  #   confidence: confidence level, e.g. 0.95
  #   expand: logical, if TRUE then use modified percentiles for better
  #              small-sample accuracy.
  #   ...: passed to Quantile (e.g. na.rm=TRUE - but this behavior may change)
  #   probs: confidence level probabilities, e.g. c(0.025, 0.975)
  # Return:
  #   matrix, rows correspond to statistics and columns to probs
  probs2 <- IfElse(expand, ExpandProbs(probs, min(x$n)), probs)
  result <- Quantile(x, probs = probs2, ...)
  dimnames(result) <- list(names(x$observed), .FormatProbs(probs))
  return(result)
}
# TODO: separate calibration for t and SE effects

ExpandProbs <- function(probs, n) {
  # Modify the probs to adjust for two factors:
  # Bootstrap distributions are too narrow by a factor of sqrt((n-1)/n)
  # Percentile interval corresponds to estimate +- z_alpha instead of t_alpha s.
  # Find probs2 such that qnorm(probs2) * sqrt((n-1)/n) = qt(probs, n-1)
  pnorm(qt(probs, n-1) * sqrt(n / (n-1)))
}

.FormatProbs <- function(probs) {
  # Format quantiles, e.g. 0.025 becomes 2.5%
  paste0(formatC(100 * probs, format = "fg", width = 1,
                 max(2L, getOption("digits"))), "%")
}

CI.t <- function(x, confidence = 0.95, expand = TRUE,
                 probs = sort(1 + c(-1, 1) * confidence) / 2) {
  # T interval with bootstrap SE
  # Args:
  #   x: a bootstrap or bootstrap2 object
  #   confidence: confidence level, e.g. 0.95
  #   probs: confidence level probabilities, e.g. c(0.025, 0.975)
  # Return:
  #   matrix, rows correspond to statistics and columns to probs
  obs <- x$observed
  n <- min(x$n) - 1
  SE <- x$stats$SE * IfElse(expand, sqrt(n/(n-1)), 1)
  result <- rep(obs, each = length(probs)) + qt(probs, df = n-1) %o% SE
  dimnames(result) <- list(.FormatProbs(probs), names(obs))
  t(result)
}


CI.bootstrapT <- function(x, confidence = 0.95,
                          probs = sort(1 + c(-1, 1) * confidence) / 2) {
  # Bootstrap t interval
  # This assumes that the first dimension of the statistic
  # is an estimate, and the second is proportional to a SE for the estimate.
  # E.g. for bootstrapping the mean, they could be the mean and s.
  # Args:
  #   x: bootstrap object
  tAlpha <- Quantile((x$replicates[, 1] - x$observed[1]) / x$replicates[, 2],
                     1-probs)
  result <- x$observed[1] - tAlpha * x$observed[2]
  matrix(result, 1, dimnames = list(names(x$observed[1]), .FormatProbs(probs)))
}
# S+Resample has bootstrapT, with a pivot argument.


CI.bca <- function(x, confidence = 0.95,
                   expand = TRUE, L = NULL,
                   probs = sort(1 + c(-1, 1) * confidence) / 2) {
  # Bootstrap BCa confidence interval
  # Args:
  #   x: a bootstrap object
  #   confidence: confidence level, e.g. 0.95
  #   probs: confidence levels
  #   L:     empirical influence function
  # Return:
  #   matrix, rows correspond to statistics and columns to probs
  if(inherits(x, "bootstrap2"))
    stop("bootstrap2 objects are not currently supported")

  skewness <- function(x) {
    x <- x - mean(x)
    mean(x^3) / mean(x^2)^1.5
  }
  if(is.null(L)){
    Call <- x$call
    Call[[1]] <- as.name("jackknife")
    Call$R <- NULL
    Call$seed <- NULL
    Call$sampler <- NULL
    Call$block.size <- NULL
    jackObject <- eval(Call, sys.parent())
    L <- -jackObject$replicates
  }
  a <- apply(L, 2, skewness) / (6 * sqrt(nrow(L)))
  w <- qnorm(colMeans(x$replicates < rep(x$observed, each = x$R)))
  probs2 <- IfElse(expand, ExpandProbs(probs, min(x$n)), probs)
  zalpha <- qnorm(probs2)
  # Statistics (a, w) in rows, probs (and zalpha) in columns.
  result <- matrix(NA, nrow = x$p, ncol = length(probs),
                   dimnames = list(names(x$observed), .FormatProbs(probs)))
  for(j in 1:x$p) {
    probs3 <- pnorm(w[j] + (w[j] + zalpha) / (1 - a[j] * (w[j] + zalpha)))
    result[j, ] <- Quantile(x$replicates[, j], probs3)
  }
  return(result)
}
# TODO: more testing.
# TODO: Support two-sample applications.



if(FALSE) {
  x9 <- sqrt(1:999)  # skewed
  xDF <- data.frame(X = x9, Y = max(x9) - x9)

  r <- bootstrap(x9, mean, R=100, seed = 0)
  CI.t(r)
  CI.percentile(r)
  CI.bca(r)

  CI.t(r, confidence = c(.90, .95))
  CI.percentile(r, confidence = c(.90, .95))
  CI.bca(r, confidence = c(.90, .95))

  CI.t(r, probs = c(.8, .9, .95))
  CI.percentile(r, probs = c(.8, .9, .95))
  CI.bca(r, probs = c(.8, .9, .95))


  r2 <- bootstrap(xDF, colMeans, R=100, seed = 0)
  CI.t(r2)
  CI.percentile(r2)
  CI.bca(r2)

  CI.t(r2, probs = c(.8, .9, .95))
  CI.percentile(r2, probs = c(.8, .9, .95))
  CI.bca(r2, probs = c(.8, .9, .95))

  r3 <- bootstrap(xDF, mean(Y), R=100, seed = 0)
  CI.t(r3, probs = c(.8, .9, .95))
  CI.percentile(r3, probs = c(.8, .9, .95))
  CI.bca(r3, probs = c(.8, .9, .95))

  all.equal(unname(CI.t(r2)),
            unname(rbind(CI.t(r), CI.t(r3))))
  all.equal(unname(CI.percentile(r2)),
            unname(rbind(CI.percentile(r), CI.percentile(r3))))
  all.equal(unname(CI.bca(r2)),
            unname(rbind(CI.bca(r), CI.bca(r3))))

  # Check multiple statistics
  r <- bootstrap2(treatment = t9, xDF, colMeans, R=100)
  CI.t(r)
  CI.percentile(r)
  CI.bca(r)

  r <- bootstrap((1:10)^2, c(mean(data), sd(data)), R=100)
  CI.bootstrapT(r)
  CI.bootstrapT(r, probs = .9)
  CI.bca(r)

}


if(FALSE){ # Development work, using examples in bootstrap/howToTeach
  for(i in list.files(path = "~/Rlang/resample/resample/resample/R",
                      pattern = "*.R$",
                      full.names = TRUE)) source(i)

  jackBasic <- jackknife(Basic, mean)
  jack2 <- jackknife(Basic, function(x) mean((x-mean(x))^2))
  jack2
  # Compare that bias to
  sd(Basic)^2 - mean((Basic-mean(Basic))^2)
  CI.t(jackBasic)
  CI.t(jack2)

  boot <- bootstrap(Basic, c(mean = mean(Basic), median = median(Basic)))
  jack <- jackknife(Basic, c(mean = mean(Basic), median = median(Basic)))
  CI.percentile(bootBasic)
  CI.percentile(boot)
  CI.percentile(bootBasic, expand = TRUE)

  CI.bca(bootBasic)
  CI.bca(boot)
}
