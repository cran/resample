# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

resample <- function(data, resampleFun, sampler, B = 1000,
                     seed = NULL,
                     statisticNames = NULL,
                     block.size = 100, trace = FALSE, ...,
                     observedIndices = 1:n)
{
  # Perform resampling. This is the workhorse function, called by others.
  #
  # Args:
  #   data:            vector, matrix, or data frame
  #   resampleFun:     function(data, ii), created by .resampleMakeFun
  #   sampler:         function(n, B, ...) to generate samples
  #   B:               number of replications
  #   seed:            old value of .Random.seed, or argument to set.seed
  #   statisticNames:  names used for printing, character vector of length 'd'
  #   block.size:      replicates are done 'block.size' at a time
  #   trace:           logical, if TRUE an indication of progress is printed
  #   ...:             additional arguments passed to sampler
  #   observedIndices  vector of indices that correspond to the observed value.
  #                    This is used by permutationTest2

  .resampleSetSeed(seed)
  oldSeed <- .Random.seed

  dimData <- dim(data)
  n <- IfElse(is.null(dimData), length(data), dimData[1])

  # Observed value
  observed <- resampleFun(data, observedIndices)
  p <- length(observed)
  if(!is.null(dim(observed)))
    dim(observed) <- NULL
  if(!is.null(statisticNames) && length(statisticNames) != p)
    stop("statisticNames has the wrong length, is ",
         length(statisticNames),
         ", but observed has length ", p)
  if(is.null(statisticNames))
    statisticNames <- names(observed)
  if(is.null(statisticNames))
    statisticNames <- paste0("stat", 1:p)
  temp <- which(statisticNames == "")
  if(length(temp))
    statisticNames[temp] <- paste0("stat", temp)
  names(observed) <- statisticNames

  replicates <- matrix(vector(storage.mode(observed), B * p),
                       nrow = B, ncol = p,
                       dimnames = list(NULL, statisticNames))
  failures <- integer(0)

  # Do calculations, mostly block.size at a time.
  # Generate indices using sampler, and statistics using resampleFun.
  nBlocks <- ceiling(B / block.size)
  B1 <- 0
  for(iBlock in 1:nBlocks) {
    B0 <- B1 # Number of replications already performed
    B1 <- min(B, iBlock * block.size) # Last replication to do in this block
    if(trace)
      cat("Calculating replications ", B0 + 1, ":", B1, "\n", sep = "")
    BB <- B1 - B0
    indices <- sampler(n, BB, ...)
    for(j in 1:BB) {
      theta <- try(resampleFun(data, indices[,j]), silent = TRUE)
      if(!is(theta, "try-error")) {
        replicates[B0 + j, ] <- theta
      } else {
        failures <- c(failures, B0 + j)
      }
    }
  }
  result <- list(observed = observed,
                 replicates = replicates,
                 n = n, p = p, B = B,
                 seed = oldSeed)
  if(!missing(observedIndices)) {
    # Working with a subset of a larger dataset.
    result$n <- length(observedIndices)
    result$nCombined <- n
  }
  class(result) <- "resample"
  return(result)
}


print.resample <- function(x, ...) {
  # This is mostly used for objects with $call and $stats, like "bootstrap"
  if(!is.null(x$call))
    cat0n("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catn("Replications:", x$B)
  if(is.null(x$ratio)) { # one-sample application
    catn("\nSummary Statistics:")
  } else {
    catn("Two samples, sample sizes are", x$n[1], x$n[2])
    catn("\nSummary Statistics for the",
         IfElse(x$ratio, "ratio", "difference"),
         "between samples 1 and 2:")
  }
  print(IfElse(is.null(x$stats),
               .resampleBootStats(x)[1:3], # Observed SE Mean, not Bias
               x$stats), ...)
  invisible(x)
}

plot.resample <- function(x, ..., resampleColumns = 1:x$p, xlim = NULL,
                          xlab = NULL, main = "") {
  statNames <- names(x$observed)
  nullXlim <- is.null(xlim)
  for(j in resampleColumns) {
    xj <- x$replicates[, j]
    if(nullXlim)
      xlim <- range(xj, x$observed[j])
    hist(xj, probability = TRUE,
         xlab = IfElse(is.null(xlab), statNames[j], xlab),
         xlim = xlim, ..., main = main)
    lines(density(xj, from = xlim[1], to = xlim[2]))
    abline(v = x$observed[j], lwd = 2)
  }
  invisible(NULL)
}

hist.resample <- plot.resample

qqnorm.resample <- function(y, ..., resampleColumns = 1:y$p,
                            ylab = NULL) {
  statNames <- names(y$observed)
  for(j in resampleColumns) {
    yj <- y$replicates[, j]
    qqnorm(yj, ...,
           ylab = IfElse(is.null(ylab), statNames[j], ylab))
  }
  invisible(NULL)
}

quantile.resample <- function(x, ..., type = 6) {
  apply(x$replicates, 2, quantile, ..., type = type)
}


if(FALSE) { # manual testing code
  r <- bootstrap((xDF), colMeans) # xDF defined in bootstrap.R
  plot(r, resampleColumns = 2)
  qqnorm(r, resampleColumns = 2)
  quantile(r, probs = c(.025, .975))
  limits.percentile(r, probs = c(.025, .975))

  r <- bootstrap2(treatment = t9, x9, mean)
  plot(r)

  r <- permutationTest2(treatment = t9, x9, mean)
  plot(r)
}
