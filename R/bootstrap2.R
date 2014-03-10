# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

bootstrap2 <- function(data, statistic,
                       treatment, data2 = NULL,
                       B = 1000, ratio = FALSE, args.stat = NULL,
                       seed = NULL,
                       sampler = samp.bootstrap,
                       label = NULL, statisticNames = NULL,
                       block.size = 100, trace = FALSE)
{
  # Two-sample nonparametric bootstrapping
  # Either specify (data, data2) or (data, treatment).
  #
  # Args:
  #   data:      vector, matrix, or data frame.
  #              Let 'nObs' = length of vector, or nrow otherwise.
  #   statistic: a function, or expression (e.g. mean(data, trim = .2)
  #              If data is a data frame, can refer to variables in it.
  #              This may be a vector; let 'd' be its length.
  #   treatment: a vector of length nObs, with two unique values.
  #   data2:     like data.
  #   B:         number of replications.
  #   ratio:     logical, if FALSE then bootstrap the difference between
  #              data and data2 (or first and second treatment).
  #              If TRUE then bootstrap the ratio.
  #   args.stat: additional arguments to pass to the function.
  #   seed:      old value of .Random.seed, or argument to set.seed.
  #   sampler:   a function for resampling, see help(samp.bootstrap).
  #   label:     used for labeling plots.
  #   statisticNames: names used for printing, character vector of length 'd'.
  #   block.size: replicates are done 'block.size' at a time.
  #   trace:     logical, if TRUE an indication of progress is printed.

  Call <- match.call()

  dimData <- dim(data)
  n <- IfElse(is.null(dimData), length(data), dimData[1])
  stopifnot(length(dimData) <= 2)

  resultsBoth <- vector("list", 2)
  if(is.null(data2)) {                  # use treatment
    if(is.data.frame(data))
      treatment <- eval(substitute(treatment),
                        envir = data, enclos = parent.frame())
    treatmentInds <- split(1:n, treatment)
    treatmentNames <- names(treatmentInds)
    if(length(treatmentInds) != 2)
      stop("treatment must have 2 unique values, observed: ",
           paste(treatmentNames, collapse = " "))
    for(k in 1:2) {
      # Create a smaller data set corresponding to this treatment
      dataK <- IfElse(length(dimData) == 2,
                      data[treatmentInds[[k]], , drop = FALSE],
                      data[treatmentInds[[k]]])
      resampleFun <-
        .resampleMakeFunction(dataK, statistic,
                              substitute(data), substitute(statistic),
                              args.stat)
      resultsBoth[[k]] <-
        resample(dataK, resampleFun, sampler = sampler, B = B,
                 seed = IfElse(k == 1, seed, .Random.seed),
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace)
      resultsBoth[[k]]$call <-
        paste("bootstrap for treatment", treatmentNames[k])
    }
  } else {  # data and data2
    if(length(dim(data2)) != length(dimData) ||
       class(data2) != class(data))
      stop("data and data2 must have similar structure")
    treatmentNames <- c(IfElse(is.name(substitute(data)),
                               substitute(data), "data"),
                        IfElse(is.name(substitute(data2)),
                               substitute(data2), "data2"))
    for(k in 1:2) {
      resampleFun <-
        .resampleMakeFunction(IfElse(k == 1, data, data2), statistic,
                              as.name("data"), substitute(statistic), args.stat)
      resultsBoth[[k]] <-
        resample(IfElse(k == 1, data, data2),
                 resampleFun, sampler = sampler, B = B,
                 seed = IfElse(k == 1, seed, .Random.seed),
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace)
      resultsBoth[[k]]$call <-
        paste("bootstrap for data set", k)
    }
  }
  for(k in 1:2) {
    resultsBoth[[k]]$stats <- .resampleBootStats(resultsBoth[[k]])
    class(resultsBoth[[k]]) <- c("bootstrap", "resample")
  }
  names(resultsBoth) <- treatmentNames

  if(!identical(names(resultsBoth[[1]]$observed),
                names(resultsBoth[[2]]$observed))) {
    warning("Statistics returned from the two samples are not compatible,",
            "will return results so far, without combining them.")
    print(all.equal(resultsBoth[[1]]$observed,
                    resultsBoth[[2]]$observed))
    return(resultsBoth)
  }

  # combine results
  op <- get(IfElse(ratio, "/", "-"))
  result <- list(observed =
                 op(resultsBoth[[1]]$observed, resultsBoth[[2]]$observed),
                 replicates =
                 op(resultsBoth[[1]]$replicates, resultsBoth[[2]]$replicates),
                 n = c(resultsBoth[[1]]$n, resultsBoth[[2]]$n),
                 p = resultsBoth[[1]]$p,
                 B = B,
                 seed = resultsBoth[[1]]$seed,
                 call = match.call(),
                 resultsBoth = resultsBoth,
                 ratio = ratio)
  result$call <- match.call()
  result$stats <- .resampleBootStats(result)
  class(result) <- c("bootstrap2", "bootstrap", "resample")
  result
}

# Note difference:
#   bootstrap(x, mean(x))
#   bootstrap2(x1, data2 = x2, mean(data))


# print.bootstrap2 <- function(x, ...) {
#   cat0n("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
#   catn("Replications:", x$B)
#   catn("Two groups, sample sizes are", x$n[1], x$n[2])
#   catn("\nSummary Statistics for the", IfElse(x$ratio, "ratio", "difference"),
#        "between samples 1 and 2:")
#   print(x$stats, ...)
#   invisible(x)
# }
# same as print.permutatonTest2


if(FALSE) {
  x9 <- 1:9
  xDF <- data.frame(a=x9, b=2*x9)
  t9 <- letters[c(1,1,1,2,1,2,1,2,1)]
  xDF1 <- xDF[t9 == 1, ]
  xDF2 <- xDF[t9 == 2, ]
  x91 <- x9[t9 == 1]
  x92 <- x9[t9 == 2]

  source("~/resample/R/bootstrap2.R")

  ##### treatment
  ### statistic by name
  # base case: data by name, statistic is function by name
  bootstrap2(treatment = t9, x9, mean)

  temp <- .Last.value
  sapply(temp$resultsBoth, class)
  temp$resultsBoth

  # data expression
  bootstrap2(treatment = t9, (x9), mean)

  # args.stat
  bootstrap2(treatment = t9, x9, mean, args.stat = list(trim = .25))

  # inline function
  bootstrap2(treatment = t9, x9, function(z) mean(z))

  # data frame
  bootstrap2(treatment = t9, xDF, colMeans)

  # data expression, data frame
  bootstrap2(treatment = t9, (xDF), colMeans)

  # data expression, matrix
  bootstrap2(treatment = t9, as.matrix(xDF), colMeans)


  ### statistic expression
  # data by name
  bootstrap2(treatment = t9, x9, mean(x9))

  # data as expression, refer to 'data'
  bootstrap2(treatment = t9, (x9), mean(data))

  # data frame
  bootstrap2(treatment = t9, xDF, mean(a))

  # data frame expression
  bootstrap2(treatment = t9, (xDF), mean(a))

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))


  ##### data and data2
  ### statistic by name
  # base case: data by name, statistic is function by name
  bootstrap2(x91, data2 = x92, mean)

  temp <- .Last.value
  sapply(temp$resultsBoth, class)
  temp$resultsBoth

  # data expression
  bootstrap2((x91), data2 = (x92), mean)

  # args.stat
  bootstrap2(x91, data2 = x92, mean, args.stat = list(trim = .25))

  # inline function
  bootstrap2(x91, data2 = x92, function(z) mean(z))

  # data frame
  bootstrap2(xDF1, data2 = xDF2, colMeans)

  # data expression, data frame
  bootstrap2((xDF1), data2 = (xDF2), colMeans)

  # data expression, matrix
  bootstrap2(as.matrix(xDF1), data2 = as.matrix(xDF2), colMeans)


  ### statistic expression
  # data by name
  bootstrap2(x91, data2 = x92, mean(data))  # use 'data' in this case

  # data as expression, refer to 'data'
  bootstrap2((x91), data2 = (x92), mean(data))

  # data frame
  bootstrap2(xDF1, data2 = xDF2, mean(a))

  # data frame expression
  bootstrap2((xDF1), data2 = (xDF2), mean(a))


  source("~/resample/R/bootstrap2.R")
}