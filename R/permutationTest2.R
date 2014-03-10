# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

permutationTest2 <- function(data, statistic,
                             treatment, data2 = NULL,
                             B = 999,
                             alternative = "two.sided",
                             ratio = FALSE, paired = FALSE,
                             args.stat = NULL,
                             seed = NULL,
                             sampler = samp.permute,
                             label = NULL, statisticNames = NULL,
                             block.size = 100, trace = FALSE)
{
  # Two-sample permutation tests
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
  #   B:         number of replications
  #   alternative: one of "two.sided", "greater", or "less"
  #   ratio:     logical, if FALSE then resample the difference between
  #              data and data2 (or first and second treatment).
  #              If TRUE then resample the ratio.
  #   paired:    logical, if TRUE then data are paired, permute within pairs.
  #   args.stat: additional arguments to pass to the function.
  #   seed:      old value of .Random.seed, or argument to set.seed.
  #   sampler:   a function for resampling, see help(samp.permute).
  #   label:     used for labeling plots.
  #   statisticNames: names used for printing, character vector of length 'd'.
  #   block.size: replicates are done 'block.size' at a time.
  #   trace:     logical, if TRUE an indication of progress is printed.

  Call <- match.call()

  dimData <- dim(data)
  stopifnot(length(dimData) <= 2)
  if(paired) stop("paired is not yet implemented") # TODO: implement it

  resultsBoth <- vector("list", 2)
  if(is.null(data2)) {                  # use treatment
    resampleFun <-
      .resampleMakeFunction(data, statistic,
                            substitute(data), substitute(statistic), args.stat)
    if(is.data.frame(data))
      treatment <- eval(substitute(treatment),
                        envir = data, enclos = parent.frame())
    n <- length(treatment)
    treatmentInds <- split(1:n, treatment)
    treatmentNames <- names(treatmentInds)
    if(length(treatmentNames) != 2)
      stop("treatment must have 2 unique values, observed: ",
           paste(treatmentNames, collapse = " "))
    nBoth <- sapply(treatmentInds, length)
    stopifnot(all(nBoth >= 0))
    for(k in 1:2) {
      samplerK <- eval(bquote(function(...)
                              sampler(..., groupSizes = .(nBoth),
                                      returnGroup = .(k))))
      resultsBoth[[k]] <-
        resample(data, resampleFun, sampler = samplerK, B = B,
                 observedIndices = treatmentInds[[k]],
                 seed = seed,
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace)
      resultsBoth[[k]]$call <-
        paste("resample for treatment", treatmentNames[k])
    }
  } else {  # data and data2
    if(length(dim(data2)) != length(dimData) ||
       class(data2) != class(data))
      stop("data and data2 must have similar structure")
    treatmentNames <- c(IfElse(is.name(substitute(data)),
                               as.character(substitute(data)), "data"),
                        IfElse(is.name(substitute(data2)),
                               as.character(substitute(data2)), "data2"))
    dataBoth <- IfElse(length(dimData) == 2, rbind(data, data2),
                       c(data, data2))
    nBoth <- IfElse(length(dimData) == 2, c(nrow(data), nrow(data2)),
                    c(length(data), length(data2)))
    stopifnot(all(nBoth >= 0))
    resampleFun <-
      .resampleMakeFunction(dataBoth, statistic,
                            as.name("data"), substitute(statistic), args.stat)
    for(k in 1:2) {
      samplerK <- eval(bquote(function(...)
                              sampler(..., groupSizes = .(nBoth),
                                      returnGroup = .(k))))
      resultsBoth[[k]] <-
        resample(dataBoth, resampleFun, sampler = samplerK, B = B,
                 observedIndices = IfElse(k == 1, 1:nBoth[1],
                   nBoth[1] + 1:nBoth[2]),
                 seed = seed,
                 statisticNames = statisticNames,
                 block.size = block.size, trace = trace)
      resultsBoth[[k]]$call <-
        paste("resample for data set", k)
    }
  }
  for(k in 1:2) {
    # Compute individual P-values. Use alternative="two.sided" because there
    # are too many possibilities for what people would want for the individual
    # results, depending on k, ratio, and which replicates are negative.
    resultsBoth[[k]]$stats <-
      .resamplePermutationStats(resultsBoth[[k]], alternative = "two.sided")
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
                 n = sum(nBoth),
                 p = resultsBoth[[1]]$p,
                 B = B,
                 seed = resultsBoth[[1]]$seed,
                 call = match.call(),
                 resultsBoth = resultsBoth,
                 ratio = ratio)
  result$call <- match.call()
  result$stats <- .resamplePermutationStats(result, alternative = alternative)
  class(result) <- c("permutationTest2", "permutationTest", "resample")
  result
}


# print.permutationTest2 <- function(x, ...) {
#   cat0n("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
#   catn("Replications:", x$B)
#   catn("Two groups, sample sizes are", x$n[1], x$n[2])
#   catn("\nSummary Statistics for the", IfElse(x$ratio, "ratio", "difference"),
#        "between samples 1 and 2:")
#   print(x$stats, ...)
#   invisible(x)
# }
# same as print.bootstrap2


if(FALSE) {
  x9 <- 1:9
  xDF <- data.frame(a=x9, b=2*x9)
  t9 <- letters[c(1,1,1,2,1,2,1,2,1)]

  x1 <- x9[t9 == "a"]
  x2 <- x9[t9 == "b"]
  xDF1 <- data.frame(a = x1)
  xDF2 <- data.frame(a = x2)


  source("~/resample/R/permutationTest2.R")

  ##### treatment
  ### statistic by name
  # base case: data by name, statistic is function by name
  permutationTest2(treatment = t9, x9, mean)

  temp <- .Last.value
  sapply(temp$resultsBoth, class)
  temp$resultsBoth

  # data expression
  permutationTest2(treatment = t9, (x9), mean)

  # args.stat
  permutationTest2(treatment = t9, x9, mean, args.stat = list(trim = .25))

  # inline function
  permutationTest2(treatment = t9, x9, function(z) mean(z))

  # data frame,
  permutationTest2(treatment = t9, xDF, colMeans)

  # data expression, data frame
  permutationTest2(treatment = t9, xDF, colMeans)

  # data expression, matrix
  permutationTest2(treatment = t9, as.matrix(xDF), colMeans)

  # data frame, treatment in data frame
  permutationTest2(treatment = tt9, cbind(xDF, tt9 = t9),
                   function(x) colMeans(x[1:2]))


  ### statistic expression
  # data by name
  permutationTest2(treatment = t9, x9, mean(x9))

  # data as expression, refer to 'data'
  permutationTest2(treatment = t9, (x9), mean(data))

  # data frame
  permutationTest2(treatment = t9, xDF, mean(a))

  # data frame expression
  permutationTest2(treatment = t9, (xDF), mean(a))

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))


  ##### data and data2
  ### statistic by name
  # base case: data by name, statistic is function by name
  permutationTest2(x1, data2 = x2, mean)

  # data expression
  permutationTest2((x1), data2 = (x2), mean)

  # args.stat
  permutationTest2(x1, data2 = x2, mean, args.stat = list(trim = .25))

  # inline function
  permutationTest2(x1, data2 = x2, function(z) mean(z))

  # data frame,
  permutationTest2(xDF1, data2 = xDF2, colMeans)

  # data expression, data frame,
  permutationTest2((xDF1), data2 = (xDF2), colMeans)

  # data expression, matrix
  permutationTest2(as.matrix(xDF1), data2 = as.matrix(xDF2), colMeans)


  ### statistic expression
  # data by name  (refer to 'data')
  permutationTest2(x1, data2 = x2, mean(data))

  # data as expression, refer to 'data'
  permutationTest2((x1), data2 = (x2), mean(data))

  # data frame
  permutationTest2(xDF1, data2 = xDF2, mean(a))

  # data frame expression
  permutationTest2((xDF1), data2 = (xDF2), mean(a))

  # See if results reproduce
  temp <- .Last.value
  .Random.seed <- temp$seed
  all.equal(temp, eval(temp$call))

  source("~/resample/R/permutationTest2.R")
}
