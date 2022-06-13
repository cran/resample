# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for permutationTest2

# do.test("~/Rlang/resample/resample/tests/permutationTest2.t")

{
  set.seed(0)
  v1 <- rnorm(20)
  v2 <- rnorm(25) + .5
  v12 <- c(v1, v2)
  df1 <- data.frame(v = v1)
  df2 <- data.frame(v = v2)
  df12 <- rbind(df1, df2)
  treat <- rep(1:2, c(20, 25))
  df12t <- cbind(df12, treatdf = treat)

  # Do resampling using vanilla R, for comparison
  Robserved <- mean(v1) - mean(v2)
  Robserved1 <- mean(v1)
  Robserved2 <- mean(v2)
  Rreplicates1 <- numeric(999)
  Rreplicates2 <- numeric(999)
  RreplicatesTrim <- numeric(999)
  set.seed(1)
  for(i in 1:999) {
    ii <- sample(45, replace = FALSE)
    Rreplicates1[i] <- mean(v12[ii[1:20]])
    Rreplicates2[i] <- mean(v12[ii[21:45]])
    RreplicatesTrim[i] <- (mean(v12[ii[1:20]], trim = .25) -
                             mean(v12[ii[21:45]], trim = .25))
  }
  Rreplicates <- Rreplicates1 - Rreplicates2
  RpValue <- (1 + sum(Rreplicates <= Robserved)) / (999 + 1)
  if(FALSE){ # for work by hand
    Robserved         # -.585
    mean(Rreplicates) # .004
    RpValue           # .025
    qqnorm(Rreplicates)
    hist(Rreplicates); abline(v = Robserved, col = "red")
  }
  compareFun <- function(r) {
    r2 <- r$resultsBoth
    allTrue(all.equal(Robserved, unname(r$observed)),
            all.equal(Robserved1, unname(r2[[1]]$observed)),
            all.equal(Robserved2, unname(r2[[2]]$observed)),
            all.equal(Rreplicates1, as.vector(r2[[1]]$replicates)),
            all.equal(Rreplicates2, as.vector(r2[[2]]$replicates)),
            all.equal(RpValue, r$stats$PValue))
  }
  TRUE
}


##### treatment
### statistic is function
{
  # base case: data by name, statistic is function by name
  r <- permutationTest2(treatment = treat, v12, mean,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression
  r <- permutationTest2(treatment = treat, (v12), mean,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # args.stat
  r <- permutationTest2(treatment = treat, v12, mean,
                        args.stat = list(trim = .25),
                        alternative = "less", seed = 1, R = 999)
  all.equal(RreplicatesTrim, as.vector(r$replicates))
}

{
  # inline function
  r <- permutationTest2(treatment = treat, v12, function(z) mean(z),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame,
  r <- permutationTest2(treatment = treat, df12, colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression, data frame
  r <- permutationTest2(treatment = treat, df12, colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression, matrix
  r <- permutationTest2(treatment = treat, as.matrix(df12), colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame, treatment in data frame
  r <- permutationTest2(treatment = treatdf, df12t, function(x) mean(x[[1]]),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # seed not supplied
  set.seed(1)
  r <- permutationTest2(treatment = treat, v12, mean,
                        alternative = "less", R = 999)
  compareFun(r)
}


### statistic is expression
{
  # data by name
  r <- permutationTest2(treatment = treat, v12, mean(v12),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data by name, but user referred to 'data'
  r <- permutationTest2(treatment = treat, v12, mean(data),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- permutationTest2(treatment = treat, (v12), mean(data),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame
  r <- permutationTest2(treatment = treat, df12, mean(v),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame expression
  r <- permutationTest2(treatment = treat, (df12), mean(v),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # See if results reproduce
  .Random.seed <- r$seed
  all.equal(r, eval(r$call))
}


##### data and data2
### statistic is function
{
  # base case: data by name, statistic is function by name
  r <- permutationTest2(v1, data2 = v2, mean,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression
  r <- permutationTest2((v1), data2 = (v2), mean,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # args.stat
  r <- permutationTest2(v1, data2 = v2, mean, args.stat = list(trim = .25),
                        alternative = "less", seed = 1, R = 999)
  all.equal(RreplicatesTrim, as.vector(r$replicates))
}

{
  # inline function
  r <- permutationTest2(v1, data2 = v2, function(z) mean(z),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame,
  r <- permutationTest2(df1, data2 = df2, colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression, data frame,
  r <- permutationTest2((df1), data2 = (df2), colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data expression, matrix
  r <- permutationTest2(as.matrix(df1), data2 = as.matrix(df2), colMeans,
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}


### statistic is expression
{
  # data by name  (refer to 'data')
  r <- permutationTest2(v1, data2 = v2, mean(data),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- permutationTest2((v1), data2 = (v2), mean(data),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame
  r <- permutationTest2(df1, data2 = df2, mean(v),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # data frame expression
  r <- permutationTest2((df1), data2 = (df2), mean(v),
                        alternative = "less", seed = 1, R = 999)
  compareFun(r)
}

{
  # seed not supplied, data and data2
  set.seed(1)
  r <- permutationTest2(v1, data2 = v2, mean,
                        alternative = "less", R = 999)
  compareFun(r)
}

{
  # See if results reproduce
  .Random.seed <- r$seed
  all.equal(r, eval(r$call))
}

{
  rm(v1, v2, v12, df1, df2, df12, treat, df12t,
     Robserved, Robserved1, Robserved2,
     Rreplicates1, Rreplicates2,
     i, ii, Rreplicates, RreplicatesTrim, compareFun, r)
  TRUE
}
