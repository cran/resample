# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for bootstrap2

# do.test("~/resample/resample/tests/bootstrap2.t")

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
  Rreplicates1 <- numeric(1000)
  Rreplicates2 <- numeric(1000)
  RreplicatesTrim1 <- numeric(1000)
  RreplicatesTrim2 <- numeric(1000)
  set.seed(1)
  for(i in 1:1000) {
    ii <- sample(20, replace = TRUE)
    Rreplicates1[i] <- mean(v1[ii])
    RreplicatesTrim1[i] <- mean(v1[ii], trim = .25)
  }
  for(i in 1:1000) {
    ii <- sample(25, replace = TRUE)
    Rreplicates2[i] <- mean(v2[ii])
    RreplicatesTrim2[i] <- mean(v2[ii], trim = .25)
  }
  Rreplicates <- Rreplicates1 - Rreplicates2
  RreplicatesTrim <- RreplicatesTrim1 - RreplicatesTrim2

  if(FALSE){ # for work by hand
    Robserved         # -.585
    mean(Rreplicates) # -.590
    qqnorm(Rreplicates)
    hist(Rreplicates); abline(v = Robserved, col = "red")
  }
  compareFun <- function(r) {
    r2 <- r$resultsBoth
    allTrue(all.equal(Robserved, unname(r$observed)),
            all.equal(Robserved1, unname(r2[[1]]$observed)),
            all.equal(Robserved2, unname(r2[[2]]$observed)),
            all.equal(Rreplicates1, as.vector(r2[[1]]$replicates)),
            all.equal(Rreplicates2, as.vector(r2[[2]]$replicates)))
  }
  TRUE
}

# This should test most combinations of:
#    statistic is function (inline or by name), or expression
#    data2 or treatment
#    data is given by name, or an expression
#    data is a data frame (special handling for expression & treatment), or not.
#    args.stat



##### treatment
### statistic is function
{
  # base case: data by name, statistic is function by name
  r <- bootstrap2(treatment = treat, v12, mean, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression
  r <- bootstrap2(treatment = treat, (v12), mean, seed = 1, R = 1000)
  compareFun(r)
}

{
  # args.stat
  r <- bootstrap2(treatment = treat, v12, mean, args.stat = list(trim = .25),
                  seed = 1, R = 1000)
  all.equal(RreplicatesTrim, as.vector(r$replicates))
}

{
  # inline function
  r <- bootstrap2(treatment = treat, v12, function(z) mean(z),
                  seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame,
  r <- bootstrap2(treatment = treat, df12, colMeans, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression, data frame
  r <- bootstrap2(treatment = treat, df12, colMeans, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression, matrix
  r <- bootstrap2(treatment = treat, as.matrix(df12), colMeans,
                  seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame, treatment in data frame
  r <- bootstrap2(treatment = treatdf, df12t, function(x) mean(x[[1]]),
                  seed = 1, R = 1000)
  compareFun(r)
}


### statistic is expression
{
  # data by name
  r <- bootstrap2(treatment = treat, v12, mean(v12), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data by name, but user referred to 'data'
  r <- bootstrap2(treatment = treat, v12, mean(data), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- bootstrap2(treatment = treat, (v12), mean(data), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame
  r <- bootstrap2(treatment = treat, df12, mean(v), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame expression
  r <- bootstrap2(treatment = treat, (df12), mean(v), seed = 1, R = 1000)
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
  r <- bootstrap2(v1, data2 = v2, mean, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression
  r <- bootstrap2((v1), data2 = (v2), mean, seed = 1, R = 1000)
  compareFun(r)
}

{
  # args.stat
  r <- bootstrap2(v1, data2 = v2, mean, args.stat = list(trim = .25),
                  seed = 1, R = 1000)
  all.equal(RreplicatesTrim, as.vector(r$replicates))
}

{
  # inline function
  r <- bootstrap2(v1, data2 = v2, function(z) mean(z), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame,
  r <- bootstrap2(df1, data2 = df2, colMeans, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression, data frame,
  r <- bootstrap2((df1), data2 = (df2), colMeans, seed = 1, R = 1000)
  compareFun(r)
}

{
  # data expression, matrix
  r <- bootstrap2(as.matrix(df1), data2 = as.matrix(df2), colMeans,
                  seed = 1, R = 1000)
  compareFun(r)
}


### statistic is expression
{
  # data by name  (refer to 'data')
  r <- bootstrap2(v1, data2 = v2, mean(data), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- bootstrap2((v1), data2 = (v2), mean(data), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame
  r <- bootstrap2(df1, data2 = df2, mean(v), seed = 1, R = 1000)
  compareFun(r)
}

{
  # data frame expression
  r <- bootstrap2((df1), data2 = (df2), mean(v), seed = 1, R = 1000)
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
     RreplicatesTrim1, RreplicatesTrim2,
     i, ii, Rreplicates, RreplicatesTrim, compareFun, r)
  TRUE
}
