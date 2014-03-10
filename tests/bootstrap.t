# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for bootstrap

# do.test("~/resample/tests/bootstrap.t")

{
  set.seed(0)
  v1 <- rnorm(20)
  df1 <- data.frame(v = v1)

  # Do resampling by hand, for comparison
  Robserved1 <- mean(v1)
  Rreplicates1 <- numeric(1000)
  RreplicatesTrim1 <- numeric(1000)
  set.seed(1)
  for(i in 1:1000) {
    ii <- sample(20, replace = TRUE)
    Rreplicates1[i] <- mean(v1[ii])
    RreplicatesTrim1[i] <- mean(v1[ii], trim = .25)
  }

  if(FALSE){ # for work by hand
    Robserved         # -.585
    mean(Rreplicates) # -.590
    qqnorm(Rreplicates)
    hist(Rreplicates); abline(v = Robserved, col="red")
  }
  compareFun <- function(r) {
    r2 <- r$resultsBoth
    allTrue(all.equal(Robserved1, unname(r$observed)),
            all.equal(Rreplicates1, as.vector(r$replicates)))
  }
  TRUE
}

# This should test most combinations of:
#    statistic is function (inline or by name), or expression
#    data is given by name, or an expression
#    data is a data frame (special handling for expression), or not.
#    args.stat


### statistic is function
{
  # base case: data by name, statistic is function by name
  r <- bootstrap(v1, mean, seed = 1)
  compareFun(r)
}

{
  # data expression
  r <- bootstrap((v1), mean, seed = 1)
  compareFun(r)
}

{
  # args.stat
  r <- bootstrap(v1, mean, args.stat = list(trim = .25), seed = 1)
  all.equal(RreplicatesTrim1, as.vector(r$replicates))
}

{
  # inline function
  r <- bootstrap(v1, function(z) mean(z), seed = 1)
  compareFun(r)
}

{
  # data frame,
  r <- bootstrap(df1, colMeans, seed = 1)
  compareFun(r)
}

{
  # data expression, data frame
  r <- bootstrap(df1, colMeans, seed = 1)
  compareFun(r)
}

{
  # data expression, matrix
  r <- bootstrap(as.matrix(df1), colMeans, seed = 1)
  compareFun(r)
}


### statistic is expression
{
  # data by name
  r <- bootstrap(v1, mean(v1), seed = 1)
  compareFun(r)
}

{
  # data by name, but user referred to 'data'
  r <- bootstrap(v1, mean(data), seed = 1)
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- bootstrap((v1), mean(data), seed = 1)
  compareFun(r)
}

{
  # data frame
  r <- bootstrap(df1, mean(v), seed = 1)
  compareFun(r)
}

{
  # data frame expression
  r <- bootstrap((df1), mean(v), seed = 1)
  compareFun(r)
}

{
  # See if results reproduce
  .Random.seed <- r$seed
  all.equal(r, eval(r$call))
}

{
  rm(v1, df1, Robserved1, Rreplicates1, RreplicatesTrim1, i, ii, compareFun, r)
  TRUE
}
