# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for jackknife

# do.test("~/resample/resample/tests/jackknife.t")

{
  set.seed(0)
  v1 <- rnorm(20)
  df1 <- data.frame(v = v1)

  # Do resampling using vanilla R, for comparison
  Robserved1 <- mean(v1)
  RobservedTrim1 <- mean(v1, trim = .25)
  Rreplicates1 <- numeric(20)
  RreplicatesTrim1 <- numeric(20)
  for(i in 1:20) {
    ii <- -i
    Rreplicates1[i] <- mean(v1[ii])
    RreplicatesTrim1[i] <- mean(v1[ii], trim = .25)
  }

  if(FALSE){ # for work by hand
    Robserved1         # -0.001778674
    mean(Rreplicates1) # -0.001778674
    qqnorm(Rreplicates1)
    hist(Rreplicates1); abline(v = Robserved1, col = "red")
  }
  compareFun <- function(r) {
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
  r <- jackknife(v1, mean)
  compareFun(r)
}

{
  # data expression
  r <- jackknife((v1), mean)
  compareFun(r)
}

{
  # args.stat
  r <- jackknife(v1, mean, args.stat = list(trim = .25))
  all.equal(RreplicatesTrim1, as.vector(r$replicates))
}

{
  # inline function
  r <- jackknife(v1, function(z) mean(z))
  compareFun(r)
}

{
  # data frame,
  r <- jackknife(df1, colMeans)
  compareFun(r)
}

{
  # data expression, data frame
  r <- jackknife(df1, colMeans)
  compareFun(r)
}

{
  # data expression, matrix
  r <- jackknife(as.matrix(df1), colMeans)
  compareFun(r)
}


### statistic is expression
{
  # data by name
  r <- jackknife(v1, mean(v1))
  compareFun(r)
}

{
  # data by name, but user referred to 'data'
  r <- jackknife(v1, mean(data))
  compareFun(r)
}

{
  # data as expression, refer to 'data'
  r <- jackknife((v1), mean(data))
  compareFun(r)
}

{
  # data frame
  r <- jackknife(df1, mean(v))
  compareFun(r)
}

{
  # data frame expression
  r <- jackknife((df1), mean(v))
  compareFun(r)
}

{
  # See if results reproduce
  .Random.seed <- r$seed
  all.equal(r, eval(r$call))
}

### Basic computations
{ # jackknife bias for mean
  r <- jackknife(v1, mean)
  all.equal(r$stats$Bias, 0)
}
{ # jackknife SE for mean
  all.equal(r$stats$SE, sd(v1)/sqrt(20))
}
{ # jackknife bias for trimmed mean
  r <- jackknife(v1, mean, args.stat = list(trim = .25))
  all.equal(r$stats$Bias, 19 * (mean(RreplicatesTrim1) - RobservedTrim1))
}


{
  rm(v1, df1, Robserved1, Rreplicates1, RreplicatesTrim1, i, ii, compareFun, r)
  TRUE
}
