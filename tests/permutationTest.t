# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for permutationTest

# do.test("~/resample/tests/permutationTest.t")

{
  n <- 40
  myX <- 1:n
  myY <- seq(from = 1, by = 3, length = n) %% 16 + myX/6
  myDF <- data.frame(x = myX, y = myY)

  # Do resampling by hand, for comparison
  Robserved <- cor(myX, myY)
  Rreplicates <- numeric(999)
  set.seed(1)
  for(i in 1:999) {
    ii <- sample(n, replace = FALSE)
    Rreplicates[i] <- cor(myX[ii], myY)
  }
  RpValueL <- (1 + sum(Rreplicates <= Robserved)) / (999 + 1)
  RpValueG <- (1 + sum(Rreplicates >= Robserved)) / (999 + 1)
  RpValue <- 2 * min(RpValueL, RpValueG)
  if(FALSE){ # for work by hand
    c(Robserved, mean(Rreplicates)) # 0.540611004 0.002724527
    c(RpValueL, RpValueG, RpValue)  # 0.993 0.008 0.016
    qqnorm(Rreplicates)
    hist(Rreplicates); abline(v = Robserved, col="red")
  }
  compareFun <- function(r) {
    allTrue(all.equal(Robserved, unname(r$observed)),
            all.equal(Rreplicates, as.vector(r$replicates)),
            all.equal(RpValue, r$stats$PValue))
  }
  myCor <- function(x, y = myY) cor(x, y)
  TRUE
}


##### resampleColumns
{ # resampleColumns = character, statistic is expression of variables
  r <- permutationTest(myDF, cor(x, y), resampleColumns = "x", seed = 1)
  compareFun(r)
}

{ # resampleColumns = integer
  r <- permutationTest(myDF, cor(x, y), resampleColumns = 1, seed = 1)
  compareFun(r)
}

{ # statistic is expression of data
  r <- permutationTest(myDF, cor(myDF$x, myDF$y), resampleColumns = 1,
                       seed = 1)
  compareFun(r)
}

{ # statistic is inline function
  r <- permutationTest(myDF, function(X) cor(X$x, X$y), resampleColumns = 1,
                       seed=1)
  compareFun(r)
}

##### data is a vector, correlation with another vector
### statistic is function
{ # data by name, statistic is function by name
  r <- permutationTest(myX, myCor, seed = 1)
  compareFun(r)
}

{ # data expression
  r <- permutationTest((myX), myCor, seed = 1)
  compareFun(r)
}

{ # args.stat
  r <- permutationTest(myX, myCor, args.stat = list(y = myY), seed = 1)
  compareFun(r)
}

{ # inline function
  r <- permutationTest(myX, function(z) cor(z, myY), seed = 1)
  compareFun(r)
}

{ # statistic is expression
  r <- permutationTest(myX, cor(myX, myY), seed = 1)
  compareFun(r)
}

{ # statistic is expression, refer to "data" (probably doesn't work)
  r <- permutationTest(myX, cor(data, myY), seed = 1)
  compareFun(r)
}

{ # statistic is expression, data is expression
  r <- permutationTest((myX), cor(data, myY), seed = 1)
  compareFun(r)
}


### statistic is expression
{ # data by name
  r <- permutationTest(myX, myCor(myX), seed = 1)
  compareFun(r)
}

{ # data by name, but user referred to 'data'
  r <- permutationTest(myX, myCor(data), seed = 1)
  compareFun(r)
}

{ # data as expression, refer to 'data'
  r <- permutationTest((myX), myCor(data), seed = 1)
  compareFun(r)
}

{ # data frame, variable name
  r <- permutationTest(myDF, myCor(x), seed = 1)
  compareFun(r)
}

{ # data frame expression, variable name
  r <- permutationTest((myDF), myCor(x), seed = 1)
  compareFun(r)
}

{ # See if results reproduce
  .Random.seed <- r$seed
  all.equal(r, eval(r$call))
}

{ # One-sided test, and update
  r2 <- update(r, alternative = "less")
  all.equal(RpValueL, r2$stats$PValue)
}

{ # Multivariate statistic
  r <- permutationTest(myDF, cor, args.stat = list(y = 1:n))
  all.equal(2, nrow(r$stats))
}

{
  rm(n, myX, myY, myDF,
     Robserved, i, ii, Rreplicates, RpValueL, RpValueG, RpValue,
     compareFun, myCor, r, r2)
  TRUE
}
