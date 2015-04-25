# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Test samp.bootstrap and samp.permute

# do.test("~/resample/tests/samplers.t")

{
  n <- 9
  R <- 3
  TRUE
}

{ # samp.bootstrap, dimensions of result
  allTrue(all.equal(c(n, R), dim(samp.bootstrap(n, R))),
          all.equal(c(7, R), dim(samp.bootstrap(n, R, size = 7))),
          all.equal(c(7, R), dim(samp.bootstrap(n, R, reduceSize = 2))))
}

{ # samp.bootstrap, values
  set.seed(0)
  temp <- matrix(sample(n, size = n * R, replace = TRUE), nrow = n)
  set.seed(0)
  all.equal(temp, samp.bootstrap(n, R))
}

{ # samp.bootstrap, values when size reduced
  set.seed(0)
  temp <- matrix(sample(n, size = 7 * R, replace = TRUE), nrow = 7)
  set.seed(0)
  all.equal(temp, samp.bootstrap(n, R, size = 7))
}

{ # samp.permute, dimensions of result
  allTrue(all.equal(c(n, R), dim(samp.permute(n, R))),
          all.equal(c(7, R), dim(samp.permute(n, R, size = 7))),
          all.equal(c(7, R), dim(samp.permute(n, R, reduceSize = 2))),
          all.equal(c(3, R),
            dim(samp.permute(n, R, groupSizes = c(3, n-3), returnGroup = 1))),
          all.equal(c(n-3, R),
            dim(samp.permute(n, R, groupSizes = c(3, n-3), returnGroup = 2))))
}

{ # samp.permute, values
  set.seed(0)
  standardPermutation <- matrix(NA, n, R)
  for(j in 1:R)
    standardPermutation[, j] <- sample(n, replace = FALSE)
  set.seed(0)
  all.equal(standardPermutation, samp.permute(n, R))
}

{ # samp.permute, values when using groups
  set.seed(0)
  all.equal(standardPermutation[1:3, ],
            samp.permute(n, R, groupSizes = c(3, n-3), returnGroup = 1))
}

{ # samp.permute, values for group 2
  set.seed(0)
  all.equal(standardPermutation[4:n, ],
            samp.permute(n, R, groupSizes = c(3, n-3), returnGroup = 2))
}

{ # samp.permute, values when size reduced
  set.seed(0)
  all.equal(standardPermutation[1:3, ],
            samp.permute(n, R, size = 3))
}

{
  rm(n, R, temp, standardPermutation)
  TRUE
}
