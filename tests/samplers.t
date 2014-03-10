# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Test samp.bootstrap and samp.permute

# do.test("~/resample/tests/samplers.t")

{
  n <- 9
  B <- 3
  TRUE
}

{ # samp.bootstrap, dimensions of result
  allTrue(all.equal(c(n, B), dim(samp.bootstrap(n, B))),
          all.equal(c(7, B), dim(samp.bootstrap(n, B, size = 7))),
          all.equal(c(7, B), dim(samp.bootstrap(n, B, reduceSize = 2))))
}

{ # samp.bootstrap, values
  set.seed(0)
  temp <- matrix(sample(n, size = n * B, replace = TRUE), nrow = n)
  set.seed(0)
  all.equal(temp, samp.bootstrap(n, B))
}

{ # samp.bootstrap, values when size reduced
  set.seed(0)
  temp <- matrix(sample(n, size = 7 * B, replace = TRUE), nrow = 7)
  set.seed(0)
  all.equal(temp, samp.bootstrap(n, B, size = 7))
}

{ # samp.permute, dimensions of result
  allTrue(all.equal(c(n, B), dim(samp.permute(n, B))),
          all.equal(c(7, B), dim(samp.permute(n, B, size = 7))),
          all.equal(c(7, B), dim(samp.permute(n, B, reduceSize = 2))),
          all.equal(c(3, B),
            dim(samp.permute(n, B, groupSizes = c(3,n-3), returnGroup = 1))),
          all.equal(c(n-3, B),
            dim(samp.permute(n, B, groupSizes = c(3,n-3), returnGroup = 2))))
}

{ # samp.permute, values
  set.seed(0)
  standardPermutation <- matrix(NA, n, B)
  for(j in 1:B)
    standardPermutation[, j] <- sample(n, replace = FALSE)
  set.seed(0)
  all.equal(standardPermutation, samp.permute(n, B))
}

{ # samp.permute, values when using groups
  set.seed(0)
  all.equal(standardPermutation[1:3, ],
            samp.permute(n, B, groupSizes = c(3,n-3), returnGroup = 1))
}

{ # samp.permute, values for group 2
  set.seed(0)
  all.equal(standardPermutation[4:n, ],
            samp.permute(n, B, groupSizes = c(3,n-3), returnGroup = 2))
}

{ # samp.permute, values when size reduced
  set.seed(0)
  all.equal(standardPermutation[1:3, ],
            samp.permute(n, B, size = 3))
}

{
  rm(n, B, temp, standardPermutation)
  TRUE
}
