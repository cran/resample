# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

# Tests for bootstrap

# do.test("~/resample/resample/tests/confidence.t")

{
  set.seed(0)
  v1 <- rnorm(20)
  df1 <- data.frame(v = v1)
  df2 <- data.frame(v = v1, w = 2 * v1)
  r1 <- bootstrap(v1, mean, seed=1, R=1000)
  r2 <- bootstrap(df2, colMeans, seed=1, R=1000)
  r3 <- bootstrap(v1, c(mean(v1), sd(v1)), seed=1, R=1000)
  CIt1 <- CI.t(r1)
  CIt2 <- CI.t(r2)
  CIp1 <- CI.percentile(r1)
  CIp2 <- CI.percentile(r2)
  CIbca1 <- CI.bca(r1)
  CIbca2 <- CI.bca(r2)
  CIbootT <- CI.bootstrapT(r3)
  TRUE
}

### Dimensions and dimnames
{ # dimensions
  allTrue(all.equal(c(1,2), dim(CIt1)),
          all.equal(c(2,2), dim(CIt2)),
          all.equal(c(1,2), dim(CIp1)),
          all.equal(c(2,2), dim(CIp2)),
          all.equal(c(1,2), dim(CIbca1)),
          all.equal(c(2,2), dim(CIbca2)),
          all.equal(c(1,2), dim(CIbootT)))
}

{ # dimnames
  dn1 <- list(names(r1$observed), c("2.5%", "97.5%"))
  dn2 <- list(names(r2$observed), c("2.5%", "97.5%"))
  dn3 <- list(names(r3$observed)[1], c("2.5%", "97.5%"))
  allTrue(identical(dn1, dimnames(CIt1)),
          identical(dn2, dimnames(CIt2)),
          identical(dn1, dimnames(CIp1)),
          identical(dn2, dimnames(CIp2)),
          identical(dn1, dimnames(CIbca1)),
          identical(dn2, dimnames(CIbca2)),
          identical(dn3, dimnames(CIbootT)))
}

{ # consistent values
  allTrue(all.equal(CIt1[1,], CIt2[1,]),
          all.equal(CIp1[1,], CIp2[1,]),
          all.equal(CIbca1[1,], CIbca2[1,]),
          all.equal(CIt1[1,], CIt2[2,]/2),
          all.equal(CIp1[1,], CIp2[2,]/2),
          all.equal(CIbca1[1,], CIbca2[2,]/2))
}

{ # approximately equal to t interval
  tCI <- mean(v1) + qt(c(.025,.975), length(v1)-1) * sd(v1)/sqrt(length(v1))
  allTrue(all.equal(tCI, c(CIt1), tol = .1),
          all.equal(tCI, c(CIp1), tol = .1),
          all.equal(tCI, c(CIbca1), tol = .1),
          all.equal(tCI, c(CIbootT), tol = .15))
}

{
  rm(v1, df1, df2, r1, r2, r3, CIt1, CIt2, CIp1, CIp2, CIbca1, CIbca2,
     CIbootT, dn1, dn2, dn3, tCI)
  TRUE
}
