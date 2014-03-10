# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

cat0 <- function(...) cat(..., sep = "")
catn <- function(...) {
  cat(...)
  cat("\n")
}
cat0n <- function(...) {
  cat(..., sep = "")
  cat("\n")
}
