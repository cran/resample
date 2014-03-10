# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

anyNA <- function(x)
  UseMethod("anyNA")

anyNA.default <- function(x){
  # Return TRUE if there are missing values.
  if (is.list(x)) {
    for (i in seq_along(x))
      if (anyNA(x[[i]]))
        return(TRUE)
    return(FALSE)
  }
  if(isS4(x))
    return(any(is.na(x)))
  .Call(anyNA_sym, x) # mapped in NAMESPACE from anyNA
}

# TODO: Would like to handle everything in .Call, including dispatching to
# S3 and S4 methods. Could do that with an implicit generic, but that
# requires that the function be .Internal.  Defer for now; might be easiest
# to do if and when this incorporated into base R.
