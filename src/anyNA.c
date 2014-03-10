/* Copyright 2013 Google Inc. All rights reserved.

   Use of this source code is governed by a BSD-style
   license that can be found in the LICENSE file or at
   https://developers.google.com/open-source/licenses/bsd
*/

/*
   Determine whether an object contains a missing value.
   Return TRUE if there is a missing value, FALSE otherwise.
*/

#include <R_ext/Arith.h>      /* R_FINITE, ISNAN, ... */
#include <Rinternals.h>

/* This is modeled after the definitions for
   do_isna   in src/main/coerce.c
   do_logic3 in src/main/logic.c
   src/main/names.c gives mappings from 'is.na' to do_isna, 'any' to do_logic3.
*/


/* Check if an object has missing values. */
SEXP anyNA(SEXP x)
{
  int i, n;
  double *xDouble;
  int *xInt;

  n = length(x);

  if (n == 0)
    return Rf_ScalarLogical(FALSE);
  switch (TYPEOF(x)) {
    case REALSXP:
      xDouble = REAL(x);
      for (i = 0; i < n; i++)
        if (ISNAN(xDouble[i])) return Rf_ScalarLogical(TRUE);
      return Rf_ScalarLogical(FALSE);
    case INTSXP:
      xInt = INTEGER(x);
      for (i = 0; i < n; i++)
        if (xInt[i] == NA_INTEGER) return Rf_ScalarLogical(TRUE);
      return Rf_ScalarLogical(FALSE);
    case LGLSXP:
      xInt = LOGICAL(x);
      for (i = 0; i < n; i++)
        if (xInt[i] == NA_LOGICAL) return Rf_ScalarLogical(TRUE);
      return Rf_ScalarLogical(FALSE);
    case CPLXSXP:
      for (i = 0; i < n; i++)
        if (ISNAN(COMPLEX(x)[i].r) ||
           ISNAN(COMPLEX(x)[i].i)) return Rf_ScalarLogical(TRUE);
      return Rf_ScalarLogical(FALSE);
    case STRSXP:
      for (i = 0; i < n; i++)
        if (STRING_ELT(x, i) == NA_STRING) return Rf_ScalarLogical(TRUE);
      return Rf_ScalarLogical(FALSE);
    case LISTSXP:
      for (i = 0; i < n; i++){
        if (LOGICAL(anyNA(CAR(x)))[0]) return Rf_ScalarLogical(TRUE);
        x = CDR(x);
      }
      return Rf_ScalarLogical(FALSE);
    case VECSXP:
      for (i = 0; i < n; i++) {
        SEXP s = VECTOR_ELT(x, i);
        if (LOGICAL(anyNA(s))[0]) return Rf_ScalarLogical(TRUE);
      }
      return Rf_ScalarLogical(FALSE);
    case RAWSXP:
      /* no such thing as a raw NA */
      return Rf_ScalarLogical(FALSE);
    default:
      error("isNA applied to non-(list or vector) of type '%s'",
            type2char(TYPEOF(x)));
      return Rf_ScalarLogical(FALSE);
  }
}
