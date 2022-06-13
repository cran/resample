\name{IfElse}
\alias{IfElse}
\title{Conditional Data Selection}
\concept{splus}
\usage{IfElse(test, yes, no, ...)}
\description{This is equivalent to \code{{if(test) yes else no}}.
  The advantages of using this function are better formatting, and a more
  natural syntax when the result is being assigned; see examples below.

  With 5 arguments, this is equivalent to
  \code{{if(test1) yes else if(test2) u else v}}
  (where arguments are given by name, not position).
}
\arguments{
\item{test}{logical value; if \code{TRUE} return \code{yes}.}

\item{yes}{any object; this is returned if \code{test} is \code{TRUE}.}

\item{no}{normally any object; this is returned if \code{test} is \code{FALSE}.
          If there are more than three arguments this should be logical.}
\item{...}{there should be 3, 5, 7, etc. arguments to this function;
arguments 1, 3, 5, etc. should be logical values; the other arguments
(even numbered, and last) are objects that may be returned.}
}

\details{
\code{test} should be a scalar logical, and only one of \code{yes} or
\code{no} is evaluated, depending on whether \code{test = TRUE} or
\code{test = FALSE}, and \code{yes} and \code{no} may be any objects.
In contrast, for
\code{ifelse}, test is normally a vector, both \code{yes} and \code{no}
are evaluated, even if
not used, and \code{yes} and \code{no} are vectors the same length as
\code{test}.
}

\value{
  with three arguments, one of \code{yes} or \code{no}.
  With k arguments, one of arguments 2, 4, ..., k-1, k.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{https://www.timhesterberg.net/bootstrap-and-resampling}
}
\seealso{
\code{\link[base]{ifelse}}, \code{\link[base]{if}}.}

\examples{
IfElse(TRUE, "cat", "dog")
IfElse(FALSE, "one", TRUE, "two", "three")
IfElse(FALSE, "one", FALSE, "two", "three")
}
\keyword{utilities}
