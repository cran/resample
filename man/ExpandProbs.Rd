\name{ExpandProbs}
\alias{ExpandProbs}
\title{
  Calculate modified probabilities for more accurate confidence intervals
}
\description{
  Compute modified quantiles levels, for more accurate confidence
  intervals. Using these levels gives sider intervals, with closer to
  desired coverage.
}
\usage{
ExpandProbs(probs, n)
}
\arguments{
  \item{probs}{
    vector of numerical values between 0 and 1.
}
  \item{n}{
    number of observations.
}
}
\details{
Bootstrap percentile confidence interval for a sample mean correspond roughly to
\deqn{\bar x \pm z_\alpha \hat\sigma}{xbar +- z_alpha sigmaHat}
instead of
\deqn{\bar x \pm t_{\alpha,n-1} s}{xbar +- t_alpha,n-1 s}
where
\deqn{\hat\sigma = \sqrt{(n-1)/n s}}{sigmaHat = sqrt((n-1)/n) s}
is like s but computed using a divisor of n instead of n-1.
Similarly for other statistics, the bootstrap percentile interval
is too narrow, typically by roughly the same proportion.

This function finds modified probability levels probs2, such that
\deqn{z_{\mbox{probs2}} \sqrt{(n-1)/n} = t_{\mbox{probs}, n-1}}%
     {z_probs2 sqrt((n-1)/n) = t_probs,n-1}
  so that for symmetric data, the bootstrap percentile interval
  approximately matches the usual $t$ confidence interval.
}
\value{
A vector like \code{probs}, but with values closer to 0 and 1.
}
\references{
  This discusses the expanded percentile interval:
  Hesterberg, Tim (2014), What Teachers Should Know about the
  Bootstrap: Resampling in the Undergraduate Statistics Curriculum,
  \url{http://arxiv.org/abs/1411.5279}.
}
\author{
  Tim Hesterberg \email{timhesterberg@gmail.com},\cr
  \url{http://www.timhesterberg.net/bootstrap}
}
\seealso{
  \code{\link{CI.percentile}},
  \code{\link{CI.bca}},
}
\examples{
probs <- c(0.025, 0.975)
n <- c(5, 10, 20, 40, 100, 200, 1000)
outer(probs, n, ExpandProbs)
}
\keyword{nonparametric}
\keyword{htest}
