\name{dip}
\alias{dip}
\alias{DIP}

\title{ Compute Hartigan's DIP Test Statistic for Unimodality }
\description{
  Computes Hartigan's DIP test statistic for testing unimodality,
  and additionally the modal interval. This function is taken from Martin Maechler's \pkg{diptest} package.
}
\usage{
dip(x, 
    full.result = FALSE, 
    debug = FALSE)
}
\arguments{
  \item{x}{numeric. Vector of observations. }
  \item{full.result}{ logical. If \code{TRUE} returns the full result
    list, see below.}
  \item{debug}{ logical. If \code{TRUE}, some tracing information is printed
    (from the C routine).}
}
\value{
  Depending on \code{full.result} either a number, the dip statistic, or
  a list with components
  \item{x}{the sorted \code{\link{unname}()}d data.}
  \item{n}{\code{length(x)}.}
  \item{dip}{the dip statistic}
  \item{lo.hi}{indices into \code{x} for lower and higher end of modal interval}
  \item{xl, xu}{lower and upper end of modal interval}
  \item{gcm, lcm}{(last used) indices for \bold{g}reatest \bold{c}onvex
    \bold{m}inorant and the \bold{l}east \bold{c}oncave \bold{m}ajorant.}
  \item{mn, mj}{index vectors of length \code{n} for the GC minorant and
    the LC majorant respectively.}
}
\note{
  For \eqn{n \le 3}{n <= 3} where \code{n <- length(x)}, the dip
  statistic is always zero, i.e., there's no possible dip test.

  Yong Lu \email{lyongu+@cs.cmu.edu} found in Oct 2003 that the code
  wasn't giving symmetric results for mirrored data (and was giving
  results of almost 1, and then found the reason, a misplaced ")" in the
  original Fortran code.  This bug has been corrected for diptest
  version 0.25-0.
}
\references{
  \itemize{
    \item Hartigan J.A. and Hartigan P.M. (1985).
    The Dip Test of Unimodality.
    \emph{Ann. Statist.}, \bold{13}:70-84.
    
    \item Hartigan P.M. (1985).
    Computation of the Dip Statistic to Test for Unimodality.
    \emph{Appl. Statist. (JRSS C)}, \bold{34}:320-325.
  }
}
\author{Martin Maechler \email{maechler@stat.math.ethz.ch}, based on
  earlier code from Dario Ringach \email{dario@wotan.cns.nyu.edu}
}
\seealso{\code{\link{isoreg}} for isotonic regression.}
\examples{
data(precip)
plot(density(precip))
str(dip(precip, full = TRUE, debug = TRUE))
}
\keyword{htest}
\keyword{distribution}
