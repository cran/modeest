\name{Chernoff}
%\alias{Chernoff}
%\alias{chernoff}
\alias{dchern}
\alias{pchern}
\alias{qchern}
\alias{rchern}

\title{The Chernoff Distribution}
\description{
  Density, distribution function, quantile function and random generation 
  for the Chernoff distribution.
}

\usage{
dchern(x, log = FALSE)
pchern(q, lower.tail = TRUE, log.p = FALSE)
qchern(p, lower.tail = TRUE, log.p = FALSE)
rchern(n)
}
\arguments{
  \item{x, q}{ vector of quantiles. }
  \item{p}{ vector of probabilities. }
  \item{n}{ number of observations. }
  \item{log, log.p}{logical. If \code{TRUE}, probabilities p are given as log(p).}
  \item{lower.tail}{logical. If \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.}
}
\value{
  \code{dchern} gives the density, \code{pchern} gives the distribution
  function, \code{qchern} gives the quantile function and \code{rchern}
  generates random deviates.
}
\details{
  The Chernoff distribution is the law of the random variable
  \deqn{Z =  \arg\max_{t} (B(t) - t^2)}{Z = argmax_t (B(t) - t^2)}
  where \eqn{B(t)}{B(t)} is two-sided Brownian motion starting from zero.
  Up to a location and a scale parameter, this distribution is the limit 
  distribution of the Chernoff mode estimator (also called the \code{\link{naive}} 
  mode estimator), see Chernoff (1964).
}

\section{Warning}{
Presently, \code{pchern}, \code{qchern}, and \code{rchern} are *very* slow.
}

\source{
  \code{dchern}, \code{pchern} and \code{qchern} are calculated
  from the C code due to Groeneboom and Wellner (2001).
}
\references{
  \itemize{
    \item Chernoff H. (1964).
    Estimation of the mode.
    \emph{Ann. Inst. Statist. Math.}, \bold{16}:31-41.
  
    \item Groeneboom P.and Wellner J.A. (2001).
    Computing Chernoff's distribution.
    \emph{J. Comput. Graph. Statist.}, \bold{10}:388-400.
  }
}
\author{ Groeneboom and Wellner for the original C code, \cr
Paul Poncet \email{paulponcet@yahoo.fr} for the R translation.
}
\seealso{
  \code{\link{mlv}} for general mode estimation;
  \code{\link{naive}} for the Chernoff mode estimate;
  \code{\link{distribMode}} for computation of the mode of the classical distributions
}
\examples{
dchern(seq(0.1,0.9,0.1))
curve(dchern(x), xlim = c(-2,2))
}
\keyword{distribution}
