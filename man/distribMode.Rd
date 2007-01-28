\name{distribMode}
\alias{distribMode}
\alias{betaMode}
\alias{cauchyMode}
\alias{chisqMode}
\alias{expMode}
\alias{fMode}
\alias{gammaMode}
\alias{normMode}
\alias{frechetMode}
\alias{gevMode}
\alias{ghMode}
\alias{gpdMode}
\alias{gumbelMode}
\alias{hypMode}
\alias{logisMode}
\alias{lnormMode}
\alias{nigMode}
\alias{stableMode}
\alias{symstbMode}
\alias{rweibullMode}
\alias{tMode}
\alias{unifMode}
\alias{weibullMode}
\alias{bernMode}
\alias{binomMode}
\alias{geomMode}
\alias{hyperMode}
\alias{nbinomMode}
\alias{poisMode}

\title{ Computing the Mode of Some Distributions }
\description{
These functions return the mode of the main distributions implemented in \R.
}

\usage{
## Continuous distributions
betaMode(shape1, shape2, ncp = 0) # Beta
cauchyMode(location = 0, ...) # Cauchy
chisqMode(df, ncp = 0) # Chisquare
expMode(...) # Exponentiel
fMode(df1, df2) # F
frechetMode(loc = 0, scale = 1, shape = 1, ...) # Fréchet (package 'evd')
gammaMode(shape, rate = 1, scale = 1/rate) # Gamma
normMode(mean = 0, ...) # Normal (Gaussian)
gevMode(loc = 1, scale = 1, shape = 1, ...) # Generalised Extreme Value (package 'evd')
ghMode(alpha = 1, beta = 0, delta = 1, mu = 0, 
       lambda = 1, ...) # Generalised Hyperbolic (package 'fBasics')
gpdMode(loc = 0, scale = 1, shape = 0, ...) # Generalised Pareto (package 'evd')
gumbelMode(loc = 0, ...) # Gumbel (package 'evd')
hypMode(alpha = 1, beta = 0, delta = 1, mu = 0, 
        pm = c(1, 2, 3, 4)) # Hyperbolic (package 'fBasics')
logisMode(location = 0, ...) # Logistic
lnormMode(meanlog = 0, sdlog = 1) # Lognormal
nigMode(alpha = 1, beta = 0, delta = 1, 
        mu = 0, ...) # Normal Inverse Gaussian (package 'fBasics')
stableMode(alpha, beta, gamma = 1, delta = 0, pm = 0, ...) # Stable (package 'fBasics')
symstbMode(...) # Symmetric stable (package 'fBasics')
rweibullMode(loc = 0, scale = 1, shape = 1, ...) # Negative Weibull (package 'evd')
tMode(df, ncp = 0) # T (Student)
unifMode(min = 0, max = 1) # Uniform
weibullMode(shape, scale = 1, ...) # Weibull

## Discrete distributions
bernMode(prob) # Bernoulli
binomMode(size, prob) # Binomial
geomMode(...) # Geometric
hyperMode(m, n, k, ...) # Hypergeometric
nbinomMode(size, prob, mu) # Negative Binomial
poisMode(lambda) # Poisson
}
\arguments{
  \item{shape1, shape2, ncp, location, df, df1, df2, loc, scale, shape,}{}
  \item{rate, mean, alpha, beta, delta, mu, lambda, pm, meanlog, sdlog,}{}
  \item{gamma, min, max, prob, size, m, n, k} {The different arguments are those of the corresponding distribution functions.}
  \item{...}{ further arguments, which will be ignored. }
}
\value{
A numeric value is returned, the (true) mode of the distribution.
}
\author{ Paul Poncet \email{paulponcet@yahoo.fr}, 
except for \code{hypMode} and \code{stableMode} written by Diethelm Wuertz, see package \pkg{fBasics}.}
\note{
Some functions like \code{normMode} or \code{cauchyMode}, which are related 
to symmetric distributions, are trivial, but are implemented for exhaustivity.
}
\seealso{ \code{\link{mlv}} for the estimation of the mode; 
the documentation of the related distributions \code{\link{Beta}}, \code{\link{GammaDist}}, etc. }
\examples{
layout(mat = matrix(1:2,1,2))

## Beta distribution
curve(dbeta(x, shape1 = 2, shape2 = 3.1), xlim = c(0,1), ylab = "Beta density")
M <- betaMode(shape1 = 2, shape2 = 3.1)
abline(v = M, col = 2)
mlv("beta", shape1 = 2, shape2 = 3.1)
 
## Lognormal distribution          
curve(dlnorm(x, meanlog = 3, sdlog = 1.1), xlim = c(0, 10), ylab = "Lognormal density")
M <- lnormMode(meanlog = 3, sdlog = 1.1)
abline(v = M, col = 2)
mlv("lnorm", meanlog = 3, sdlog = 1.1) 

## Poisson distribution
poisMode(lambda = 6)
poisMode(lambda = 6.1)
mlv("poisson", lambda = 6.1)

layout(mat = matrix(1,1,1)) 
}
\keyword{ univar }
\keyword{ distribution }
