\name{naive}
%\alias{mlv.chernoff}
%\alias{mlv.naive}
\alias{naive}

\title{ The Chernoff Mode Estimator }
\description{
  This estimator, also called the *naive* mode estimator, is defined as the 
  center of the interval of given length containing the most observations. 
  It is identical to Parzen's kernel mode estimator, when the kernel is chosen 
  to be the uniform kernel.
}
\usage{
naive(x, 
      bw = 1/2)
}

\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. Should belong to (0, 1). See below. }
}
\value{
A numeric vector is returned, the mode estimate, which is the center of the 
interval of length \code{2*bw} containing the most observations.
}
\note{
The user should preferentially call \code{naive} through 
 \code{mlv(x, method = "naive", bw)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{ 
    \item Chernoff H. (1964).
    Estimation of the mode.
    \emph{Ann. Inst. Statist. Math.}, \bold{16}:31-41.
    
    \item Leclerc J. (1997).
    Comportement limite fort de deux estimateurs du mode : le shorth et l'estimateur naïf.
    \emph{C. R. Acad. Sci. Paris, Série I}, \bold{325}(11):1207-1210.
    
    \item Leclerc J. (2000).
    Strong limiting behavior of two estimates of the mode: the shorth and the naive estimator. 
    \emph{Statistics and Decisions},\bold{18}(4).
  }
}

\author{ Paul Poncet \email{paulponcet@yahoo.fr} }
\seealso{ \code{\link{mlv}} for general mode estimation;
\code{\link{parzen}} for Parzen's kernel mode estimation }
\examples{
# Unimodal distribution
x <- rf(10000, df1 = 40, df2 = 30)
## True mode
fMode(df1 = 40, df2 = 30)
## Estimate of the mode
mean(naive(x, bw = 1/4))
M <- mlv(x, method = "naive", bw = 1/4)
print(M)
plot(M, xlim = c(0,2))
}
\keyword{ univar }
