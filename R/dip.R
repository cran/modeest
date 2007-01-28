# Author : Martin Maechler <maechler@stat.math.ethz.ch> (author of package 'diptest'),
#          based on earlier code from Dario Ringach <dario@wotan.cns.nyu.edu>
dip <-
function(x,                  # sample (the data)
         full.result = FALSE,# if TRUE, a comprehensive result is returned
         debug = FALSE)      # for C code
{
  nx <- length(x)
  x <- sort(unname(x), method = "quick")
  r <- .C("diptst",
          x = as.double(x),
          n = nx,
          dip = double(1),
          lo.hi = integer(2),
          ifault = integer(1),
          gcm = integer(nx),
          lcm = integer(nx),
          mn = integer(nx),
          mj = integer(nx),
          debug = as.logical(debug),
          DUP = FALSE,
          PACKAGE = "modeest")[if (full.result) TRUE else "dip"]
  if (full.result)
    c(r, { u <- x[r$lo.hi] ; list(xl = u[1], xu = u[2]) })
  else r[[1]]
}
