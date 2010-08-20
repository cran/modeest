# Author: Paul Poncet

vieu <-
function(x,                       # sample (the data)
         bw = NULL,               # bandwidth
         kernel = "gaussian",     # kernel used
         abc = FALSE,             # if FALSE, 'optim' is used
         ...)
{
###################################################
# Mode estimator based on the kernel estimate 
#  of the derivative of the density 
###################################################

  if (pmatch(tolower(kernel), "normal", no = 0)) {
    kernel <- "gaussian"
  } else kernel <- match.arg(tolower(kernel), c(.kernelList, "uniform")) # '.kernelList' is defined in 'K.R'
    
  ## Initialization
  nx <- length(x)
  if (is.null(bw)) bw <- bw.SJ(x)
  
  fn <-
  function(z)
  {
    z <- (rep(z,nx) - x)/rep(bw,nx)   #! enlever 'rep' : on ne peut pas appliquer 'fn' à un vecteur 'z' de longueur > 1 !
                                      #! cela permettrait d'accelerer la methode 'abc'  
    z <- do.call(paste(".kernel.d", kernel, sep = ""), list(z))$k
    return(sum(z))   #! à modifier
  }

  if (!abc) {
    r <- uniroot(f = fn, interval = c(min(x), max(x)), ...)
    M <- r$root
  } else {
    f <- abs(sapply(x,FUN=fn)) #! pas tres rapide... mieux vaut dans ce cas utiliser 'density' je pense
    M <- x[f == min(f)]
  }
  
  ## Output
  return(M)   
}
