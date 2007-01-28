# Author: Paul Poncet

parzen <-
function(x,                       # sample (the data)
         bw = NULL,               # bandwidth
         kernel = "gaussian",     # kernel used
         biau = FALSE,            # if FALSE, 'optim' is used
         par = shorth(x),         # initial value used in 'optim'
         optim.method = "BFGS",   # method used in 'optim'
         ...)
{
###################################################
# Mode estimator based on kernel density estimation
###################################################

  if (pmatch(tolower(kernel), "normal", no = 0)) {
    kernel <- "gaussian"
  } else kernel <- match.arg(tolower(kernel), c(.kernelList, "uniform")) # '.kernelList' is defined in 'K.R'

  if (kernel == "uniform") {
    if (is.null(bw)) bw <- 1/2
    if (bw <= 0 | bw >= 1) stop("argument 'bw' must belong to (0, 1)")    
    Fn <- ecdf(x)
    fn <- Fn(x+bw) - Fn(x-bw) # the estimate of the density is (Fn(x+bw) - Fn(x-bw))/(2*bw)
    ## Remark : optimization fails since fn is not regular enough
    ## (any point is viewed as a local maxima). The following is the only solution:
    M <- x[fn == max(fn)]
  
  } else {  
    
    ## Initialization
    nx <- length(x)
    if (is.null(bw)) bw <- bw.SJ(x)
    
    fn <-
    function(z)
    {
      z <- (rep(z,nx) - x)/rep(bw,nx)   #! enlever 'rep' : on ne peut pas appliquer 'fn' à un vecteur 'z' de longueur > 1 !
                                        #! cela permettrait d'accélérer la méthode 'biau'  
      z <- do.call(paste(".kernel.", kernel, sep = ""), list(z))$k
      return(sum(z))   #! à modifier
    }
  
    if (!biau) {
      maxi <- optim(par, fn, method = optim.method, control=list(fnscale=-1), ...)
      M <- maxi$par
      attr(M, "value") <- maxi$value
      attr(M, "counts") <- maxi$counts
      attr(M, "convergence") <- maxi$convergence
      attr(M, "message") <- maxi$message
    } else {
      f <- sapply(x,FUN=fn) #! pas très rapide... mieux vaut dans ce cas utiliser 'density' je pense
      M <- x[f == max(f)]
    }
  }
  
  ## Output
  return(M)   
}

#Parzen <- parzen

naive <-
function(x,
         bw = 1/2)
{  
  parzen(x = x, bw = bw, kernel = "uniform", biau = TRUE)
}

#pp <- list()
#for (i in 1:500) {
#  j <- 1/i
#  pp[[i]] <- paste("function(x) { (Fn(x + ", j, ") - Fn(x - ", j, "))/(2 * ", j, ") }", sep = "")
#}
#fn <- lapply(pp, function(x) eval(parse(text = x)))

#v <- seq(-10,10,0.1)
#plot(v, fn[[1]](v), col = cc[1])
#for (i in 2:500) {
#  lines(v, fn[[i]](v), col = cc[i])
#}
