##' Create psi_func_cached object using cached numerical integration for
##' E... slots.
##'
##' The E... slots will not be fully functional: they just return the
##' value for the current defaults and ignore their arguments.
##'
##' @title psiFuncCached constructor
##' @param rho rho-function
##' @param psi psi-function
##' @param wgt wgt-function
##' @param Dwgt derivative of weight function
##' @param Dpsi derivative of psi
##' @param name descriptor of this function family
##' @param ... default values for tuning constants
##' @return psi_func_cached-class object
##' @export
psiFuncCached <- function(rho,psi,wgt,Dwgt,Dpsi,name=NULL, ...) {
    lent <- length(dotsargs <- list(...))
    ## '...'  must contain all tuning parameters and their defaults:
    stopifnot(length(nt <- names(dotsargs)) == lent,
              all(nchar(nt)) >= 1)
    if(lent >= 1) {
        ## rho, psi,... checking: must have argument names
        argn <- c("x", nt)
        for(fnam in list("rho", "psi", "wgt", "Dwgt", "Dpsi")) {
            f <- get(fnam, inherits = FALSE)
            ef <- environment(f)
            nf <- names(ff <- formals(f)) # "x" and "k" for Huber's
            if(!identical(nf, argn))
                stop("arguments of function '",fnam,"' are (",
                     paste(nf,  collapse=","),") but should be (",
                     paste(argn,collapse=","),").")
            
            formals(f)[-1] <- dotsargs
            environment(f) <- ef
            assign(fnam, f, inherits = FALSE)
        }
    }

    Erho.val <- integrate(function(x) rho(x)*dnorm(x),-Inf, Inf,
                          rel.tol = .Machine$double.eps^0.5)$value
    Epsi2.val <- integrate(function(x) psi(x)^2*dnorm(x),-Inf, Inf,
                           rel.tol = .Machine$double.eps^0.5)$value
    EDpsi.val <- integrate(function(x) Dpsi(x)*dnorm(x),-Inf, Inf,
                           rel.tol = .Machine$double.eps^0.5)$value
    
    new("psi_func_cached",
        rho = new("functionX", rho),
        psi = new("functionX", psi),
        wgt = new("functionX", wgt),
        Dpsi= new("functionX", Dpsi),
        Dwgt= new("functionX", Dwgt),
        ## tNams = if(lent) nt else character(0),
        tDefs = if(lent) unlist(dotsargs) else numeric(0),
        Erho= Erho <- new("functionXal", function(arg=1) rep(Erho.val, length(arg))),
        Epsi2= Epsi2 <- new("functionXal", function(arg=1) rep(Epsi2.val, length(arg))),
        EDpsi= EDpsi <- new("functionXal", function(arg=1) rep(EDpsi.val, length(arg))),
        name= name
        )
}

##' Change the default arguments for a psi_func_cached object
##'
##' @title Change default arguments
##' @param ... arguments to change
##' @export
setMethod("chgDefaults", signature("psi_func_cached"),
          function(object, ...) {
              ##cat("~~~~ chgDefaults of psi_func_cached ~~~~~\n")
              lent <- length(dotsargs <- list(...))
              ## '...'  must contain all tuning parameters and their defaults:
              stopifnot(length(nt <- names(dotsargs)) == lent,
                        all(nchar(nt)) >= 1)
              if(lent >= 1) {
                  ## rho "..." must conform to rho, etc:
                  nf <- names(ff <- formals(object@rho))
                  if(!identical(nf[-1], nt))
                     stop("invalid tuning parameter names: ",
                          paste(nt,    collapse=",")," instead of ",
                          paste(nf[-1],collapse=","),".")

                  for(fnam in list("rho", "psi", "wgt", "Dwgt", "Dpsi")) {
                      f <- slot(object, fnam)
                      ef <- environment(f)
                      formals(f)[-1] <- dotsargs
                      environment(f) <- ef
                      ## lowlevel {faster than}: slot(..) <- new("functionX", f)
                      slot(object, fnam)@.Data <- f
                  }
                  object@tDefs <- unlist(dotsargs)
              }

              Erho.val <- integrate(function(x) object@rho(x)*dnorm(x),-Inf, Inf,
                                    rel.tol = .Machine$double.eps^0.5)$value
              Epsi2.val <- integrate(function(x) object@psi(x)^2*dnorm(x),-Inf, Inf,
                                     rel.tol = .Machine$double.eps^0.5)$value
              EDpsi.val <- integrate(function(x) object@Dpsi(x)*dnorm(x),-Inf, Inf,
                                     rel.tol = .Machine$double.eps^0.5)$value
              object@Erho <- new("functionXal", function(arg=1) rep(Erho.val, length(arg)))
              object@Epsi2 <- new("functionXal", function(arg=1) rep(Epsi2.val, length(arg)))
              object@EDpsi <- new("functionXal", function(arg=1) rep(EDpsi.val, length(arg)))
              
              object
          })

## from example(psiFunc)
F0 <- function(x=1, .) rep.int(0, length(x))
F1 <- function(x=1, .) rep.int(1, length(x))
FF1 <- function(.) rep.int(1, length(.))
FF1.2 <- function(.) rep.int(1/2, length(.))
##' Classical psi function
##'
##' Use this psi function to get a classical fit.
##' @export
cPsi <- psiFunc(rho = function(x, .) x^2 / 2, psi = function(x, .) x,
                 wgt = F1, Dwgt = F0, Dpsi = F1, Erho = FF1.2,
                 Epsi2 = FF1, EDpsi = FF1,
                 name = "classic (x^2/2)", . = Inf)

##' Smoothed Huber Function
##'
##' Instead of a sharp bend like the huber function,
##' this function bends smoothly. The first tuning
##' contant, k, can be compared to the tuning constant
##' of the original Huber function. The second tuning
##' constant, s, determines the smoothness of the bend.
##' @title smoothPsi
##' @examples
##'   curve(smoothPsi@@psi(x, 1.345, 10), -3, 3, col="red")
##'   curve(huberPsi@@psi(x, 1.345), -3, 3, add=TRUE)
##' @export
smoothPsi <- psiFuncCached(rho = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, x^2/2, c^2/2 + k*(ax-c) -
                                       ((ax-d)^(1-s) - a^(1-s))/(1-s))
                            },
                            psi = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, x, sign(x)*(k - (ax-d)^(-s)))
                            },
                            Dpsi = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 1, s*(ax-d)^(-s-1))
                            },
                            wgt = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 1, (k - (ax-d)^(-s))/ax)
                            },
                            Dwgt = function(x, k, s) {
                                a <- s^(1/(s+1))
                                c <- k - a^(-s)
                                d <- c - a
                                ax <- abs(x)
                                ifelse(ax <= c, 0,
                                       (ax - d)^(-s-1)*s/x -
                                       (k - (ax-d)^(-s))/(x*ax))
                            },
                            k = 1.345, s = 10,
                            name = "smoothed Huber")
