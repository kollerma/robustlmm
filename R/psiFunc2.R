##' Constructor for psiFunc2
##'
##' Simple constructor for psiFunc2 objects
##' @param name name of psiFunc
##' @param ... passed to psiFunc (psi_func constructor)
##' @export
psiFunc2 <- function(name = NULL, ...) {
    do.call(chgDefaults, c(list(new("psi_func2",
                                    t <- psiFunc(...),
                                    name= name)),
                           t@tDefs))  
}

##' Change the default arguments for a psi_func2 object
##'
##' @title Change default arguments
##' @param ... arguments to change
##' @export
setMethod("chgDefaults", signature("psi_func2"),
          function(object, ...) {
              ##cat("~~~~ chgDefaults of psi_func2 ~~~~~\n")
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

                  for(fnam in list("rho", "psi", "wgt", "Dpsi", "Erho",
                                   "Epsi2", "EDpsi", "lambda")) {
                      f <- slot(object, fnam)
                      ef <- environment(f)
                      if (is(f, "functionXal"))
                          formals(f) <- dotsargs else formals(f)[-1] <- dotsargs
                      environment(f) <- ef
                      ## lowlevel {faster than}: slot(..) <- new("functionX", f)
                      slot(object, fnam)@.Data <- f
                  }
                  object@tDefs <- unlist(dotsargs)
              }
              object
          })


##' Create psi_func2_cached object using cached numerical integration for
##' E... slots.
##'
##' The E... slots will not be fully functional: they just return the
##' value for the current defaults and ignore their arguments.
##'
##' @title psiFuncCached constructor
##' @param rho rho-function
##' @param psi psi-function
##' @param wgt wgt-function
##' @param Dpsi derivative of psi
##' @param name descriptor of this function family
##' @param ... default values for tuning constants
##' @return psi_func2_cached-class object
##' @export
psiFunc2Cached <- function(rho,psi,wgt,Dpsi,name=NULL, ...) {
    lent <- length(dotsargs <- list(...))
    ## '...'  must contain all tuning parameters and their defaults:
    stopifnot(length(nt <- names(dotsargs)) == lent,
              all(nchar(nt)) >= 1)
    if(lent >= 1) {
        ## rho, psi,... checking: must have argument names
        argn <- c("x", nt)
        for(fnam in list("rho", "psi", "wgt", "Dpsi")) {
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
    
    new("psi_func2_cached",
        rho = new("functionX", rho),
        psi = new("functionX", psi),
        wgt = new("functionX", wgt),
        Dpsi= new("functionX", Dpsi),
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
setMethod("chgDefaults", signature("psi_func2_cached"),
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

                  for(fnam in list("rho", "psi", "wgt", "Dpsi")) {
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

.sprintPsiFunc2 <- function(x, short=FALSE) {
    v <- x@tDefs
    n <- names(v)
    name <- x@name
    if (short) name <- gsub('\\s?(psi|function|\\(.*\\))', '', name)
    if (length(v) >= 1) {
        if (short)
            paste(name, paste(n, round(v, 3), sep = "=", collapse = "\n"),
                  sep = "\n")
        else 
            paste(name, " (",
                  paste(n, round(v, 3), sep = " = ", collapse = ", "), ")",
                  sep="")
    } else name
}

##' Print a psi_func2 object
##'
##' @title Print
##' @param x psi_func2 object
##' @param ... ignored
##' @method print psi_func2
##' @S3method print psi_func2
print.psi_func2 <- function(x, ...) print(.sprintPsiFunc2(x))
##' @S3method print psi_func2_cached
print.psi_func2_cached <- function(x, ...) print(.sprintPsiFunc2(x))

##' Huber psi function (psi_func2-class)
##'
##' Same as huberPsi from robustbase package but psi_func2 not psi_func class.
##' @export
huberPsi <- new("psi_func2",
                huberPsi,
                name= "Huber psi function")
                  
## from example(psiFunc)
F1 <- function(x=1) rep.int(1, length(x))
F1.2 <- function(x=1) rep.int(1/2, length(x))
##' Classical psi function
##'
##' Use this psi function to get a classical fit.
##' @export
cPsi <- psiFunc2(rho = function(x) x^2 / 2, psi = function(x) x,
                 wgt = F1, Dpsi = F1, Erho = F1.2,
                 Epsi2 = F1, EDpsi = F1,
                 name = "classic (x^2/2)")

## TODO???: setup conversion from psiFunc to psiFunc2

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
smoothPsi <- psiFunc2Cached(rho = function(x, k, s) {
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
                            k = 1.345, s = 10,
                            name = "smoothed Huber psi function")
