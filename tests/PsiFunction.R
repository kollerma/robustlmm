require(robustlmm)
require(robustbase)

setClass("psi_func_cached", contains = c("psi_func"))

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

F0 <- function(x=1, .) rep.int(0, length(x))
F1 <- function(x=1, .) rep.int(1, length(x))
cPsi2 <- psiFuncCached(rho = function(x, .) x^2 / 2,
                       psi = function(x, .) x,
                       wgt = F1, Dwgt = F0, Dpsi = F1, 
                       name = "classic (x^2/2)",
                       . = Inf ## dummy, need at least one parameter
)
stopifnot(all.equal(cPsi$Erho(), cPsi2@Erho()),
          all.equal(cPsi$Epsi2(), cPsi2@Epsi2()),
          all.equal(cPsi$EDpsi(), cPsi2@EDpsi()))

smoothPsiOld <- 
  psiFuncCached(rho = function(x, k, s) {
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

chgDefaults.old <- function(object, ...) {
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
}

x <- seq(-10, 10, length.out = 1001)
stopifnot(all.equal(smoothPsi$wgt(x), smoothPsiOld@wgt(x)),
  all.equal(smoothPsi$Dwgt(x), smoothPsiOld@Dwgt(x)),
  all.equal(smoothPsi$psi(x), smoothPsiOld@psi(x)),
  all.equal(smoothPsi$Dpsi(x), smoothPsiOld@Dpsi(x)),
  all.equal(smoothPsi$rho(x), smoothPsiOld@rho(x)),
  all.equal(smoothPsi$Erho(), smoothPsiOld@Erho()),
  all.equal(smoothPsi$Epsi2(), smoothPsiOld@Epsi2()),
  all.equal(smoothPsi$EDpsi(), smoothPsiOld@EDpsi()))

sP <- chgDefaults(smoothPsi, k = 2.0, s = 9.0)
sPOld <- chgDefaults.old(smoothPsiOld, k = 2.0, s = 9.0)
stopifnot(all.equal(sP$wgt(x), sPOld@wgt(x)),
          all.equal(sP$Dwgt(x), sPOld@Dwgt(x)),
          all.equal(sP$psi(x), sPOld@psi(x)),
          all.equal(sP$Dpsi(x), sPOld@Dpsi(x)),
          all.equal(sP$rho(x), sPOld@rho(x)),
          all.equal(sP$Erho(), sPOld@Erho()),
          all.equal(sP$Epsi2(), sPOld@Epsi2()),
          all.equal(sP$EDpsi(), sPOld@EDpsi()))

stopifnot(all.equal(huberPsiRcpp$wgt(x), robustbase::huberPsi@wgt(x)),
          all.equal(huberPsiRcpp$Dwgt(x), robustbase::huberPsi@Dwgt(x)),
          all.equal(huberPsiRcpp$psi(x), robustbase::huberPsi@psi(x)),
          all.equal(huberPsiRcpp$Dpsi(x), robustbase::huberPsi@Dpsi(x) + 0),
          all.equal(huberPsiRcpp$rho(x), robustbase::huberPsi@rho(x)),
          all.equal(huberPsiRcpp$Erho(), robustbase::huberPsi@Erho()),
          all.equal(huberPsiRcpp$Epsi2(), robustbase::huberPsi@Epsi2()),
          all.equal(huberPsiRcpp$EDpsi(), robustbase::huberPsi@EDpsi()))

.psi2propII <- function(object, ...) {
  ## do not do anything for cPsi
  if (identical(object, cPsi)) return(object)
  
  ## Convert a regular psi-function into a proposal II psi function
  ## (with squared weights)
  f <- formals(object@psi)
  nf <- names(f)
  args <- paste(nf, collapse=",")
  x <- nf[1]
  
  ## wgt
  fun <- paste("function(",args,") object@wgt(", args, ")^2")
  wgt <- eval(parse(text=fun))
  formals(wgt) <- f
  ## Dwgt
  fun <- paste("function(",args,") 2*object@wgt(", args, ")*object@Dwgt(",args,")")
  Dwgt <- eval(parse(text=fun))
  formals(Dwgt) <- f
  ## psi
  fun <- paste("function(",args,") object@wgt(", args, ")*object@psi(",args,")")
  psi <- eval(parse(text=fun))
  formals(psi) <- f
  ## Dpsi
  fun <- paste("function(",args,") object@wgt(", args, ")*object@Dpsi(",args,
               ") + object@Dwgt(", args, ")*object@psi(",args,")")
  Dpsi <- eval(parse(text=fun))
  formals(Dpsi) <- f
  ## rho
  intRho <- function(psi, x, ...) {
    ret <- x
    for (i in seq_along(x)) {
      if (is.infinite(x[i])) next
      ret[i] <- integrate(psi, 0, x[i], ..., rel.tol = .Machine$double.eps^0.5)$value
    }
    ret
  }
  fun <- paste("function(",args,") intRho(psi,",args,")")
  rho <- eval(parse(text=fun))
  formals(rho) <- f
  
  ret <- do.call(psiFuncCached, c(list(wgt=wgt, Dwgt=Dwgt, psi=psi, Dpsi=Dpsi, rho=rho),
                                  f[-1], name=paste(object@name, ", Proposal II", sep="")))
  ## if ... is given: pass it to chgDefaults
  chgArgs <- list(...)
  if (length(chgArgs) > 0) {
    if (is.null(names(chgArgs))) stop("Extra arguments in ... need to be named")
    ## extend list, all arguments need to be passed to chgDefaults
    for (name in names(ret@tDefs))
      if (is.null(chgArgs[[name]])) chgArgs[[name]] <- ret@tDefs[[name]]
      ret <- do.call("chgDefaults", c(list(ret), chgArgs))
  }
  return(ret)
}

sP2 <- psi2propII(smoothPsi)
sPOld2 <- .psi2propII(smoothPsiOld)
stopifnot(all.equal(sP2$wgt(x), sPOld2@wgt(x)),
          all.equal(sP2$Dwgt(x), sPOld2@Dwgt(x)),
          all.equal(sP2$psi(x), sPOld2@psi(x)),
          all.equal(sP2$Dpsi(x), sPOld2@Dpsi(x)),
          all.equal(sP2$rho(x), sPOld2@rho(x)),
          all.equal(sP2$Erho(), sPOld2@Erho()),
          all.equal(sP2$Epsi2(), sPOld2@Epsi2()),
          all.equal(sP2$EDpsi(), sPOld2@EDpsi()))

sP <- chgDefaults(smoothPsi, k = 2.0, s = 9.0)
sPOld <- chgDefaults.old(smoothPsiOld, k = 2.0, s = 9.0)
sP2 <- psi2propII(sP)
sPOld2 <- .psi2propII(sPOld)
stopifnot(all.equal(sP2$wgt(x), sPOld2@wgt(x)),
          all.equal(sP2$Dwgt(x), sPOld2@Dwgt(x)),
          all.equal(sP2$psi(x), sPOld2@psi(x)),
          all.equal(sP2$Dpsi(x), sPOld2@Dpsi(x)),
          all.equal(sP2$rho(x), sPOld2@rho(x)),
          all.equal(sP2$Erho(), sPOld2@Erho()),
          all.equal(sP2$Epsi2(), sPOld2@Epsi2()),
          all.equal(sP2$EDpsi(), sPOld2@EDpsi()))


