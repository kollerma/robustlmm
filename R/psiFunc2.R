loadModule("psi_function_module", TRUE)

if (FALSE) {
  
  classes <- c("Rcpp_PsiFunction", "Rcpp_HuberPsi", "Rcpp_SmoothPsi",
               "Rcpp_PsiFunctionToPropIIPsiFunctionWrapper")
  
  setMeth <- function(f, def, ...) {
    for (class in classes) {
      cat("setLoadAction(function(ns) setMethod(\"", f, "\", signature(\"", class, "\"), ", 
          substitute(def), "))\n", sep = "")
    }
  }

  
}


## TODO:
## to avoid test is(..., "refClass")
## but this fails:
## setLoadAction(setIs("Rcpp_SmoothPsi", "Rcpp_PsiFunction"))
## setLoadAction(setIs("Rcpp_PsiFunctionToPropIIPsiFunctionWrapper", "Rcpp_PsiFunction"))

##' Change the default arguments for a psi_func_cached object
##'
##' @title Change default arguments
##' @param ... arguments to change
##' @keywords utilities
##' @examples
##' sPsi <- chgDefaults(smoothPsi, k=2)
##' curve(smoothPsi$psi(x), 0, 3)
##' curve(sPsi$psi(x), 0, 3, color="blue", add=TRUE)
##' @exportMethod chgDefaults Rcpp_SmoothPsi
##' @exportMethod chgDefaults Rcpp_PsiFunction
##' @exportMethod chgDefaults Rcpp_PsiFunctionToPropIIPsiFunctionWrapper

.chgDefaults <- function(object, ...) {
  if (identical(object, cPsi))
    return(cPsi)
  clone <- object$copy()
  if (length(list(...)) > 0)
    clone$chgDefaults(c(...))
  return(clone)
}
## setMeth("chgDefaults", .chgDefaults)
setLoadAction(function(ns) setMethod("chgDefaults", signature("Rcpp_PsiFunction"), .chgDefaults))
setLoadAction(function(ns) setMethod("chgDefaults", signature("Rcpp_HuberPsi"), .chgDefaults))
setLoadAction(function(ns) setMethod("chgDefaults", signature("Rcpp_SmoothPsi"), .chgDefaults))

if (FALSE) {
  ## disable this for the time being, once enabled, also enable test in PsiFunction.R
  .chgDefaultsPropII <- function(object, ...) {
    base <- object$base()
    return(psi2propII(base, ...))
  }
  setLoadAction(function(ns) setMethod("chgDefaults", signature("Rcpp_PsiFunctionToPropIIPsiFunctionWrapper"), .chgDefaultsPropII))
}
  
.sprintPsiFunc <- function(x, short=FALSE) {
    v <- x$tDefs()
    n <- names(v)
    ## do not print a single dummy parameter "."
    if (length(n) == 1 && n == ".") {
        v <- numeric(0)
        n <- character(0)
    }
    name <- x$name()
    if (short) name <- gsub('\\s?(psi|function|\\(.*\\))', '', name)
    if (length(v) >= 1) {
        paste(name, " (",
              paste(n, round(v, 3), sep = " = ", collapse = ", "), ")",
              sep="")
    } else name
}

##' \eqn{\psi}{Psi}-functions are used by \code{\link{rlmer}}
##' in the estimating equations and to compute robustness
##' weights. Change tuning parameters using \code{\link{chgDefaults}}
##' and convert to squared robustness weights using the
##' \code{\link{psi2propII}} function.
##' 
##' The \bold{\dQuote{classical} \eqn{\psi}{psi}-function \code{cPsi}}
##' can be used to get a non-robust, i.e., classical, fit.
##' The \code{psi} slot equals the identity function, and
##' the \code{rho} slot equals quadratic function. Accordingly,
##' the robustness weights will always be 1 when using \code{cPsi}.
##'
## The \bold{Huber \eqn{\psi}{psi}-function \code{huberPsi}} is identical to
## the one in the package \code{robustbase}. The \code{psi} slot equals
## the identity function within \eqn{\pm k}{+-k} (where \eqn{k}{k} is
## the tuning parameter). Outside this interval it is equal to
## \eqn{\pm k}{+-k}. The \code{rho} slot equals the quadratic
## function within \eqn{\pm k}{+-k} and a linear function outside.
##
##' The \bold{smoothed Huber \eqn{\psi}{psi}-function} is very similar to
##' the regular Huber \eqn{\psi}{psi}-function.
##' Instead of a sharp bend like the Huber function,
##' the smoothe Huber function bends smoothly. The first tuning
##' contant, k, can be compared to the tuning constant
##' of the original Huber function. The second tuning
##' constant, s, determines the smoothness of the bend.
##'
##' @title Classical, smoothed Huber psi- and rho-functions
##' @name psi-functions
##' @rdname psi-functions
##' @aliases cPsi smoothPsi SmoothPsi PsiFunction
##' @usage ## see examples
##' @seealso \code{\link{chgDefaults}} and \code{\link{psi2propII}}
##' for changing tuning parameters;
##' \code{\link{PsiFunction}} and
##' \code{\link{SmoothPsi}} for a more detailed description of the
##' slots; 
##' @examples
##' plot(cPsi)
##' plot(huberPsiRcpp)
##' plot(smoothPsi)
##' curve(cPsi$psi(x), -3, 3)
##' curve(smoothPsi$psi(x, 1.345, 10), -3, 3, add=TRUE, col="red")
##' curve(huberPsiRcpp$psi(x, 1.345), -3, 3, add=TRUE, col="blue")
##' @export cPsi
setLoadAction(function(ns) assign("cPsi", new(PsiFunction), envir = ns))

##' @export huberPsiRcpp
setLoadAction(function(ns) assign("huberPsiRcpp", new(HuberPsi), envir = ns))

##' @export smoothPsi
setLoadAction(function(ns) assign("smoothPsi", new(SmoothPsi), envir = ns))

##' Converts the psi_func object into a function that corresponds
##' to Proposal II, i.e., a function of the squared weights.
##' The other elements of the psi_func object are adapted accordingly.
##'
##' @title Convert to Propsal II weight function
##' @param object instance of Rcpp_PsiFunction class to convert
##' @param ... optional, new default arguments passed to chgDefaults.
##' @aliases psi2propII,Rcpp_SmoothPsi
##' @keywords utilities
##' @examples
##' par(mfrow=c(2,1))
##' plot(smoothPsi)
##' plot(psi2propII(smoothPsi))
##' @export
setGeneric("psi2propII", function(object, ...) standardGeneric("psi2propII"))
##' @exportMethod psi2propII Rcpp_PsiFunction
##' @exportMethod psi2propII Rcpp_SmoothPsi
##' @exportMethod psi2propII Rcpp_PsiFunctionToPropIIPsiFunctionWrapper

.psi2propII <- function(object, ...) {
  if (identical(object, cPsi))
    return(cPsi)
  clone <- object$copy()
  clone$chgDefaults(object$tDefs())
  if (length(list(...)) > 0)
    clone$chgDefaults(c(...))
  return(new(PsiFunctionToPropIIPsiFunctionWrapper, clone))
}
## setMeth("psi2propII", .psi2propII)
setLoadAction(function(ns) setMethod("psi2propII", signature("Rcpp_PsiFunction"), .psi2propII))
setLoadAction(function(ns) setMethod("psi2propII", signature("Rcpp_HuberPsi"), .psi2propII))
setLoadAction(function(ns) setMethod("psi2propII", signature("Rcpp_SmoothPsi"), .psi2propII))
setLoadAction(function(ns) setMethod("psi2propII", signature("Rcpp_PsiFunctionToPropIIPsiFunctionWrapper"), .psi2propII))

## set show method
##' @exportMethod show Rcpp_PsiFunction
##' @exportMethod show Rcpp_SmoothPsi
##' @exportMethod show Rcpp_PsiFunctionToPropIIPsiFunctionWrapper
.show <- function(object) cat(object$show(), "\n")
## setMeth("show", .show)
setLoadAction(function(ns) setMethod("show", signature("Rcpp_PsiFunction"), .show))
setLoadAction(function(ns) setMethod("show", signature("Rcpp_HuberPsi"), .show))
setLoadAction(function(ns) setMethod("show", signature("Rcpp_SmoothPsi"), .show))
setLoadAction(function(ns) setMethod("show", signature("Rcpp_PsiFunctionToPropIIPsiFunctionWrapper"), .show))

## copied here from robustbase, version 0.92-5
matplotPsi <- function(x, m.psi, psi, par, main = "full",
                       col = c("black", "red3", "blue3", "dark green"),
                       leg.loc = "right", lty = 1, ...) {
  ## Original Author: Martin Maechler, Date: 13 Aug 2010, 10:17
  ## Modified by Manuel Koller, Date: 7 Jan 2013
  fExprs <- quote(list(rho(x), psi(x), {psi*minute}(x),
                       w(x) == psi(x)/x, {w*minute}(x)))
  ## build legend
  map <- if (is.null(colnames(m.psi))) {
    1:(ncol(m.psi)+1)
  } else {
    c(1, c(rho=2, psi=3, Dpsi=4, wgt=5, Dwgt=6)[colnames(m.psi)])
  }
  fExprs <- fExprs[map]
  ## ... title
  if(is.character(main)) {
    shortMain <- (main == "short")
    elist <- list(FF = if(shortMain) fExprs[[2]] else fExprs,
                  PSI = psi, PPP = paste(formatC(par), collapse=","))
    tit <- if(shortMain)
      substitute(FF ~ "etc, with"  ~ psi*"-type" == PSI(PPP), elist)
    else
      substitute(FF ~~ ~~ " with "~~ psi*"-type" == PSI(PPP), elist)
  } else tit <- NULL
  ## plot
  matplot(x, m.psi, col=col, lty=lty, type="l", main = tit,
          ylab = quote(f(x)), xlab = quote(x), ...)
  abline(h=0,v=0, lty=3, col="gray30")
  fE <- fExprs; fE[[1]] <- as.name("expression")
  legend(leg.loc, inset=.02, eval(fE), col=col, lty=lty, bty="n")
  invisible(cbind(x=x, m.psi))
}

plotPsi <- function(x, y, which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
                    main = "full",
                    col = c("black", "red3", "blue3", "dark green", "light green"),
                    leg.loc = "right", ...) {
  ## x: psi_func object
  ## y: points to plot at (x-Axis in plot)
  which <- match.arg(which, several.ok = TRUE)
  if(missing(y)) y <- seq(-5, 10, length=1501)
  tmp <- lapply(which, function(name) 
    eval(substitute(`$`(x, ..name)(y), list(..name = name))))
  m.psi <- do.call(cbind, tmp)
  colnames(m.psi) <- which
  matplotPsi(y, m.psi, x$name(), x$tDefs(),
             main=main, col=col, leg.loc=leg.loc, ...)
}

##' @importFrom robustbase matplotPsi
##' @exportMethod plot
##' @exportMethod plot Rcpp_PsiFunction
##' @exportMethod plot Rcpp_SmoothPsi
##' @exportMethod plot Rcpp_PsiFunctionToPropIIPsiFunctionWrapper
## setMeth("plot", plotPsi)
setLoadAction(function(ns) setMethod("plot", signature("Rcpp_PsiFunction"), plotPsi))
setLoadAction(function(ns) setMethod("plot", signature("Rcpp_HuberPsi"), plotPsi))
setLoadAction(function(ns) setMethod("plot", signature("Rcpp_SmoothPsi"), plotPsi))
setLoadAction(function(ns) setMethod("plot", signature("Rcpp_PsiFunctionToPropIIPsiFunctionWrapper"), plotPsi))