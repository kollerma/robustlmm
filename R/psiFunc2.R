loadModule("psi_function_module", TRUE)

#' @importFrom robustbase psiFunc
#' @importClassesFrom robustbase psi_func
setClass("psi_func_rcpp",
         slots = c(getRcppClass = "function",
                   getInstanceWithOriginalDefaults = "function"),
         contains = "psi_func")

fixTDefs <- function(..., defaultTDefs) {
  tDefs <- c(...)
  if (length(tDefs) > 0) {
    if (is.null(names(tDefs))) {
      if (length(tDefs) > length(defaultTDefs)) {
        stop("Expected only ", length(defaultTDefs), " arguments but got ",
             length(tDefs), ".")
      }
      names(tDefs) <- names(defaultTDefs)[seq_along(tDefs)]
    } else {
      if (any(sapply(names(tDefs), is.null))) {
        stop("Either all parameters need to be named or none of them.")
      }
    }
    unknownTDefs <- setdiff(names(tDefs), names(defaultTDefs))
    if (length(unknownTDefs) > 0) {
      stop("Found tuning parameter names not present in defaults: ",
           paste(unknownTDefs, collapse = ", "), ".")
    }
  }
  missingTDefs <- setdiff(names(defaultTDefs), names(tDefs))
  if (length(missingTDefs) > 0) {
    tDefs <- c(tDefs, defaultTDefs[missingTDefs])
  }
  if (length(tDefs) == 0) {
    tDefs <- c(`.` = NA_real_)
  } else {
    tDefs <- tDefs[names(defaultTDefs)]
  }
  return(tDefs)
}

psiFuncRcpp <- function(rcppClass, ...) {
  ..RcppClass.. <- rcppClass
  createInstance <- function() {
    if (length(..RcppClass..) == 1) {
      return(new(get(..RcppClass..)))
    } else if (length(..RcppClass..) == 2) {
      innerInstance <- new(get(..RcppClass..[2]))
      return(new(get(..RcppClass..[1]), innerInstance))
    } else {
      stop("Can only have rcppClass of length 1 or 2.")
    }
  }
  ..Instance.. <- createInstance()
  ..TDefs.. <- fixTDefs(..., defaultTDefs = ..Instance..$tDefs())
  ..Instance..$chgDefaults(..TDefs..)
  argNames <- paste(names(..TDefs..), collapse = ", ")
  getPointer <- function() {
    return(get(".pointer", envir = as.environment(..Instance..)))
  }
  getInstance <- function(args) {
    if (isnull(getPointer())) {
      ..Instance.. <<- createInstance()
    }
    ..Instance..$chgDefaults(args)
    return(..Instance..)
  }
  createFun <- function(method, envir = parent.frame()) {
    text <- paste0("function(x, ", argNames, ") {\n",
                   "getInstance(c(", argNames, "))$", method, "(x) }")
    fun <- eval(parse(text = text), envir = envir)
    formals(fun)[-1] <- ..TDefs..
    return(fun)
  }
  createEFun <- function(method, envir = parent.frame()) {
    if (length(..TDefs..) == 1) {
      text <- paste0("function(", argNames, ") {\n",
                     "if (missing(", argNames, ")) {\n",
                     "return(getInstance(", argNames, ")$", method, "())\n",
                     "} else {\n",
                     "return(sapply(", argNames, ", function(arg)\n",
                     "getInstance(arg)$", method, "())) } }")
    } else {
      text <- paste0("function(", argNames, ") {\n",
                     "getInstance(c(", argNames, "))$", method, "() }")
    }
    fun <- eval(parse(text = text), envir = envir)
    formals(fun) <- ..TDefs..
    return(fun)
  }
  rho <- createFun("rho")
  psi <- createFun("psi")
  wgt <- createFun("wgt")
  Dpsi <- createFun("Dpsi")
  Dwgt <- createFun("Dwgt")
  Erho <- createEFun("Erho")
  Epsi2 <- createEFun("Epsi2")
  EDpsi <- createEFun("EDpsi")
  func <- do.call(psiFunc,
                  c(list(rho, psi, wgt, Dpsi, Dwgt, Erho, Epsi2, EDpsi,
                         ..Instance..$name()), ..TDefs..))
  args <- list("psi_func_rcpp")
  for(slotName in slotNames("psi_func")) {
     args[[slotName]] <- slot(func, slotName)
  }
  args[["getRcppClass"]] <- function() {
    return(..RcppClass..)
  }
  args[["getInstanceWithOriginalDefaults"]] <- function() {
    instance <- createInstance()
    instance$chgDefaults(..TDefs..)
    return(instance)
  }
  funcRcpp <- do.call("new", args)
  return(funcRcpp)
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
##' the smoothed Huber function bends smoothly. The first tuning
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
##' curve(cPsi@psi(x), -3, 3)
##' curve(smoothPsi@psi(x, 1.345, 10), -3, 3, add=TRUE, col="red")
##' curve(huberPsiRcpp@psi(x, 1.345), -3, 3, add=TRUE, col="blue")
##' @export cPsi
setLoadAction(function(ns) assign("cPsi", psiFuncRcpp("PsiFunction"), envir = ns))

##' @export huberPsiRcpp
setLoadAction(function(ns) assign("huberPsiRcpp", psiFuncRcpp("HuberPsi"), envir = ns))

##' @export smoothPsi
setLoadAction(function(ns) assign("smoothPsi", psiFuncRcpp("SmoothPsi"), envir = ns))

.chgDefaults <- function(object, ...) {
  if (identical(object, cPsi))
    return(cPsi)
  clone <- do.call(psiFuncRcpp,
                   c(list(object@getRcppClass()),
                     fixTDefs(..., defaultTDefs = object@tDefs)))
  return(clone)
}
##' Change the default arguments for a psi_func_rcpp object
##'
##' @title Change default arguments
##' @param ... arguments to change
##' @keywords utilities
##' @examples
##' sPsi <- chgDefaults(smoothPsi, k=2)
##' curve(smoothPsi@@psi(x), 0, 3)
##' curve(smoothPsi@@psi(x), 0, 3, color="blue", add=TRUE)
##' @export
setMethod("chgDefaults", signature("psi_func_rcpp"), .chgDefaults)

.sprintPsiFunc <- function(x, short=FALSE) {
  v <- x@tDefs
  n <- names(v)
  ## do not print a single dummy parameter "."
  if (length(n) == 1 && n == ".") {
    v <- numeric(0)
    n <- character(0)
  }
  name <- x@name
  if (short) name <- gsub('\\s?(psi|function|\\(.*\\))', '', name)
  if (length(v) >= 1) {
    paste(name, " (",
          paste(n, round(v, 3), sep = " = ", collapse = ", "), ")",
          sep="")
  } else name
}

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
##' @exportMethod psi2propII psi_func_rcpp

.psi2propII <- function(object, ...) {
  if (identical(object, cPsi))
    return(cPsi)
  if (object@getRcppClass()[1] == "PsiFunctionToPropIIPsiFunctionWrapper") {
    stop("Cannot apply psi2propII multiple times.")
  }
  func <- do.call(psiFuncRcpp,
                  c(list(c("PsiFunctionToPropIIPsiFunctionWrapper",
                           object@getRcppClass()),
                         fixTDefs(..., defaultTDefs = object@tDefs))))
  return(func)
}

##' @rdname psi2propII
setMethod("psi2propII", signature("psi_func_rcpp"), .psi2propII)
