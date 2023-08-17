loadModule("rlmerMatrixUtils_module", TRUE)

##' @importFrom Matrix .bdiag
##' @importFrom methods new setAs setClass setRefClass setMethod
##' @importFrom robustbase summarizeRobWeights
##' @importMethodsFrom robustbase chgDefaults plot
##' @importMethodsFrom Matrix diag solve determinant t crossprod tcrossprod
##'   as.vector drop rowSums rowMeans colSums colMeans chol which


## This is basically a copy of the merPredD-class
##
## The one thing missing is Ptr (and the lock).
##
## REVISION: 1611
##
## @title rlmerPredD
## @name rlmerPredD-class
## @slot all see rlmerPredD class
##' @importMethodsFrom Matrix isDiagonal isTriangular
setRefClass("rlmerPredD",
            fields =
                list(U_e     = "ddiMatrix",
                     v_e     = "numeric",
                     U_b     = "sparseMatrix",
                     .Lambdat = "dgCMatrix",
                     .Lind    = "integer", ## for .Lambdat
                     Lind    = "integer", ## for U_b
                     X       = "matrix",
                     Zt      = "dgCMatrix",
                     beta    = "numeric",
                     b.s     = "numeric",
                     b.r     = "numeric",
                     theta   = "numeric",
                     sigma   = "numeric",
                     n       = "numeric",
                     p       = "numeric",
                     q       = "numeric",
                     lower   = "numeric",
                     zeroB   = "logical",
                     ghz     = "numeric",
                     ghw     = "numeric",
                     ghZ     = "matrix",
                     ghZt    = "matrix",
                     ghW     = "matrix",
                     D_e     = "ddiMatrix",    ## Diagonal Matrix D_e
                     D_b     = "ddiMatrix",    ## Diagonal Matrix D_b
                     Lambda_b = "ddiMatrix",   ## solve(D_b) lambda_e
                     Lambda_bD_b = "ddiMatrix", ## lambda_e I
                     sqrtD_e = "ddiMatrix",    ## sqrt(D_e)
                     U_eX   = "Matrix",       ## U_e %*% X
                     U_eZ   = "Matrix",       ## U_e %*% Z
                     M_XX    = "Matrix",       ## M_XX
                     U_eZU_b = "sparseMatrix",   ## U_e %*% Z %*% U_b
                     calledInit = "logical",   ## whether initMatrices has been called yet
                     M_XZ0   = "Matrix",       ## M_XZ = M_XZ0 %*% U_b
                     M_ZZ0   = "Matrix",       ## M_ZZ = U_b\tr %*% M_ZZ0 %*% U_b
                     M_XX.M_ZZ0 = "Matrix",    ## M_XX^-1 %*% M_XZ0
                     M_ZX0M_XX.M_ZZ0 = "Matrix", ## crossprod(M_XZ0, M_XX^M_ZZ0)
                     rho_e   = "psi_func_rcpp",
                     rho_sigma_e = "psi_func_rcpp",
                     rho_b   = "list",
                     rho_sigma_b = "list",
                     dim     = "numeric",
                     Epsi2_e = "numeric",      ## \Erw[\psi^2_e]
                     Epsi2_b = "numeric",      ## \Erw[\psi^2_b]
                     Epsi_bbt = "Matrix",      ## \Erw[\psi_b(u)u\tr]
                     Epsi_bpsi_bt = "Matrix",  ## \Erw[\psi_b(u)\psi_b(u)\tr]
                     set.M   = "logical",      ## cache
                     cache.M = "list",         ## cache
                     set.unsc = "logical",     ## cache
                     cache.unsc = "matrix",    ## cache
                     MAT1 = "Matrix",          ## used in fitEffects
                     MAT2 = "Matrix"           ## used in fitEffects
                     ),
                methods =
                list(
                     initialize =
                     function(X, Zt, Lambdat, Lind, theta, beta, b.s, lower,
                              sigma, v_e, ...) {
                         set.unsc <<- set.M <<- calledInit <<- FALSE
                         initGH(...)
                         if (!nargs()) return()
                         X <<- as(X, "matrix")
                         Zt <<- as(Zt, "CsparseMatrix")
                         setLambdat(Lambdat, Lind)
                         theta <<- as.numeric(theta)
                         n <<- nrow(X)
                         p <<- ncol(X)
                         q <<- nrow(Zt)
                         stopifnot(length(theta) > 0L,
                                   length(Lind) > 0L,
                                   all(sort(unique(Lind)) == seq_along(theta)),
                                   length(sigma) == 1,
                                   sigma > 0,
                                   length(lower) == length(theta),
                                   length(v_e) == n
                                   )
                         v_e <<- if (is(v_e, "numeric")) v_e else diag(v_e)
                         U_e <<- Diagonal(x=sqrt(if (is(v_e, "numeric")) v_e else diag(v_e)))
                         beta <<- beta
                         b.s <<- b.s
                         b.r <<- as(U_b %*% b.s, "numeric")
                         sigma <<- sigma
                         lower <<- lower
                         setZeroB()
                     },
                     btapply = function(x, fun, ..., rep, add.s=TRUE) {
                         'little helper function to apply a function to all blocktypes'
                         stopifnot(length(x) == length(dim))
                         x <- as.list(x)
                         ret <- list()
                         for (bt in seq_along(dim)) {
                             tmp <- if (add.s) fun(x[[bt]], dim[bt], ...) else fun(x[[bt]], ...)
                             ret <- c(ret, list(tmp))
                         }
                         simplify2array(if (missing(rep)) ret else ret[rep])
                     },
                     setSigma = function(value) {
                         sigma <<- value
                     },
                     setU = function(value) {
                         b.s <<- value
                         b.r <<- as(U_b %*% value, "numeric")
                     },
                     setB.s = function(value) setU(value),
                     setB = function(value) {
                         b.r <<- value
                         b.s <<- stdB(value)
                     },
                     stdB = function(vector) {
                         idx <- diag(U_b) == 0.0
                         if (any(idx)) {
                             lU_b <- U_b
                             diag(lU_b)[idx] <- NA_real_
                             ret <- solve(lU_b, vector)
                             ret[idx] <- 0
                         } else {
                             ret <- solve(U_b, vector)
                         }
                         ret
                     },
                     setLambdat = function(value, Lind) {
                         .Lambdat <<- value
                         .Lind <<- as.integer(Lind)
                         ## get transposed Lind
                         LambdaIdx <- .Lambdat
                         LambdaIdx@x <- as.double(Lind)
                         LambdaIdx <- t(LambdaIdx)
                         Lind <<- as.integer(LambdaIdx@x)
                         ## convert U_b
                         if (isDiagonal(value)) {
                             U_b <<- Diagonal(x=diag(value))
                         } else if (isTriangular(value))  {
                             U_b <<- as(t(value), "triangularMatrix")
                         } else {
                             U_b <<- as(t(value), "triangularMatrix")
                         }
                     },
                     setTheta = function(value) {
                         theta <<- value
                         U_b@x <<- value[Lind]
                         ## resetting factors list
                         if (.hasSlot(U_b, "factors"))
                             U_b@factors <<- list()
                         ## update zeroB
                         setZeroB()
                         ## update U_eZU_b
                         U_eZU_b <<- U_eZ %*% U_b
                         ## clear cache
                         set.unsc <<- set.M <<- FALSE
                     },
                     Lambdat = function() {
                         .Lambdat@x <<- theta[.Lind]
                         .Lambdat@factors <<- list()
                         .Lambdat
                     },
                     setZeroB = function() {
                         zeroB <<- if (any(theta[lower >= 0] == 0)) {
                             diag(U_b) == 0 ## this works because we only allow all-non-zero blocks
                         } else {
                             rep(FALSE, q)
                         }
                     },
                     initMatrices = function(object) {
                         if (isTRUE(calledInit)) return()
                         ## initialize often required matrices and other stuff
                         dim <<- object@dim
                         U_eX <<- as(U_e %*% X, "denseMatrix")
                         U_eZ <<- tcrossprod(U_e, Zt)
                         U_eZU_b <<- U_eZ %*% U_b
                         MAT1 <<- Matrix(0, n+q, p+q)
                         MAT1[1:n,1:p] <<- U_eX
                         MAT2 <<- Matrix(0, n, p+q)
                         MAT2[,1:p] <<- U_eX
                         initRho(object)
                         calledInit <<- TRUE
                     },
                     M = function() {
                         if (set.M) return(cache.M)
                         cache.M <<- list()
                         if (any(!zeroB)) {
                             ## M_bb. := M_bb\inv
                             cache.M$M_bb. <<- as(crossprod(U_b,(M_ZZ0 - M_ZX0M_XX.M_ZZ0)) %*% U_b +
                                 Lambda_bD_b, "denseMatrix")
                             cache.M$M_XZ <<- M_XZ <- as(M_XZ0 %*% U_b, "denseMatrix")
                             M_ZX.M_XX <- t(solve(M_XX, M_XZ))
                             cache.M$M_bB <<- as(-1*solve(cache.M$M_bb., M_ZX.M_XX), "denseMatrix")
                         } else { ## all random effects dropped
                             cache.M$M_bb. <<- as(Lambda_bD_b, "denseMatrix")
                             cache.M$M_XZ <<- M_XZ <- as(Matrix(0, p, q), "denseMatrix")
                             cache.M$M_bB <<- matrix(0, q, p)
                         }
                         cache.M$M_BB <<- as(solve(M_XX, Diagonal(p) - M_XZ %*% cache.M$M_bB), "denseMatrix")
                         cache.M$M_bb <<- as(solve(cache.M$M_bb.), "denseMatrix")
                         set.M <<- TRUE
                         return(cache.M)
                     },
                     initRho = function(object) {
                         rho_e <<- object@rho.e
                         rho_sigma_e <<- object@rho.sigma.e
                         rho_b <<- object@rho.b
                         rho_sigma_b <<- object@rho.sigma.b
                         D_e <<- Diagonal(x=rep.int(rho_e@EDpsi(), n))
                         tmp <- btapply(rho_b, .calcE.D.re2, rep=object@ind[object@k])
                         D_b <<- Diagonal(x=tmp)
                         Lambda_b <<- Diagonal(x=rho_e@EDpsi() / tmp)
                         Lambda_bD_b <<- Diagonal(x=rep(rho_e@EDpsi(), q))
                         sqrtD_e <<- Diagonal(x=sqrt(D_e@x))
                         M_XX <<- crossprod(sqrtD_e %*% U_eX)
                         tmp1 <- sqrtD_e %*% U_eZ
                         M_XZ0 <<- crossprod(sqrtD_e %*% U_eX, tmp1)
                         M_ZZ0 <<- crossprod(tmp1)
                         M_XX.M_ZZ0 <<- solve(M_XX, M_XZ0)
                         M_ZX0M_XX.M_ZZ0 <<- crossprod(M_XZ0, M_XX.M_ZZ0)
                         Epsi2_e <<- rho_e@Epsi2()
                         ## calculate Epsi_bbt
                         Epsi_bbt <<- .bdiag(btapply(rho_b, .calcE.psi_bbt, rep=object@ind))
                         ## calculate Epsi_bpsi_bt
                         Epsi_bpsi_bt <<- .bdiag(btapply(rho_b, .calcE.psi_bpsi_bt, rep=object@ind))
                         Epsi2_b <<- diag(Epsi_bpsi_bt) ## for easier computation in diagonal case
                     },
                     initGH = function(numpoints=13) {
                         ### ::: required to make roxygen2 work
                         gh <- robustlmm:::ghq(numpoints, FALSE)
                         ghz <<- gh$nodes * sqrt(2)
                         ghw <<- gh$weights / sqrt(pi)
                         ghZ <<- matrix(ghz, numpoints, numpoints)
                         ghZt <<- t(ghZ)
                         ghW <<- ghw %o% ghw
                     },
                     updateMatrices = function() {
                         ## does not do anything
                         invisible(TRUE)
                     },
                     unsc = function() {
                         'the unscaled variance-covariance matrix of the fixed-effects parameters'
                         if (set.unsc) return(cache.unsc)
                         r <- M()
                         tmp <- if (all(zeroB)) { ## all theta == 0
                             Epsi2_e / rho_e@EDpsi() *
                                 tcrossprod(solve(M_XX, t(sqrtD_e %*% U_eX)))
                         } else {
                             Epsi2_e / rho_e@EDpsi() *
                                 with(r, M_BB - crossprod(M_bB, Lambda_bD_b %*% M_bB)) +
                                 if (isDiagonal(U_b)) {
                                     crossprod(Diagonal(x=sqrt(Epsi2_b)) %*% Lambda_b %*% r$M_bB)
                                 } else {
                                     with(r, crossprod(Lambda_b %*% M_bB,
                                                       Epsi_bpsi_bt %*% Lambda_b %*% M_bB))
                                 }
                         }
                         cache.unsc <<- as.matrix(tmp)
                         ## test matrix for symmetry
                         ## if its almost symmetric, then return symmetrized matrix
                         ## otherwise summary.merMod() and chol() will complain.
                         if (!isSymmetric(cache.unsc, 0)) {
                             cache.unsc <<- symmpart(cache.unsc)
                         }
                         set.unsc <<- TRUE
                         return(cache.unsc)
                     },
                     ## functions for compatibility with lme4
                     b = function(fac) {
                         stopifnot(isTRUE(all.equal(fac, 1.)))
                         b.r
                     }
                     )
                )

setRefClass("rlmerPredD_DAS",
            fields =
            list(diagA   = "numeric",       ## Digonal of matrix A
                 diagAAt = "numeric",       ## Diagonal of matrix A %*% t(A)
                 groupsA = "numeric",       ## group ids for groups of equal values in diagA
                 Kt      = "Matrix",        ## Matrix B = lfrac * Kt, K = t(Kt)
                 L       = "Matrix",        ## Matrix L
                 kappa_e = "numeric",       ## kappa_e^(sigma) (only required for DASvar and DAStau)
                 kappa_b = "numeric",       ## kappa_b^(sigma) (-------------- " ------------------)
                 EDpsi_e = "numeric",
                 method  = "character",
                 .setTau_e = "logical",     ## check whether we already have computed tau_e
                 .tau_e  = "numeric",       ## tau_e used in updateSigma
                 .setTbk = "logical",
                 .Tbk    = "list",        ## cache for T_{b,k}
                 nblocks = "numeric",
                 blocks  = "list",
                 idx     = "list"
                 ),
            contains = "rlmerPredD",
            methods =
            list(
                 initialize = function(...) {
                     callSuper(...)
                     .setTau_e <<- FALSE
                     .tau_e <<- numeric(0)
                     .setTbk <<- FALSE
                     .Tbk <<- list()
                     ## the other slots are initialized in initMatrices()
                 },
                 initMatrices = function(object) {
                    callSuper(object)
                    blocks <<- object@blocks
                    idx <<- object@idx
                    nblocks <<- sum(sapply(idx, ncol))
                    groupsA <<- 0:(n-1)
                 },
                 initRho = function(object) {
                     callSuper(object)
                     EDpsi_e <<- rho_e@EDpsi()
                     kappa_e <<- calcKappaTau(rho_sigma_e, 1)
                     kappa_b <<- calcKappaTauB(object)
                 },
                 setTheta = function(value) {
                     callSuper(value)
                     ## update Matrices
                     updateMatrices()
                     ## need to recompute tau and co.
                     .setTau_e <<- FALSE
                 },
                 B = function() {
                     Kt %*% Lambda_b
                 },
                 K = function() {
                     t(Kt)
                 },
                 TbList = function() {
                     tmp <- as.matrix(L %*% Epsi_bbt)
                     Tfull <- diag(q) - tmp - t(tmp) +
                         as.matrix(Epsi2_e * crossprod(Kt) +
                                       L %*% crossprod(Epsi_bpsi_bt, L))
                     Ts <- vector("list", nblocks)
                     i <- 1
                     ## cycle blocks
                     for (type in seq_along(blocks)) {
                         bidx <- idx[[type]]
                         for (k in 1:ncol(bidx)) { ## 1:K
                             Tblock <- symmpart(Tfull[bidx[,k], bidx[,k],
                                                      drop=FALSE])
                             Ts[[i]] <- Tblock
                             i <- i + 1
                         }
                     }
                     Ts
                 },
                 Tb = function() {
                     .bdiag(TbList())
                 },
                 setTList = function(Tbk) {
                     .Tbk <<- Tbk
                     .setTbk <<- TRUE
                 },
                 TList = function() {
                     if (.setTbk) .Tbk else TbList() ## fallback to TbList()
                 },
                 A = function() {
                     if (!isTRUE(calledInit)) initMatrices()
                     if (any(!zeroB)) {
                         r <- M()
                         A <- tcrossprod(U_eX, ## X M_Bb U_b Z U_e
                                         U_eZU_b %*% r$M_bB)
                         A <- A + t(A) + tcrossprod(U_eX, U_eX %*% r$M_BB) +
                             tcrossprod(U_eZU_b %*% r$M_bb, U_eZU_b)
                     } else {
                         ## no random effects
                         A <- U_eX %*% solve(M_XX, t(U_eX)) ## just the hat matrix
                     }
                     return(A)
                 },
                 updateMatrices = function() {
                     if (!isTRUE(calledInit)) initMatrices()
                     if (any(!zeroB)) {
                         r <- M()
                         tmp1 <- U_eZU_b %*% r$M_bb ## U_e Z U_b M_bb
                         result <- calculateA(U_eX, U_eZU_b, tmp1, r$M_bB,
                                              r$M_BB, as.integer(groupsA))
                         diagA <<- result[["diagA"]]
                         diagAAt <<- result[["diagAAt"]]
                         groupsA <<- result[["groupsA"]]
                         Kt <<- -1*(tcrossprod(U_eX, r$M_bB) + tmp1)
                         L <<- r$M_bb %*% Lambda_b
                     } else {
                         ## no random effects
                         ## FIXME also do this in c++?
                         tmp <- solve(M_XX, t(U_eX))
                         diagAAt <<- diagA <<- numeric(n)
                         for (i in 1:n) {
                             Arow <- U_eX[i, ] %*% tmp
                             diagA[i] <<- Arow[i]
                             diagAAt[i] <<- sum(Arow * Arow)
                         }
                         Kt <<- as(Matrix(0, n, q), "unpackedMatrix")
                         L <<- solve(D_b)
                     }
                 },
                 tau_e = function() {
                     if (isTRUE(.setTau_e)) return(.tau_e)
                     Btmp <- B()
                     if (method == "DASvar" || length(.tau_e) == 0) {
                         tmp <- tcrossprod(Epsi_bpsi_bt, Btmp)
                         tau2 <-
                             v_e - EDpsi_e * 2 * diagA + Epsi2_e * diagAAt +
                             computeDiagonalOfProduct(as(Btmp, "unpackedMatrix"),
                                                      as(tmp, "unpackedMatrix"))
                         tooSmall <- tau2 < 0.01
                         if (any(tooSmall)) {
                            tau2[tooSmall] <- 0.01
                            obsString <- createObservationsString(tooSmall)
                            warning("Detected very small values for tau^2 for ",
                                    obsString, ". Using 0.01 instead.")
                         }
                         .tau_e <<- sqrt(tau2)
                         if (any(is.na(.tau_e))) browser()
                     }
                     if (method == "DAStau") {
                         stmp <- .s(theta = FALSE, pp = .self, B = Btmp)
                         .tau_e <<-
                             calcTau(diagA, stmp, rho_e, rho_sigma_e, .self, kappa_e, .tau_e)
                     }
                     .setTau_e <<- TRUE
                     return(.tau_e)
                 }
                )
            )

setRefClass("rlmerResp",
            fields =
                list(mu      = "numeric",
                     offset  = "numeric",
                     sqrtrwt = "numeric",
                     weights = "numeric",
                     wtres   = "numeric",
                     y       = "numeric"
                ),
            methods =
                list(
                    initialize = function(...) {
                        if (!nargs()) return()
                        ll <- list(...)
                        if (is.null(ll$y)) stop("y must be specified")
                        y <<- as.numeric(ll$y)
                        n <- length(y)
                        mu <<- if (!is.null(ll$mu))
                            as.numeric(ll$mu) else numeric(n)
                        offset <<- if (!is.null(ll$offset))
                            as.numeric(ll$offset) else numeric(n)
                        weights <<- if (!is.null(ll$weights))
                            as.numeric(ll$weights) else rep.int(1,n)
                        sqrtrwt <<- if (!is.null(ll$sqrtrwt))
                            as.numeric(ll$sqrtrwt) else sqrt(weights)
                        wtres   <<- sqrtrwt * (y - mu)
                    },
                    updateMu = function(lmu) {
                        mu <<- lmu
                        wtres <<- sqrtrwt * (y - mu)
                    })
)

##' Class "rlmerMod" of Robustly Fitted Mixed-Effect Models
##'
##' A robust mixed-effects model as returned by \code{\link{rlmer}}.
##' @title rlmerMod Class
##' @name rlmerMod-class
##' @aliases rlmerMod-class coef.rlmerMod deviance.rlmerMod extractAIC.rlmerMod
##'   family.rlmerMod fitted.rlmerMod fixef.rlmerMod formula.rlmerMod
##'   isGLMM.rlmerMod isLMM.rlmerMod isNLMM.rlmerMod isREML.rlmerMod
##'   logLik.rlmerMod model.frame.rlmerMod model.matrix.rlmerMod nobs.rlmerMod
##'   predict.rlmerMod print.rlmerMod print.summary.rlmer print.VarCorr.rlmerMod
##'   ranef.rlmerMod resid.rlmerMod rlmerMod-class sigma.rlmerMod show.rlmerMod
##'   show,rlmerMod-method show.summary.rlmerMod summary.rlmerMod
##'   summary.summary.rlmerMod terms.rlmerMod update.rlmerMod VarCorr.rlmerMod
##'   VarCorr.summary.rlmerMod vcov.rlmerMod vcov.summary.rlmerMod
##'   weights.rlmerMod
##' @docType class
##' @section Objects from the Class: Objects are created by calls to
##'   \code{\link{rlmer}}.
##' @section Methods: Almost all methods available from objects returned from
##'   \code{\link{lmer}} are also available for objects returned by
##'   \code{\link{rlmer}}. They usage is the same.
##'
##'   It follows a list of some the methods that are exported by this package:
##'
##' \itemize{
##' \item \code{\link{coef}}
##' \item \code{\link{deviance}} (disabled, see below)
##' \item \code{\link{extractAIC}} (disabled, see below)
##' \item \code{\link{family}}
##' \item \code{\link{fitted}}
##' \item \code{\link[=fixef.merMod]{fixef}}
##' \item \code{\link{formula}}
##' \item \code{\link{getInfo}}
##' \item \code{\link{isGLMM}}
##' \item \code{\link{isLMM}}
##' \item \code{\link{isNLMM}}
##' \item \code{\link{isREML}}
##' \item \code{\link{logLik}} (disabled, see below)
##' \item \code{\link{model.frame}}
##' \item \code{\link{model.matrix}}
##' \item \code{\link{nobs}}
##' \item \code{\link[=plot.rlmerMod]{plot}}
##' \item \code{\link[=predict.merMod]{predict}}
##' \item \code{\link[=ranef.merMod]{ranef}} (only partially implemented)
##' \item \code{\link[=residuals.rlmerMod]{residuals}}
##' \item \code{\link{sigma}}
##' \item \code{\link{summary}}
##' \item \code{\link{terms}}
##' \item \code{\link{update}}
##' \item \code{\link[=VarCorr.merMod]{VarCorr}}
##' \item \code{\link{vcov}}
##' \item \code{\link{weights}}
##' }
##' @section Disabled methods: A log likelihood or even a pseudo log likelihood
##'   is not defined for the robust estimates returned by \code{\link{rlmer}}.
##'   Methods that depend on the log likelihood are therefore not available. For
##'   this reason the methods \code{deviance}, \code{extractAIC} and
##'   \code{logLik} stop with an error if they are called.
##' @seealso \code{\link{rlmer}}; corresponding class in package \code{lme4}:
##'   \code{\link{merMod}}
##' @keywords classes
##' @examples
##'
##' showClass("rlmerMod")
##'
##' ## convert an object of type 'lmerMod' to 'rlmerMod'
##' ## to use the methods provided by robustlmm
##' fm <- lmer(Yield ~ (1|Batch), Dyestuff)
##' rfm <- as(fm, "rlmerMod")
##' compare(fm, rfm)
##'
##' @export
setClass("rlmerMod",
         representation(resp    = "rlmerResp",
                        Gp      = "integer",
                        call    = "call",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        flist   = "list",
                        cnms    = "list",
                        lower   = "numeric",
                        theta   = "numeric",
                        beta    = "numeric",
                        b.s     = "numeric",
                        devcomp = "list",
                        pp      = "refClass",
                        optinfo = "list",
                        ## from rlmerResp:
                        rho.e   = "psi_func_rcpp",
                        rho.sigma.e = "psi_func_rcpp",
                        method  = "character",
                        ## from rreTrms:
                        b.r     = "numeric",
                        rho.b   = "list",
                        rho.sigma.b = "list",
                        blocks  = "list",
                        ind     = "numeric",
                        idx     = "list",
                        dim     = "numeric",
                        q       = "numeric",
                        k       = "numeric"
                        ),
         validity = function(object) {
             v <- validObject(object@rho.e)
             if (!is.logical(v) || ! v)
                 return(v)
             v <- validObject(object@rho.sigma.e)
             if (!is.logical(v) || ! v)
                 return(v)
             if (length(object@b.r) != length(object@b.s))
                 return("b and u have to be of the same length")
             if (length(object@blocks) != max(object@ind))
                 return("number of blocks and maximum index in ind must coincide")
             if (length(object@blocks) != length(object@dim))
                 return("number of blocks and length of dim must coincide")
             if (! identical(object@dim, sapply(object@idx, NROW)))
                 return("number of rows in idx not equal to dim")
             nblocks <- sapply(object@idx, NCOL)
             if (! identical(object@q, nblocks * object@dim))
                 return("number of columns in idx not equal to number of subjects")
             if (sum(object@q) != length(object@b.s))
                 return("sum of q must be identical to length of u")
             if (length(object@k) != length(object@b.s))
                 return("length of k must be identical to length of u")
             if (length(object@ind) != max(object@k))
                 return("length of ind must be equal to max k")
             if (max(object@k) != sum(nblocks) || any(object@k <= 0))
                 return("block index in k must be positive and not exceed the number of blocks")
             if (FALSE && any(unlist(lapply(object@blocks, function(x) unique(as.vector(x)))) !=
                     unique(unlist(lapply(object@blocks, as.vector)))))
                 return("parameters are not allowed to be in multiple blocks simultaneously")
             if (length(object@rho.b) != length(object@dim))
                 return("the length of list rho.b must be the same as the number ob block types")
             if (length(object@rho.b) != length(object@rho.sigma.b))
                 return("the lists rho.b and rho.sigma.b must be of the same length")
             if (!all(sapply(object@rho.b, inherits, "psi_func_rcpp")))
                 return("the entries of rho.b must be of class psi_func_rcpp or inherited")
             if (!all(sapply(object@rho.sigma.b, inherits, "psi_func_rcpp")))
                 return("the entries of rho.sigma.b must be of class psi_func_rcpp or inherited")
             TRUE
         },
         S3methods = TRUE)

#######################################################
## Coercing methods                                  ##
#######################################################

.convLme4Rlmer <- function(from) {
    X <- getME(from, "X")
    Zt <- getME(from, "Zt")
    if (is(from, "merMod")) {
        Lambdat <- getME(from, "Lambdat")
        Lind <- getME(from, "Lind")
        u <- getME(from, "u")
        lower <- getME(from, "lower")
        devcomp <- getME(from, "devcomp")
        theta <- getME(from, "theta")
        mu <- getME(from, "mu")
        cnms <- from@cnms
        y <- getME(from, "y")
        offset <- getME(from, "offset")
        weights <- weights(from)
    } else stop("Unsupported object of class", class(from))

    resp <- new("rlmerResp",
                mu = mu,
                offset = offset,
                weights = weights,
                y = y)
    pp <- new("rlmerPredD",
              X=X,
              Zt=Zt,
              Lambdat=Lambdat,
              Lind=Lind,
              theta=theta,
              beta=getME(from, "beta"),
              b.s=u,
              lower=lower,
              v_e=resp$weights,
              sigma=sigma(from),
              numpoints=13)
    b <- findBlocks(Lambdat=Lambdat, Lind=Lind)
    to <- new("rlmerMod",
              resp=resp,
              Gp=getME(from, "Gp"),
              call=from@call,
              frame=from@frame,
              flist=getME(from, "flist"),
              cnms=cnms,
              lower=lower,
              theta=theta,
              beta=getME(from, "beta"),
              b.s=u,
              devcomp=devcomp,
              pp=pp,
              optinfo=list(),
              rho.e=cPsi,
              rho.sigma.e=cPsi,
              method="DAS",
              b.r=pp$b.r,
              rho.b=rep.int(list(cPsi),length(b$dim)),
              rho.sigma.b=rep.int(list(cPsi),length(b$dim)),
              blocks=b$blocks,
              ind=b$ind,
              idx=b$idx,
              dim=b$dim,
              q=b$q,
              k=b$k)
    to@pp$initMatrices(to)
    to
}

setAs("lmerMod", "rlmerMod", .convLme4Rlmer)

setAs("rlmerPredD", "rlmerPredD_DAS", function(from) {
    to <- new("rlmerPredD_DAS")
    fields <- names(getRefClass(class(from))$fields())
    Map(function(field) to$field(field, from$field(field)), fields)
    to
})

#######################################################
## Update methods                                    ##
#######################################################

## Update robustness weights after changing any of the relevant
## parts in a rlmerMod object
##
## @title Update robustness weights
## @param object rlmerMod object to update
## @return rlmerMod object
updateWeights <- function(object) {
    ## copy estimates from pp to object
    object@theta <- theta(object)
    object@beta <- .fixef(object)
    object@b.s <- b.s(object)
    object@b.r <- .b(object)
    dd <- object@devcomp$dims
    object@devcomp$cmp[ifelse(dd["REML"], "sigmaREML", "sigmaML")] <- .sigma(object)
    ## Set slots to NA
    object@devcomp$cmp[c("ldL2", "ldRX2", "wrss", "ussq", "pwrss", "drsum", "REML", "dev")] <- NA

    object
}

## update pp slot from values in object
##
## (this does the reverse of updateWeights
##
## @title set values in pp
## @param object rlmerMod object to update
## @return nothing
updatePp <- function(object) {
    object@pp$lower <- object@lower
    dd <- object@devcomp$dims
    object@pp$setSigma(object@devcomp$cmp[[ifelse(dd["REML"], "sigmaREML", "sigmaML")]])
    object@pp$setB.s(object@b.s)
    invisible(object)
}
