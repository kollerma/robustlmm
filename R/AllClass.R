loadModule("rlmerPredD_module", TRUE)
loadModule("rlmerResp_module", TRUE)
loadModule("fitEffects_module", TRUE)

##' @importFrom Matrix bdiag
##' @importFrom methods new setAs setClass setRefClass setMethod
##' @importFrom robustbase summarizeRobWeights
##' @importMethodsFrom robustbase chgDefaults plot
##' @importMethodsFrom Matrix diag solve determinant t crossprod tcrossprod as.vector drop rowSums rowMeans colSums colMeans chol which


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
                     V_e     = "ddiMatrix",        ## crossprod(U_e)
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
                     .U_eX   = "Matrix",       ## U_e\inv X
                     .U_eZ   = "Matrix",       ## U_e\inv Z
                     M_XX    = "Matrix",       ## M_XX
                     U_btZt.U_et = "Matrix",   ## U_b\tr Z\tr U_e\inv\tr
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
                     cache.unsc = "matrix"        ## cache
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
                         Zt <<- as(Zt, "dgCMatrix")
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
                         U_e <<- Diagonal(x=sqrt(if (is(v_e, "numeric")) v_e else diag(v_e)))
                         V_e <<- Diagonal(x=if (is(v_e, "numeric")) v_e else diag(v_e))
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
                     setU = function(value) setB.s(value),
                     setB.s = function(value) {
                         b.s <<- value
                         b.r <<- as(U_b %*% value, "numeric")
                     },
                     setB = function(value) {
                         b.r <<- value
                         b.s <<- as(stdB(1, Matrix(value)), "numeric")
                     },
                     stdB = function(sigma = sigma, matrix, drop=TRUE, t=FALSE) {
                         ## matrix not given: return sigma-standardized b.s
                         if (missing(matrix)) return(b.s / sigma)
                         ## else standardize matrix
                         ret <- solve(if (t) t(U_b) else U_b, matrix)
                         ## set infinite, NA, NaN terms to 0 (they correspond to dropped vc)
                         if (any(idx <- !is.finite(ret@x)))
                             ret@x[idx] <- 0
                         if (drop) ret <- drop(ret)
                         ret/sigma
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
                             U_b <<- as(t(value), "dtCMatrix")
                         } else {
                             U_b <<- as(t(value), "dgCMatrix")
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
                         ## update U_btZt.U_et
                         U_btZt.U_et <<- t(.U_eZ %*% U_b)
                         ## clear cache
                         ##cat("set theta =", theta, "\n")
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
                         .U_eX <<- solve(U_e, X)
                         .U_eZ <<- solve(U_e, t(Zt))
                         U_btZt.U_et <<- t(.U_eZ %*% U_b)
                         initRho(object)
                         calledInit <<- TRUE
                     },
                     M = function() {
                         if (set.M) return(cache.M)
                         cache.M <<- list()
                         if (any(!zeroB)) {
                             ## M_bb. := M_bb\inv
                             cache.M$M_bb. <<- crossprod(U_b,(M_ZZ0 - M_ZX0M_XX.M_ZZ0)) %*% U_b +
                                 Lambda_bD_b
                             cache.M$M_XZ <<- M_XZ <- M_XZ0 %*% U_b
                             M_ZX.M_XX <- t(solve(M_XX, M_XZ))
                             cache.M$M_bB <<- -1*solve(cache.M$M_bb., M_ZX.M_XX)
                         } else { ## all random effects dropped
                             cache.M$M_bb. <<- Lambda_bD_b
                             cache.M$M_XZ <<- M_XZ <- Matrix(0, p, q)
                             cache.M$M_bB <<- Matrix(0, q, p)
                         }
                         cache.M$M_BB <<- solve(M_XX, Diagonal(p) - M_XZ %*% cache.M$M_bB)
                         cache.M$M_bb <<- solve(cache.M$M_bb.)
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
                         M_XX <<- crossprod(sqrtD_e %*% .U_eX)
                         tmp1 <- sqrtD_e %*% .U_eZ
                         M_XZ0 <<- crossprod(sqrtD_e %*% .U_eX, tmp1)
                         M_ZZ0 <<- crossprod(tmp1)
                         M_XX.M_ZZ0 <<- solve(M_XX, M_XZ0)
                         M_ZX0M_XX.M_ZZ0 <<- crossprod(M_XZ0, M_XX.M_ZZ0)
                         Epsi2_e <<- rho_e@Epsi2()
                         ## calculate Epsi_bbt
                         Epsi_bbt <<- bdiag(btapply(rho_b, .calcE.psi_bbt, rep=object@ind))
                         ## calculate Epsi_bpsi_bt
                         Epsi_bpsi_bt <<- bdiag(btapply(rho_b, .calcE.psi_bpsi_bt, rep=object@ind))
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
                     J = function() {
                         ## if (!isTRUE(calledInit))
                         ##     stop("Call initMatrices() first")
                         J <- Matrix(0, p + q, p + q)
                         J[1:p, 1:p] <- M_XX
                         if (any(!zeroB)) {
                             J[p+(1:q), 1:p] <-
                                 t(J[1:p, p+(1:q)] <- crossprod(crossprod(.U_eZ, D_e %*% X), U_b))
                             idx <- !zeroB
                             t.idx <- (p+(1:q))[idx]
                             La <- U_b[idx, idx]
                             J[t.idx, t.idx] <-
                                 crossprod(La, M_ZZ0[idx, idx] %*% La) +
                                     Lambda_bD_b[idx, idx]
                         }
                         return(J)
                     },
                     gh.int2d = function(fun, ...) {
                         lfun <- Vectorize(function(x, y, ...) fun(c(x, y), ...), c("x", "y"))
                         sum(lfun(ghZ, ghZt, ...)*ghW)
                     },
                     gh.int3d = function(fun, ...) {
                         lfun <- Vectorize(function(x, y, z, ...)
                                           fun(c(x, y, z)),
                                           c("x", "y"))
                         int1 <- function(z, ...) sum(lfun(ghZ, ghZt, z=z, ...)*ghW)
                         drop(crossprod(sapply(ghz, int1, ...),ghw))
                     },
                     gh.int4d = function(fun, ...) {
                         lfun <- Vectorize(function(x, y, z, u, ...)
                                           fun(c(x, y, z, u)),
                                           c("x", "y"))
                         int1 <- Vectorize(function(u, z, ...)
                                           sum(lfun(ghZ, ghZt, z=z, u=u, ...)*ghW),
                                           c("u", "z"))
                         sum(int1(ghZ, ghZt, ...)*ghW)
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
                                 tcrossprod(solve(M_XX, t(sqrtD_e %*% .U_eX)))
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
                             ## ## warn if default isSymmetric fails
                             ## if (!isSymmetric(cache.unsc)) {
                             ##     tol <- eval(formals(isSymmetric.matrix)$tol)
                             ##     warning("isSymmetric() failed: ",
                             ##             all.equal(cache.unsc, t(cache.unsc), tolerance=tol))
                             ## }
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
            list(A       = "Matrix",        ## Matrix A
                 Kt      = "Matrix",        ## Matrix B = lfrac * Kt, K = t(Kt)
                 L       = "Matrix",        ## Matrix L
                 kappa_e = "numeric",       ## kappa_e^(sigma) (only required for DASvar and DAStau)
                 kappa_b = "numeric",       ## kappa_b^(sigma) (-------------- " ------------------)
                 EDpsi_e = "numeric",
                 method  = "character",
                 .setTau_e = "logical",     ## check whether we already have computed tau_e
                 .tau_e  = "numeric",       ## tau_e used in updateSigma
                 .setTbk = "logical",
                 .Tbk    = "Matrix",        ## cache for T_{b,k}
                 blocks  = "list",
                 idx     = "list"
                 ),
            contains = "rlmerPredD",
            methods =
            list(
                 initialize = function(...) {
                     callSuper(...)
                     .setTau_e <<- FALSE
                     .tau_e <<- 0
                     .setTbk <<- FALSE
                     .Tbk <<- Matrix(0, 0, 0)
                     ## the other slots are initialized in initMatrices()
                 },
                 initMatrices = function(object) {
                    callSuper(object)
                    blocks <<- object@blocks
                    idx <<- object@idx
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
                 Tb = function() {
                     tmp <- L %*% Epsi_bbt
                     Tfull <- diag(q) - tmp - t(tmp) + Epsi2_e * crossprod(Kt) + L %*% crossprod(Epsi_bpsi_bt, L)
                     Ts <- list()
                     ## cycle blocks
                     for (type in seq_along(blocks)) {
                       bidx <- idx[[type]]
                       for (k in 1:ncol(bidx)) ## 1:K
                         Ts <- c(Ts, list(Tfull[bidx[,k],bidx[,k]]))
                     }
                     return(bdiag(Ts))
                 },
                 setT = function(T) {
                     .Tbk <<- T
                     .setTbk <<- TRUE
                 },
                 T = function() {
                     if (.setTbk) .Tbk else Tb() ## fallback to Tb()
                 },
                 updateMatrices = function() {
                     if (!isTRUE(calledInit)) initMatrices()
                     if (any(!zeroB)) {
                         r <- M()
                         tmp2 <- tcrossprod(.U_eX, ## X M_Bb U_b Z U_e\inv
                                            crossprod(U_btZt.U_et, r$M_bB))
                         tmp3 <- crossprod(U_btZt.U_et, r$M_bb) ## U_e\inv Z U_b M_bb

                         A <<- tcrossprod(.U_eX %*% r$M_BB, .U_eX) +
                             tmp2 + t(tmp2) + tmp3 %*% U_btZt.U_et
                         Kt <<- -1*(tcrossprod(.U_eX, r$M_bB) + tmp3)
                         L <<- r$M_bb %*% Lambda_b
                     } else {
                         ## no random effects
                         A <<- .U_eX %*% solve(M_XX, t(.U_eX)) ## just the hat matrix
                         Kt <<- Matrix(0, n, q)
                         L <<- solve(D_b)
                     }
                 },
                 tau_e = function() {
                     if (isTRUE(.setTau_e)) return(.tau_e)
                     .tau_e <<- calcTau(diag(A), .s(theta=FALSE, pp=.self),
                                        rho_e, rho_sigma_e, .self, kappa_e,
                                        .tau_e, method)
                     .setTau_e <<- TRUE
                     return(.tau_e)
                 }
                )
            )

## This is basically a copy of the rmerPredD-class,
## but with fields replaced by functions.
##
## @title rlmerPredD_Rcpp
## @slot all see rlmerPredD class
##' @importClassesFrom Rcpp Module
rlmerPredD_Rcpp <-
  setRefClass("rlmerPredD_Rcpp",
              fields =
                list(input = "list",
                     n = "numeric",
                     p = "numeric",
                     q = "numeric",
                     X = "matrix",
                     Zt = "dgCMatrix",
                     method = "character",
                     Ptr = "refClass"
                ),
              methods =
                list(
                  initialize = function(...) {
                    ## TODO test that we're not creating unnecessary copies!
                    args <- list(...)
                    input <<- args$input
                    X <<- as(args$X, "matrix")
                    Zt <<- as(args$Zt, "dgCMatrix")
                    n <<- nrow(X)
                    p <<- ncol(X)
                    q <<- nrow(Zt)
                    method <<- args$method
                    initializePtr()
                  },
                  copy = function(shallow = FALSE) {
                    def <- .refClassDef
                    selfEnv <- as.environment(.self)
                    vEnv    <- new.env(parent=emptyenv())

                    for (field in setdiff(names(def@fieldClasses), "Ptr")) {
                      if (shallow)
                        assign(field, get(field, envir = selfEnv), envir = vEnv)
                      else {
                        current <- get(field, envir = selfEnv)
                        if (is(current, "envRefClass"))
                          current <- current$copy(FALSE)
                        ## FIXME is this enough to make sure elements are copied as well?
                        ## (e.g., beta inside input list?)
                        assign(field, forceCopy(current), envir = vEnv)
                      }
                    }
                    do.call(rlmerPredD_Rcpp$new, c(as.list(vEnv), n=nrow(vEnv$V), Class=def))
                  },
                  ptr = function() {
                    'returns the external pointer, regenerating if necessary'
                    if (is.null(Ptr) || isnull(Ptr$.pointer)) initializePtr()
                    Ptr
                  },
                  initializePtr = function() {

                    rho_e_instance <- input$rho_e_src@getInstanceWithOriginalDefaults();
                    rhoSigma_e_instance <- input$rhoSigma_e_src@getInstanceWithOriginalDefaults()
                    rho_b_instance <- lapply(input$rho_b_src, function(x) x@getInstanceWithOriginalDefaults())
                    rhoSigma_b_instance <- lapply(input$rhoSigma_b_src, function(x) x@getInstanceWithOriginalDefaults())

                    input$rho_e_instance <<- rho_e_instance
                    input$rho_e <<- rho_e_instance$.pointer
                    input$rhoSigma_e_instance <<- rhoSigma_e_instance
                    input$rhoSigma_e <<- rhoSigma_e_instance$.pointer
                    input$rho_b_instance <<- rho_b_instance
                    input$rho_b <<- lapply(rho_b_instance, function(x) x$field(".pointer"))
                    input$rhoSigma_b_instance <<- rhoSigma_b_instance
                    input$rhoSigma_b <<- lapply(rhoSigma_b_instance, function(x) x$field(".pointer"))

                    Ptr <<- switch(tolower(method),
                                   dasvar = new(rlmerPredD_DAS, input, X, Zt),
                                   dastau = new(rlmerPredD_DAStau, input, X, Zt),
                                   stop("Unknown method ", method, ". Valid method arguments are: DASvar, DAStau."))
                    Ptr$initMatrices()
                  },
                  setSigma = function(value) {
                    ptr()$setSigma(value)
                  },
                  setU = function(value) {
                    ptr()$setU(value)
                  },
                  setB_s = function(value) {
                    ptr()$setB_s(value)
                  },
                  setB_r = function(value) {
                    ptr()$setB_r(value)
                  },
                  setBeta = function(value) {
                    ptr()$setBeta(value)
                  },
                  setTheta = function(value) {
                    ptr()$setTheta(value)
                  },
                  Lambdat = function() {
                    ptr()$Lambdat()
                  },
                  M = function() {
                    ptr()$M()
                  },
                  unsc = function() {
                    ptr()$unsc()
                  },
                  b = function() {
                    ptr()$b()
                  },
                  zeroB = function() {
                    ptr()$zeroB()
                  },
                  U_b = function() {
                    ptr()$U_b()
                  },
                  U_e = function() {
                    ptr()$U_e()
                  },
                  V_e = function() {
                    ptr()$V_e()
                  },
                  theta = function() {
                    ## take from input?
                    ptr()$theta()
                  },
                  beta = function() {
                    ## take from input?
                    ptr()$beta()
                  },
                  b_s = function() {
                    ## take from input?
                    ptr()$b_s()
                  },
                  effects = function() {
                    ptr()$effects()
                  },
                  sigma = function() {
                    ptr()$sigma()
                  },
                  M_XZ0 = function() {
                    ptr()$M_XZ0()
                  },
                  M_ZZ0 = function() {
                    ptr()$M_ZZ0()
                  },
                  M_ZZ0_sub_M_ZX0invM_XXMZZ0 = function() {
                    ptr()$M_ZZ0_sub_M_ZX0invM_XXMZZ0()
                  },
                  Lambda_b = function() {
                    ptr()$Lambda_b()
                  },
                  Lambda_bD_b = function() {
                    ptr()$Lambda_bD_b()
                  },
                  Epsi_bbt = function() {
                    ptr()$Epsi_bbt()
                  },
                  Epsi_bpsi_bt = function() {
                    ptr()$Epsi_bpsi_bt()
                  },
                  Epsi2_b = function() {
                    ptr()$Epsi2_b()
                  },
                  invU_btZtU_et = function() {
                    ptr()$invU_btZtU_et()
                  },
                  invU_eX = function() {
                    ptr()$invU_eX()
                  },
                  mu = function() {
                    ptr()$mu()
                  },
                  distB = function() {
                    ptr()$distB()
                  },
                  bBlockMap = function() {
                    ptr()$bBlockMap()
                  },
                  ## methods from rlmerPredD_DAS class:
                  s_e = function() {
                    ptr()$s_e()
                  },
                  S_b = function() {
                    ptr()$S_b()
                  },
                  kappa_e = function() {
                    ptr()$kappa_e()
                  },
                  kappa_b = function() {
                    ptr()$kappa_b()
                  },
                  Tb = function() {
                    ptr()$Tb()
                  },
                  tau_e = function() {
                    ptr()$tau_e()
                  },
                  A = function() {
                    ptr()$A()
                  },
                  Kt = function() {
                    ptr()$Kt()
                  },
                  L = function() {
                    ptr()$L()
                  }
                )
  )
rlmerPredD_Rcpp$lock("n", "p", "q", "X", "Zt", "method")

## This is basically a copy of the lmerResp-class
##
## REVISION: 1611
##
## @title rlmerResp_Rcpp
## @name rlmerResp_Rcpp-class
## @slot all see lmerResp class
rlmerResp_Rcpp <-
  setRefClass("rlmerResp_Rcpp",
              fields =
                list(Ptr     = "refClass",
                     mu      = "numeric",
                     offset  = "numeric",
                     ##sqrtXwt = "numeric",
                     sqrtrwt = "numeric",
                     weights = "numeric",
                     wtres   = "numeric",
                     y       = "numeric"),
              methods =
                list(initialize = function(...) {
                  if (!nargs()) return()
                  ll <- list(...)
                  if (is.null(ll$y)) stop("y must be specified")
                  y <<- as.numeric(ll$y)
                  n <- length(y)
                  mu <<- if (!is.null(ll[["mu"]]))
                    as.numeric(ll[["mu"]]) else numeric(n)
                  offset <<- if (!is.null(ll$offset))
                    as.numeric(ll$offset) else numeric(n)
                  weights <<- if (!is.null(ll$weights))
                    as.numeric(ll$weights) else rep.int(1,n)
                  ## sqrtXwt <<- if (!is.null(ll$sqrtXwt))
                  ##   as.numeric(ll$sqrtXwt) else sqrt(weights)
                  sqrtrwt <<- if (!is.null(ll$sqrtrwt))
                    as.numeric(ll$sqrtrwt) else sqrt(weights)
                  wtres   <<- sqrtrwt * (y - mu)
                  initializePtr()
                },
                copy         = function(shallow = FALSE) {
                  def <- .refClassDef
                  selfEnv <- as.environment(.self)
                  vEnv    <- new.env(parent=emptyenv())
                  for (field in setdiff(names(def@fieldClasses), "Ptr")) {
                    if (shallow)
                      assign(field, get(field, envir = selfEnv), envir = vEnv)
                    else {
                      current <- get(field, envir = selfEnv)
                      if (is(current, "envRefClass"))
                        current <- current$copy(FALSE)
                      ## deep-copy hack +0
                      assign(field, forceCopy(current), envir = vEnv)
                    }
                  }
                  do.call(new, c(as.list(vEnv), Class=def))
                },
                ptr = function() {
                  'returns the external pointer, regenerating if necessary'
                  if (is.null(Ptr) || isnull(Ptr$.pointer)) initializePtr()
                  Ptr
                },
                initializePtr = function() {
                  Ptr <<- new(rlmerResp, y, weights, offset, mu, ## sqrtXwt,
                              sqrtrwt, wtres)
                  Ptr$updateMu(mu - offset)
                },
                allInfo = function() {
                  'return all the information available on the object'
                  data.frame(y=y, offset=offset, weights=weights, mu=mu,
                             rwt=sqrtrwt, wres=wtres)
                },
                # setOffset  = function(oo) {
                #   'change the offset in the model (used in profiling)'
                #   ptr()$setOffset(as.numeric(oo))
                # },
                # setResp    = function(rr) {
                #   'change the response in the model, usually after a deep copy'
                #   ptr()$setResp(as.numeric(rr))
                # },
                # setWeights = function(ww) {
                #   'change the prior weights in the model'
                #   ptr()$setWeights(as.numeric(ww))
                # },
                updateMu  = function(gamma) {
                  'update mu and wtres from the linear predictor'
                  ptr()$updateMu(gamma)
                })
  )

rlmerResp_Rcpp$lock("mu", "offset", ## "sqrtXwt",
                    "sqrtrwt", "weights", "wtres")#, "y")

##' Class "rlmerMod" of Robustly Fitted Mixed-Effect Models
##'
##' A robust mixed-effects model as returned by \code{\link{rlmer}}.
##' @title rlmerMod Class
##' @name rlmerMod-class
##' @aliases rlmerMod-class coef.rlmerMod deviance.rlmerMod
##' extractAIC.rlmerMod family.rlmerMod fitted.rlmerMod fixef.rlmerMod
##' formula.rlmerMod isGLMM.rlmerMod isLMM.rlmerMod isNLMM.rlmerMod
##' isREML.rlmerMod logLik.rlmerMod model.frame.rlmerMod
##' model.matrix.rlmerMod nobs.rlmerMod predict.rlmerMod
##' print.rlmerMod print.summary.rlmer print.VarCorr.rlmerMod
##' ranef.rlmerMod resid.rlmerMod rlmerMod-class sigma.rlmerMod
##' show.rlmerMod show,rlmerMod-method show.summary.rlmerMod
##' summary.rlmerMod summary.summary.rlmerMod terms.rlmerMod
##' update.rlmerMod VarCorr.rlmerMod VarCorr.summary.rlmerMod
##' vcov.rlmerMod vcov.summary.rlmerMod weights.rlmerMod
##' @docType class
##' @section Objects from the Class: Objects are created by calls to
##' \code{\link{rlmer}}.
##' @section Methods: Almost all methods available from objects
##' returned from \code{\link{lmer}} are also available for objects
##' returned by \code{\link{rlmer}}. They usage is the
##' same.
##'
##' It follows a list of some the methods that are exported by this
##' package:
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
##' @section Disabled methods: A log likelihood or even a pseudo log
##' likelihood is not defined for the robust estimates returned by
##' \code{\link{rlmer}}. Methods that depend on the log likelihood are
##' therefore not available. For this reason the methods
##' \code{deviance}, \code{extractAIC} and \code{logLik} stop with an
##' error if they are called.
##' @seealso \code{\link{rlmer}}; corresponding class in package
##' \code{lme4}: \code{\link{merMod}}
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
         representation(resp    = "refClass",
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

## setIs(class1 = "rlmerPredD", class2 = "merPredD",
##       coerce = function(from) {
##           theta <- from$theta
##           if (any(idx <- theta == 0)) theta[idx] <- 1
##           pp <- new("merPredD", X=from$X, Zt=from$Zt,
##                     Lambdat=from$Lambdat(),
##                     Lind=from$.Lind,
##                     theta=theta, n=from$n, beta0=from$beta,
##                     u0=from$b.s)
##           if (any(idx)) pp$setTheta(from$theta)
##           pp
##       },
##       replace = function(from, value) {
##           if (!missing(value)) {
##               from$X <- value$X
##               from$Zt <- value$Zt
##               from$setLambdat(value$Lambdat, value$Lind)
##               from$theta <- value$theta
##               from$n <- nrow(value$V)
##               from$beta <- value$delb
##               from$b.s <- value$delu
##               warning("Now run updatePp on the rlmerMod object.")
##           }
##           from
##       })

## ## @name rlmerMod-class
## ### Define inheritance from lmerMod to rlmerMod
## setIs(class1 = "rlmerMod", class2 = "lmerMod",
##       coerce = function(from) {
##           ##cat("~~~~ coerce rlmerMod to lmerMod ~~~~~\n")
##           new("lmerMod",
##               resp=from@resp, Gp = from@Gp,
##               call=from@call, frame=from@frame,
##               flist=from@flist, cnms=from@cnms,
##               lower=from@lower, theta=from@theta,
##               beta=from@beta, u=from@b.s,
##               devcomp=from@devcomp,
##               pp=as(from@pp, "merPredD"))
##       },
##       replace = function(from, value) {
##           ##cat("~~~~ replace lmerMod with rlmerMod ~~~~~\n")
##           if (! missing(value)) {
##               from@resp <- value@resp
##               from@Gp <- value@Gp
##               from@call <- value@call
##               from@frame <- value@frame
##               from@flist <- value@flist
##               from@cnms <- value@cnms
##               from@lower <- value@lower
##               from@theta <- value@theta
##               from@beta <- value@beta
##               from@b.s <- value@u
##               from@devcomp = value@devcomp
##               dd <- value@devcomp$dims
##               cmp <- value@devcomp$cmp
##               from@pp <- new("rlmerPredD",
##                              X=value@pp$X,Zt=value@pp$Zt,
##                              Lambdat=value@pp$Lambdat, Lind=value@pp$Lind,
##                              theta=value@theta, beta=value@beta, b.s=value@u,
##                              lower=value@lower, v_e=value@resp$weights,
##                              sigma=cmp[[ifelse(dd["REML"], "sigmaREML", "sigmaML")]],
##                              numpoints = 13)
##               from@rho.e <- cPsi
##               from@rho.sigma.e = from@rho.e
##               from@method <- "DAStau"
##               from@b.r <- from@pp$b.r
##               b <- findBlocks(value@pp)
##               from@rho.b <- rep(list(cPsi),length(b$dim))
##               from@rho.sigma.b = from@rho.b
##               from@blocks <- b$blocks
##               from@ind <- b$ind
##               from@idx <- b$idx
##               from@dim <- b$dim
##               from@q <- b$q
##               from@k <- b$k
##               from@pp$initMatrices(from)
##           }
##           from
##       })

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

    resp <- new("rlmerResp_Rcpp", y = y, offset = offset,
                weights = weights, mu = mu)
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
    if (inherits(object@pp, "rlmerPredD")) {
        object@pp$lower <- object@lower
    } else if (!identical(object@pp$lower, object@lower)) {
          stop("Can't set lower value in Rcpp implementation at the moment.")
    }
    dd <- object@devcomp$dims
    object@pp$setSigma(object@devcomp$cmp[[ifelse(dd["REML"], "sigmaREML", "sigmaML")]])
    object@pp$setB.s(object@b.s)
    invisible(object)
}
