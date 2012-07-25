##' define new class psi_func2 which is basically
##' just psi_func with added slot name
##'
##' @title psi_func2 class
##' @exportClass psi_func2 
setClass("psi_func2", contains = c("psi_func"),
         representation(name = "character"))

##' @rdname psi_func2
##' @exportClass psi_func2_cached
setClass("psi_func2_cached", contains = c("psi_func2"))


##' This is basically a copy of the merPredD-class
##'
##' The one thing missing is Ptr (and the lock).
##'
##' REVISION: 1611
##' 
##' @title rlmerPredD
##' @name rlmerPredD-class
##' @slot all see rlmerPredD class
setRefClass("rlmerPredD",
            fields =
                list(U_e     = "ddiMatrix",
                     U_b     = "sparseMatrix",
                     .Lambdat = "dgCMatrix",
                     .Lind    = "integer", ## for .Lambdat
                     Lind    = "integer", ## for U_b
                     X       = "matrix",
                     Zt      = "dgCMatrix",
                     beta    = "numeric",
                     b.s     = "numeric",
                     b       = "numeric",                    
                     theta   = "numeric",
                     sigma   = "numeric",
                     deviance = "numeric",
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
                     lfrac   = "numeric",      ## \lambda_e / \lambda_b
                     sqrtD_e = "ddiMatrix",    ## sqrt(D_e)
                     .U_eX   = "Matrix",       ## U_e\inv X
                     .U_eZ   = "Matrix",       ## U_e\inv Z
                     M_XX    = "Matrix",       ## M_XX
                     U_btZt.U_et = "Matrix",   ## U_b\tr Z\tr U_e\inv\tr
                     Echi    = "numeric",      ## used in .s for calculating var(g)
                     ## remlK   = "Matrix",       ## used to calculate deviance in REML case
                     ## remlKcorr = "numeric",    ## ditto
                     RZX     = "matrix",       ## used to calculate deviance in REML case
                     calledInit = "logical",   ## whether initMatrices has been called yet
                     M_XZ0   = "Matrix",       ## M_XZ = M_XZ0 %*% U_b
                     M_ZZ0   = "Matrix",       ## M_ZZ = U_b\tr %*% M_ZZ0 %*% U_b
                     M_XX.M_ZZ0 = "Matrix",    ## M_XX^-1 %*% M_XZ0
                     M_ZX0M_XX.M_ZZ0 = "Matrix", ## crossprod(M_XZ0, M_XX^M_ZZ0)
                     Epsi2_e = "numeric",      ## \Erw[\psi^2_e]
                     Epsi2_b = "numeric",      ## \Erw[\psi^2_b]
                     Epsi_bpsi_bt = "Matrix",   ## \Erw[\psi_b(u)\psi_b(u)\tr]
                     set.M   = "logical",      ## cache
                     cache.M = "list",         ## cache
                     set.unsc = "logical",     ## cache
                     cache.unsc = "matrix"        ## cache
                     ),
                methods =
                list(
                     initialize =
                     function(X, Zt, RZX, Lambdat, Lind, theta, beta, b.s, lower,
                              sigma, deviance, u_e, ...) {
                         set.unsc <<- set.M <<- calledInit <<- FALSE
                         initGH(...)
                         if (!nargs()) return()
                         X <<- as(X, "matrix")
                         Zt <<- as(Zt, "dgCMatrix")
                         RZX <<- as(RZX, "matrix")
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
                                   length(deviance) == 1,
                                   length(lower) == length(theta),
                                   length(u_e) == n
                                   )
                         U_e <<- Diagonal(x=if (is(u_e, "numeric")) u_e else diag(u_e))
                         beta <<- beta
                         b.s <<- b.s
                         b <<- as(U_b %*% b.s, "numeric")
                         sigma <<- sigma
                         deviance <<- deviance
                         lower <<- lower
                         setZeroB()
                     },
                     setSigma = function(value) {
                         sigma <<- value
                     },
                     setDeviance = function(value) {
                         deviance <<- value
                     },
                     setU = function(value) setB.s(value),
                     setB.s = function(value) {
                         b.s <<- value
                         b <<- as(U_b %*% value, "numeric")
                     },
                     setB = function(value) {
                         b <<- value
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
                        ## initialize often required matrices
                         .U_eX <<- solve(U_e, X)
                         .U_eZ <<- solve(U_e, t(Zt))
                         U_btZt.U_et <<- t(.U_eZ %*% U_b)
                         ## if (.isREML(object)) {
                         ##     ## generate matrix to be used for REML
                         ##     ## P %*% X == 0, but has too many rows
                         ##     H <- .U_eX %*% solve(M_XX, t(.U_eX))
                         ##     P <- Diagonal(n) - D_e %*% H
                         ##     qrP <- qr(P)
                         ##     remlK <<- as(t(qr.Q(qrP)[,1:qrP$rank]), "Matrix") ## K %*% X == 0
                         ##     remlKcorr <<- -2 * sum(log(abs(diag(qr.R(qrP))[1:qrP$rank])))
                         ## }
                         initRho(object)
                         calledInit <<- TRUE
                     },
                     M = function() {
                         if (set.M) return(cache.M)
                         cache.M <<- list()
                         if (any(!zeroB)) {
                             ## M_bb. := M_bb\inv
                             cache.M$M_bb. <<- crossprod(U_b,(M_ZZ0 - M_ZX0M_XX.M_ZZ0)) %*% U_b +
                                 lfrac * D_b
                             cache.M$M_XZ <<- M_XZ <- M_XZ0 %*% U_b
                             M_ZX.M_XX <- t(solve(M_XX, M_XZ))
                             cache.M$M_bB <<- -1*solve(cache.M$M_bb., M_ZX.M_XX)
                             cache.M$M_BB <<- solve(M_XX, Diagonal(p) - M_XZ %*% cache.M$M_bB)
                             cache.M$M_bb <<- solve(cache.M$M_bb.)
                         }
                         set.M <<- TRUE
                         return(cache.M)
                     },
                     initRho = function(object) {
                         D_e <<- Diagonal(x=rep(object@rho.e@EDpsi(), n))
                         tmp <- if (object@wExp.b == 0) rep(object@rho.b@EDpsi(),q) else {
                             exps <- sapply(object@dim, .calcE.D.re, rho = object@rho.b)
                             exps[object@ind[object@k]]
                         }
                         D_b <<- Diagonal(x=tmp)
                         lfrac <<- object@rho.e@EDpsi() / object@rho.b@EDpsi()
                         sqrtD_e <<- Diagonal(x=sqrt(D_e@x))
                         Echi <<- if (object@wExp.b == 0) rep(object@rho.b@Epsi2(), q) else {
                             exps <- sapply(object@dim, .calcEchi, rho = object@rho.b)
                             exps[object@ind[object@k]]
                         }
                         M_XX <<- crossprod(sqrtD_e %*% .U_eX)
                         tmp1 <- sqrtD_e %*% .U_eZ
                         M_XZ0 <<- crossprod(sqrtD_e %*% .U_eX, tmp1)
                         M_ZZ0 <<- crossprod(tmp1)
                         M_XX.M_ZZ0 <<- solve(M_XX, M_XZ0)
                         M_ZX0M_XX.M_ZZ0 <<- crossprod(M_XZ0, M_XX.M_ZZ0)
                         Epsi2_e <<- object@rho.e@Epsi2()
                         Epsi2_b <<- object@rho.b@Epsi2()
                         ## calculate Epsi_bpsi_bt
                         Eblks <- list()
                         for (s_k in object@dim) {
                             if (s_k == 1) {
                                 tmp <- Matrix(Epsi2_b)
                             } else if (s_k == 2) {
                                 wgt <- .wgtTau(object@rho.b, object@wExp.b)
                                 A <- gh.int2d(function(b) wgt(sqrt((b[1]*b[1]+b[2]*b[2])/2))*b[1]*b[1])
                                 B <- gh.int2d(function(b) wgt(sqrt((b[1]*b[1]+b[2]*b[2])/2))*b[1]*b[2])
                                 tmp <- Matrix(c(A, B, B, A), 2)
                             } else {
                                 stop("block size > 2 not implemented yet")
                                 ## FIXME: like s_k == 2 but with third term of
                                 ##        chi^2 distribution with s_k - 2 df
                             }
                             Eblks <- c(Eblks, tmp)
                         }
                         Epsi_bpsi_bt <<- bdiag(Eblks[object@ind])
                     },
                     initGH = function(numpoints=13) {
                         gh <- robustbase:::ghq(numpoints)
                         ghz <<- gh$nodes
                         ghw <<- gh$weights*dnorm(gh$nodes)
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
                                     lfrac * D_b[idx, idx]
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
                         tmp <- if (is.null(r$M_bb)) { ## all theta == 0
                             Epsi2_e * tcrossprod(solve(M_XX, t(sqrtD_e %*% .U_eX)))
                         } else {
                             Epsi2_e * with(r, M_BB - crossprod(M_bB, lfrac * D_b %*% M_bB)) +
                                 lfrac^2 * 
                                 if (isDiagonal(U_b)) {
                                     Epsi2_b * crossprod(r$M_bB)
                                 } else {
                                     with(r, crossprod(M_bB, Epsi_bpsi_bt %*% M_bB))
                                 }
                         }
                         cache.unsc <<- as.matrix(tmp)
                         ## test matrix for symmetry
                         ## if its almost symmetric, then return symmetrized matrix
                         ## otherwise summary.merMod() and chol() will complain.
                         if (!isSymmetric(cache.unsc, 0)) {
                             ## warn if default isSymmetric fails
                             if (!isSymmetric(cache.unsc)) {
                                 tol <- eval(formals(isSymmetric.matrix)$tol)
                                 warning("isSymmetric() failed: ",
                                         all.equal(cache.unsc, t(cache.unsc), tolerance=tol))
                             }
                             cache.unsc <<- symmpart(cache.unsc)
                         }
                         set.unsc <<- TRUE
                         return(cache.unsc)
                     }
                     )
                )

setRefClass("rlmerPredD_DAS",
            fields =
            list(A       = "Matrix",        ## Matrix A
                 Kt      = "Matrix",        ## Matrix B = lfrac * Kt, K = t(Kt)
                 L       = "Matrix"         ## Matrix L
                 ),
            contains = "rlmerPredD",
            methods =
            list(
                 initialize = function(...) {
                     callSuper(...)
                     ## the other slots are initialized in initMatrices()
                 },
                 setTheta = function(value) {
                     callSuper(value)                     
                     ## update Matrices
                     updateMatrices()
                 },
                 B = function() {
                     lfrac * Kt
                 },
                 K = function() {
                     t(Kt)
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
                         L <<- lfrac * r$M_bb
                     } else {
                         ## no random effects
                         A <<- .U_eX %*% solve(M_XX, t(.U_eX)) ## just the hat matrix
                         Kt <<- Matrix(0, n, q)
                         L <<- Matrix(0, q, q) ## FIXME: shouldn't it be a diag 1 matrix?
                     }
                 })
            )

##' This is similar to lmerMod from lme4
##'
##' @title rlmerMod
##' @name rlmerMod-class
##' @slot pp rlmerPredD class
##' @slot cache environment used as cache.
##' @slot others as in lmerMod class
setClass("rlmerMod",
         representation(resp    = "lmerResp",
                        Gp      = "integer",
                        call    = "call",
			frame   = "data.frame", # "model.frame" is not S4-ized yet
                        flist   = "list",
                        cnms    = "list",
                        lower   = "numeric",
                        theta   = "numeric",
                        beta    = "numeric",
                        b.s       = "numeric",
                        devcomp = "list",
                        pp      = "rlmerPredD",
                        cache   = "environment",
                        ## from rlmerResp:
                        rho.e   = "psi_func2",
                        wExp.e  = "numeric",
                        rho.sigma.e = "psi_func2",
                        use.laplace = "logical",
                        method  = "character",
                        method.effects = "character",
                        ## from rreTrms:
                        b       = "numeric",
                        rho.b   = "psi_func2",
                        rho.sigma.b = "psi_func2",
                        wExp.b  = "numeric",
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
             if (length(object@wExp.e) != 1 || object@wExp.e < 0)
                 return("wExp.e has to be positive and be of length 1")
             v <- validObject(object@rho.sigma.e)
             if (!is.logical(v) || ! v)
                 return(v)
             v <- validObject(object@rho)
             if (!is.logical(v) || ! v)
                 return(v)
             if (length(object@wExp.b) != 1 || object@wExp.b < 0)
                 return("wExp.b has to be positive and be of length 1")
             v <- validObject(object@rho.sigma.b)
             if (!is.logical(v) || ! v)
                 return(v)
             if (length(object@b) != length(object@b.s)) 
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
             if (length(object@k) != length(object@b.S))
                 return("length of k must be identical to length of u")
             if (length(object@ind) != max(object@k))
                 return("length of ind must be equal to max k")
             if (max(object@k) != sum(nblocks) || any(object@k <= 0))
                 return("block index in k must be positive and not exceed the number of blocks")
             if (any(unlist(lapply(object@blocks, function(x) unique(as.vector(x)))) !=
                     unique(unlist(lapply(object@blocks, as.vector)))))
                 return("parameters are not allowed to be in multiple blocks simultaneously")
             TRUE
         },
         S3methods = TRUE)

#######################################################
## Coercing methods                                  ##
#######################################################

setIs(class1 = "rlmerPredD", class2 = "merPredD",
      coerce = function(from) {
          theta <- from$theta
          if (any(idx <- theta == 0)) theta[idx] <- 1
          pp <- new("merPredD", X=from$X, Zt=from$Zt,
                    Lambdat=from$Lambdat(), RZX=from$RZX,
                    Lind=from$.Lind,
                    theta=theta, n=from$n, beta0=from$beta,
                    u0=from$b.s)
          if (any(idx)) pp$setTheta(from$theta)
          pp          
      },
      replace = function(from, value) {
          if (!missing(value)) {
              from$X <- value$X
              from$Zt <- value$Zt
              from$setLambdat(value$Lambdat, value$Lind)
              from$theta <- value$theta
              from$n <- nrow(value$V)
              from$beta <- value$delb
              from$b.s <- value$delu
              warning("Now run updatePp on the rlmerMod object.")
          }
          from
      })

##' @name rlmerMod-class
### Define inheritance from lmerMod to rlmerMod
setIs(class1 = "rlmerMod", class2 = "lmerMod",
      coerce = function(from) {
          ##cat("~~~~ coerce rlmerMod to lmerMod ~~~~~\n")
          new("lmerMod",
              resp=from@resp, Gp = from@Gp,
              call=from@call, frame=from@frame,
              flist=from@flist, cnms=from@cnms,
              lower=from@lower, theta=from@theta,
              beta=from@beta, u=from@b.s,
              devcomp=from@devcomp,
              pp=as(from@pp, "merPredD"))
      },
      replace = function(from, value) {
          ##cat("~~~~ replace lmerMod with rlmerMod ~~~~~\n")
          if (! missing(value)) {
              from@resp <- value@resp
              from@Gp <- value@Gp
              from@call <- value@call
              from@frame <- value@frame
              from@flist <- value@flist
              from@cnms <- value@cnms
              from@lower <- value@lower
              from@theta <- value@theta
              from@beta <- value@beta
              from@b.s <- value@u
              from@devcomp = value@devcomp
              dd <- value@devcomp$dims
              cmp <- value@devcomp$cmp
              from@pp <- new("rlmerPredD",
                             X=value@pp$X,Zt=value@pp$Zt, RZX=value@pp$RZX,
                             Lambdat=value@pp$Lambdat, Lind=value@pp$Lind,
                             theta=value@theta, beta=value@beta, b.s=value@u,
                             lower=value@lower, u_e=value@resp$weights,
                             sigma=cmp[[ifelse(dd["REML"], "sigmaREML", "sigmaML")]],
                             deviance=cmp[[ifelse(dd["REML"], "REML", "dev")]],
                             numpoints = 13)
              from@cache <- new.env()
              from@rho.e <- cPsi
              from@wExp.e <- 0
              from@rho.sigma.e = from@rho.e
              from@use.laplace <- TRUE ## default for Opt method
              from@method <- "Opt"
              from@method.effects <- "IRWLS"
              from@b <- from@pp$b
              from@rho.b <- cPsi
              from@wExp.b <- 0
              from@rho.sigma.b = from@rho.b
              b <- findBlocks(value@pp)
              from@blocks <- b$blocks
              from@ind <- b$ind
              from@idx <- b$idx
              from@dim <- b$dim
              from@q <- b$q
              from@k <- b$k
              from@pp$initMatrices(from)
          }
          from
      })

setAs("rlmerPredD", "rlmerPredD_DAS", function(from) {
    to <- new("rlmerPredD_DAS")
    fields <- names(getRefClass(class(from))$fields())
    Map(function(field) to$field(field, from$field(field)), fields)
    to
})

#######################################################
## Update methods                                    ##
#######################################################

##' Update robustness weights after changing any of the relevant
##' parts in a rlmerMod object
##'
##' @title Update robustness weights
##' @param object rlmerMod object to update
##' @return rlmerMod object
updateWeights <- function(object) {
    ## copy estimates from pp to object
    object@theta <- object@pp$theta
    object@beta <- object@pp$beta
    object@b.s <- object@pp$b.s
    object@b <- object@pp$b
    dd <- object@devcomp$dims
    object@devcomp$cmp[ifelse(dd["REML"], "sigmaREML", "sigmaML")] <- object@pp$sigma
    object@devcomp$cmp[ifelse(dd["REML"], "REML", "dev")] <- object@pp$deviance
    ## Set slots to NA
    object@devcomp$cmp[c("ldL2", "ldRX2", "wrss", "ussq", "pwrss", "drsum")] <- NA
    
    object
}

##' update pp slot from values in object
##'
##' (this does the reverse of updateWeights
##' 
##' @title set values in pp
##' @param object rlmerMod object to update
##' @return nothing
updatePp <- function(object) {
    object@pp$lower <- object@lower
    dd <- object@devcomp$dims
    object@pp$setSigma(object@devcomp$cmp[[ifelse(dd["REML"], "sigmaREML", "sigmaML")]])
    object@pp$setDeviance(object@devcomp$cmp[[ifelse(dd["REML"], "REML", "dev")]])
    object@pp$setB.s(object@b.s)
    invisible(object)
}
