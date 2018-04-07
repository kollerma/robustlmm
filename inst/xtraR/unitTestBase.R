if (!exists("unitTestBaseWasRun")) {
  unitTestBaseWasRun <- TRUE

  require(robustlmm)
  require(lme4)
  doBench <- require(microbenchmark)

  xtraR <- system.file("xtraR", package="robustlmm")
  source(file.path(xtraR, "DAS-scale-Rimpl.R"))

  clear <- function() rm("unitTestBaseWasRun", "fms", "objs",
                         "resps", "rPDASs", "rPDAStaus", "rPDASTests",
                         "rPDs", "rPDTests", envir = .GlobalEnv)

  cat("Starting...\n")

  set.seed(1)

  ########################################################
  ## old classes

  setRefClass("rlmerPredD_test",
              fields =
                list(U_e     = "ddiMatrix",
                     V_e     = "ddiMatrix",        ## crossprod(U_e) FIXME test this
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
                      V_e <<- Diagonal(x=if (is(v_e, "numeric")) v_e^2 else diag(v_e)^2)
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
                  initMatrices = function(list) {
                    if (isTRUE(calledInit)) return()
                    ## initialize often required matrices and other stuff
                    dim <<- list$dim
                    .U_eX <<- solve(U_e, X)
                    .U_eZ <<- solve(U_e, t(Zt))
                    U_btZt.U_et <<- t(.U_eZ %*% U_b)
                    initRho(list)
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
                  initRho = function(list) {
                    rho_e <<- list$rho_e_src
                    rho_sigma_e <<- list$rhoSigma_e_src
                    rho_b <<- list$rho_b_src
                    rho_sigma_b <<- list$rhoSigma_b_src
                    D_e <<- Diagonal(x=rep.int(rho_e@EDpsi(), n))
                    tmp <- btapply(rho_b, calcE.D.re2, rep=list$ind[list$k])
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
                    Epsi_bbt <<- bdiag(btapply(rho_b, calcE.psi_bbt, rep=list$ind))
                    ## calculate Epsi_bpsi_bt
                    Epsi_bpsi_bt <<- bdiag(btapply(rho_b, calcE.psi_bpsi_bt, rep=list$ind))
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

  setRefClass("rlmerPredD_DAS_test",
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
                     .Tbk     = "Matrix",       ## cache for T_{b,k}
                     blocks  = "list",
                     idx     = "list"
                ),
              contains = "rlmerPredD_test",
              methods =
                list(
                  initialize = function(...) {
                    callSuper(...)
                    .setTau_e <<- FALSE
                    .tau_e <<- 0
                    .setTbk <<- FALSE
                    .Tbk <<- Matrix(0, 0, 0)
                    method <<- "DASvar"
                    ## the other slots are initialized in initMatrices()
                  },
                  initMatrices = function(object) {
                    callSuper(object)
                    blocks <<- object$blocks
                    idx <<- object$idx
                  },
                  initRho = function(list) {
                    callSuper(list)
                    EDpsi_e <<- rho_e@EDpsi()
                    kappa_e <<- calcKappaTau(rho_sigma_e, 1)
                    kappa_b <<- calcKappaTauB(list, .self, list[["rhoSigma_b_src"]])
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

  setClass("rlmerMod_test",
           representation(resp    = "refClass",
                          pp      = "ANY",
                          cnms    = "list",
                          rho.e   = "psi_func_rcpp",
                          rho.sigma.e = "psi_func_rcpp",
                          rho.b   = "list",
                          rho.sigma.b = "list",
                          blocks  = "list",
                          ind     = "numeric",
                          idx     = "list",
                          dim     = "numeric",
                          q       = "numeric",
                          k       = "numeric"
           ))

  createInput <- function(fm, rhoFun) {
    Lambdat <- getME(fm, "Lambdat")
    Lind <- getME(fm, "Lind")
    b <- robustlmm:::findBlocks(Lambdat = Lambdat, Lind = Lind)
    ## convert to prop2 always here and use different tuning constants for rho_b.
    prop2 <- psi2propII(rhoFun, 1.8)
    rho_b = sapply(seq_along(b$blocks), function(i) chgDefaults(rhoFun, 1.2 * i))
    rho_sigma_b = sapply(seq_along(b$blocks), function(i) {
      if (NROW(b$blocks[[i]]) > 1) {
        return(chgDefaults(rhoFun, 1.5 * i))
      } else {
        return(psi2propII(chgDefaults(rhoFun, 1.5 * i)))
      }
      })
    rho_e_instance <- rhoFun@getInstanceWithOriginalDefaults();
    rhoSigma_e_instance <- prop2@getInstanceWithOriginalDefaults()
    rho_b_instance <- lapply(rho_b, function(x) x@getInstanceWithOriginalDefaults())
    rhoSigma_b_instance <- lapply(rho_sigma_b, function(x) x@getInstanceWithOriginalDefaults())
    input <- list(X = getME(fm, "X"),
                  Zt = getME(fm, "Zt"),
                  Lambdat = Lambdat,
                  lower = getME(fm, "lower"),
                  dim = as.integer(b$dim),
                  v_e = weights(fm),
                  Lind = as.integer(Lind),
                  ind = as.integer(b$ind),
                  k = as.integer(b$k),
                  idx = b$idx,
                  blockBMap = unlist(unlist(lapply(b$idx, apply, 2, list), recursive = FALSE), recursive = FALSE),
                  thetaBlockMap = lapply(b$blocks, function(b) b[lower.tri(b, diag = TRUE)]),
                  rho_e_src = rhoFun,
                  rho_e_instance = rho_e_instance,
                  rho_e = rho_e_instance$.pointer,
                  rhoSigma_e_src = prop2,
                  rhoSigma_e_instance = rhoSigma_e_instance,
                  rhoSigma_e = rhoSigma_e_instance$.pointer,
                  rho_b_src = rho_b,
                  rho_b_instance = rho_b_instance,
                  rho_b = lapply(rho_b_instance, function(x) x$field(".pointer")),
                  rhoSigma_b_src = rho_sigma_b,
                  rhoSigma_b_instance = rhoSigma_b_instance,
                  rhoSigma_b = lapply(rhoSigma_b_instance, function(x) x$field(".pointer")),
                  theta = getME(fm, "theta"),
                  beta = getME(fm, "beta"),
                  b.s = getME(fm, "u"),
                  sigma = getME(fm, "sigma"),
                  blocks = b$blocks,
                  idx = b$idx,
                  ## for DAS
                  maxOperations = 200,
                  relativeTolerance = 1e-6,
                  ## for DAStau
                  nodes = 13,
                  ## for use in tests below
                  b = b
    )
    return(input)
  }

  convToRlmerPredD <- function(fm, rhoFun) {
    cat(paste("Converting", deparse(fm@call$formula), "..."))
    input <- createInput(fm, rhoFun)
    cat("...created input...\n")
    obj <- new(robustlmm:::rlmerPredD, input, input[["X"]], input[["Zt"]])
    cat("created rlmerPredD")
    attr(obj, "input") <- input
    obj$initMatrices()
    obj
  }

  convToRlmerPredD_test <- function(rlmerPredD) {
    input <- attr(rlmerPredD, "input")
    obj <- new("rlmerPredD_test", input[["X"]], input[["Zt"]],
               input[["Lambdat"]], input[["Lind"]], input[["theta"]],
               input[["beta"]], input[["b.s"]], input[["lower"]],
               input[["sigma"]], input[["v_e"]])
    obj$initMatrices(input)
    obj
  }

  convToRlmerPredD_DAS <- function(fm, rhoFun) {
    input <- createInput(fm, rhoFun)
    obj <- new(robustlmm:::rlmerPredD_DAS, input, input[["X"]], input[["Zt"]])
    attr(obj, "input") <- input
    obj$initMatrices()
    obj
  }

  convToRlmerPredD_DAS_test <- function(rlmerPredD) {
    input <- attr(rlmerPredD, "input")
    obj <- new("rlmerPredD_DAS_test", input[["X"]], input[["Zt"]],
               input[["Lambdat"]], input[["Lind"]], input[["theta"]],
               input[["beta"]], input[["b.s"]], input[["lower"]],
               input[["sigma"]], input[["v_e"]])
    obj$initMatrices(input)
    obj$updateMatrices()
    obj
  }

  convToRlmerPredD_DAStau <- function(fm, rhoFun) {
    input <- createInput(fm, rhoFun)
    obj <- new(robustlmm:::rlmerPredD_DAStau, input, input[["X"]], input[["Zt"]])
    attr(obj, "input") <- input
    obj$initMatrices()
    obj
  }

  convToRlmerResp <- function(fm) {
    y <- getME(fm, "y")
    offset <- getME(fm, "offset")
    weights <- weights(fm)
    mu <- getME(fm, "mu")
    sqrtrwt <- sqrt(weights)
    wtres <- sqrtrwt * (y - mu)
    obj <- new(robustlmm:::rlmerResp, y, weights, offset, mu, ## sqrtXwt,
               sqrtrwt, wtres)
    obj
  }

  conv2rlmerObj_test <- function(i) {
    pp <- rPDASs[[i]]
    ppTest <- rPDASTests[[i]]
    b <- attr(pp, "input")[["b"]]
    new("rlmerMod_test",
        resp = resps[[i]],
        pp = pp,
        cnms = fms[[i]]@cnms,
        rho.e = ppTest$rho_e,
        rho.sigma.e = ppTest$rho_sigma_e,
        rho.b   = ppTest$rho_b,
        rho.sigma.b = ppTest$rho_sigma_b,
        blocks=b$blocks,
        ind=b$ind,
        idx=b$idx,
        dim=b$dim,
        q=b$q,
        k=b$k)
  }

  ########################################################
  ## some helper functions

  testForAllPsiFunc <- function(test) {
    test(cPsi)
    gc()
    test(huberPsiRcpp)
    gc()
    test(smoothPsi)
    gc()
  }

  testBlocks <- function(object, test, setTheta) {
    blocks <- robustlmm:::findBlocks(Lambdat = object@pp$Lambdat, Lind = object@pp$Lind)
    for (block in blocks$blocks) {
      thetaInd <- block[lower.tri(block, diag = TRUE)]
      theta <- abs(rnorm(length(getME(object, "theta"))))
      theta[thetaInd] <- 0
      setTheta(theta)
      test()
      theta[thetaInd] <- abs(rnorm(length(thetaInd)))
      setTheta(theta)
      test()
    }
  }
}
