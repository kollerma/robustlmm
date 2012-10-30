## test calculation of matrices for the classes
## rlmerPredD and rlmerPred_...
require(robustlmm)
##.calcEchi <- robustlmm:::.calcEchi
.calcE.D.re <- robustlmm:::.calcE.D.re
setTheta <- robustlmm:::setTheta
len <- robustlmm:::len

## old variants
## (basically a concatenation of initMatrices and
##  updateMatrices)
## should produce the same matrices A, B, K, L as
## the rlmerPred_DAS class
calcMatricesDiagOnly <- function(object, numpoints=13, envir=new.env()) {
    evalq({
        object <- get("object", envir = parent.frame(2))
        numpoints <- get("numpoints", envir = parent.frame(2))
        
        X <- object@pp$X
        Zt <- object@pp$Zt
        rho.resp <- object@rho.e
        rho.re <- object@rho.b
        lfrac <- object@rho.e@EDpsi() / object@rho.b@EDpsi()
        method <- object@method

        p <- object@pp$p
        q <- object@pp$q
        n <- object@pp$n

        ## D_e and D_b matrices
        D_e <- Diagonal(x=rep(object@rho.e@EDpsi(), n))
        tmp <- if (object@wExp.b == 0) rep(object@rho.b@EDpsi(),q) else {
            exps <- sapply(object@dim, .calcE.D.re, rho = object@rho.b)
            exps[object@ind[object@k]]
        }
        D_b <- Diagonal(x=tmp)
        
        ## Matrices we need for calculating the Jacobian
        DX <- D_e %*% X
        ZtD <- Zt %*% D_e
        XtDX <- crossprod(X, DX)
        ZtDX <- ZtD %*% X
        ZtDZ <- tcrossprod(ZtD, Zt)

        ## Initialize Jacobian Matrix (complete for given theta)
        J0 <- Matrix(0, p + q, p + q)
        J0[1:p, 1:p] <- XtDX
        
        ## CXt = C X\tr = solve(X\tr D.resp X, X\tr)
        CXt <- solve(XtDX, t(X))
        ## H = X CXt
        H <- X %*% CXt
        
        ## lfrac = la = \lambda_e / \lambda_b
        laD.re <- lfrac * D_b
        
        I <- Matrix(diag(n))
        ## P = I - D.resp H
        P <- I - D_e %*% H
        
        ZtPDZ <- tcrossprod(Zt %*% P, ZtD)
        CXtD <- CXt %*% D_e 
        
        ## Now we can assume that there are the matrices
        ## CXt, H, I, P, ZtPDZ, CXtD are in the environment

        ## get index of non-zero U
        idx <- !object@pp$zeroB
        
        if (any(idx)) {
            ## La = \Lambda_\theta
            Lat <- t(La <- object@pp$U_b)
            LtZt <- Lat %*% Zt
            
            ## ## Lai = \Lambda_\theta^{-1}
            ## ## lLtDL <- object@pp$lfrac * crossprod(Lai, object@pp$D_b %*% Lai)
            ## lLtDL <- object@pp$lfrac * solve(Lat, t(solve(Lat, object@pp$D_b[idx, idx])))
            
            ## Ms = M_\theta^* = solve(L\tr Z\tr D_resp P Z L + lD)
            ## Mst = Ms since D_resp is diagonal
            ## MsLtZt = solve(L\tr Z\tr P D_resp Z L + lD), L\tr Z\tr)
            ## LZMs = MsLtZt\tr
            Ms <- solve(Lat %*% ZtPDZ %*% La +
                        lfrac * D_b)
            MsLtZt <- Ms %*% LtZt
            
            ## Q = Q_\theta = CXt D.resp Z L Ms
            Qt <- tcrossprod(MsLtZt, CXtD)
            ## S = S_\theta = X Q
            St <- tcrossprod(Qt, X)
            ## T = Z L S\tr
            T <- crossprod(LtZt, St)

            ## K = K_\theta = Q\tr X\tr - Ms Z\tr
            K <- St - MsLtZt
            ## L = L_\theta = object@pp$lfrac * Ms
            L <- lfrac * Ms
            
            ## A = H - T - T\tr P + Z L Ms L\tr Z\tr
            A <- H - T - crossprod(T, P) + crossprod(MsLtZt, LtZt)
            ## B = lambda_e / lambda_b *(S - Z L Ms) 
            B <- lfrac * t(K)

            ## Complete Jacobian
            J <- J0
            J[p+(1:q), 1:p] <- t(J[1:p, p+(1:q)] <- crossprod(ZtDX, object@pp$U_b))
            t.idx <- (p+(1:q))[idx]
            J[t.idx, t.idx] <- crossprod(La[idx, idx], ZtDZ[idx, idx] %*% La[idx, idx]) + laD.re[idx, idx]
        } else {
            Ms <- solve(lfrac * D_b)
            A <- H
            ## Fixing Dimnames slot
            A@Dimnames <- list(A@Dimnames[[1]], NULL)
            B <- Matrix(0, object@pp$n, object@pp$q)
            K <- Matrix(0, object@pp$q, object@pp$n)
            L <- lfrac * Ms
            J <- J0
            if (method == "DASe2")
                expZero <- numeric(len(object, "theta"))
        }
        
        ## cat("    rho.resp:", rho.resp@name, "\n")
        ## cat("    rho.re:", rho.re@name, "\n")
        ## cat("   D_e[1,1]:", D_e[1,1], " vs. ", object@pp$D_e[1,1], "\n")
        ## cat("   D_b[1,1]:", D_b[1,1], " vs. ", object@pp$D_b[1,1], "\n")
        ## cat("      XtDX:", XtDX[1,1], " vs. ", object@pp$M_XX[1,1], "\n")
        ## if (exists("Ms")) cat("   Ms[1,1]:", Ms[1,1], " vs. ", object@pp$M()$M_bb[1,1], "\n")
        ## cat("    K[1,1]:", K[1,1], " vs. ", object@pp$K()[1,1], "\n")
        ## cat("    L[1,1]:", L[1,1], " vs. ", object@pp$L[1,1], "\n")
        ## cat("    A[1,1]:", A[1,1], " vs. ", object@pp$A[1,1], "\n")
        ## cat("    B[1,1]:", B[1,1], " vs. ", object@pp$B()[1,1], "\n")

        return(list(A = A, B = B, K = K, L = L, J = J,
                    D_e = D_e, D_b = D_b, Lambda_b = Diagonal(x=rep(lfrac, q)),
                    laD.re = laD.re))
    }, envir)
}

calcMatrices <- function(object, numpoints=13, envir=new.env()) {
    evalq({
        object <- get("object", envir = parent.frame(2))
        numpoints <- get("numpoints", envir = parent.frame(2))
        
        X <- object@pp$X
        Zt <- object@pp$Zt
        rho.resp <- object@rho.e
        rho.re <- object@rho.b
        method <- object@method

        p <- object@pp$p
        q <- object@pp$q
        n <- object@pp$n

        ## D_e and D_b matrices
        D_e <- Diagonal(x=rep(object@rho.e@EDpsi(), n))
        tmp <- if (object@wExp.b == 0) rep(object@rho.b@EDpsi(),q) else {
            exps <- sapply(object@dim, .calcE.D.re, rho = object@rho.b)
            exps[object@ind[object@k]]
        }
        D_b <- Diagonal(x=tmp)
        Lambda_b <- solve(D_b) * object@rho.e@EDpsi()

        ## Matrices we need for calculating the Jacobian
        DX <- D_e %*% X
        ZtD <- Zt %*% D_e
        XtDX <- crossprod(X, DX)
        ZtDX <- ZtD %*% X
        ZtDZ <- tcrossprod(ZtD, Zt)

        ## Initialize Jacobian Matrix (complete for given theta)
        J0 <- Matrix(0, p + q, p + q)
        J0[1:p, 1:p] <- XtDX
        
        ## CXt = C X\tr = solve(X\tr D.resp X, X\tr)
        CXt <- solve(XtDX, t(X))
        ## H = X CXt
        H <- X %*% CXt
        
        ## lfrac = la = \lambda_e / \lambda_b
        laD.re <- Lambda_b %*% D_b
        
        I <- Matrix(diag(n))
        ## P = I - D.resp H
        P <- I - D_e %*% H
        
        ZtPDZ <- tcrossprod(Zt %*% P, ZtD)
        CXtD <- CXt %*% D_e 
        
        ## Now we can assume that there are the matrices
        ## CXt, H, I, P, ZtPDZ, CXtD are in the environment

        ## get index of non-zero U
        idx <- !object@pp$zeroB
        
        if (any(idx)) {
            ## La = \Lambda_\theta
            Lat <- t(La <- object@pp$U_b)
            LtZt <- Lat %*% Zt
            
            ## Ms = M_\theta^* = solve(L\tr Z\tr D_resp P Z L + lD)
            ## Mst = Ms since D_resp is diagonal
            ## MsLtZt = solve(L\tr Z\tr P D_resp Z L + lD), L\tr Z\tr)
            ## LZMs = MsLtZt\tr
            Ms <- solve(Lat %*% ZtPDZ %*% La +
                        laD.re)
            MsLtZt <- Ms %*% LtZt
            
            ## Q = Q_\theta = CXt D.resp Z L Ms
            Qt <- tcrossprod(MsLtZt, CXtD)
            ## S = S_\theta = X Q
            St <- tcrossprod(Qt, X)
            ## T = Z L S\tr
            T <- crossprod(LtZt, St)

            ## K = K_\theta = Q\tr X\tr - Ms Z\tr
            K <- St - MsLtZt
            ## L = L_\theta = object@pp$lfrac * Ms
            L <- Ms %*% Lambda_b
            
            ## A = H - T - T\tr P + Z L Ms L\tr Z\tr
            A <- H - T - crossprod(T, P) + crossprod(MsLtZt, LtZt)
            ## B = lambda_e / lambda_b *(S - Z L Ms) 
            B <- t(K) %*% Lambda_b

            ## Complete Jacobian
            J <- J0
            J[p+(1:q), 1:p] <- t(J[1:p, p+(1:q)] <- crossprod(ZtDX, object@pp$U_b))
            t.idx <- (p+(1:q))[idx]
            J[t.idx, t.idx] <- crossprod(La[idx, idx], ZtDZ[idx, idx] %*% La[idx, idx]) + laD.re[idx, idx]
        } else {
            Ms <- solve(laD.re)
            A <- H
            ## Fixing Dimnames slot
            A@Dimnames <- list(A@Dimnames[[1]], NULL)
            B <- Matrix(0, object@pp$n, object@pp$q)
            K <- Matrix(0, object@pp$q, object@pp$n)
            L <- Ms %*% Lambda_b
            J <- J0
            if (method == "DASe2")
                expZero <- numeric(len(object, "theta"))
        }

        return(list(A = A, B = B, K = K, L = L, J = J,
                    D_e = D_e, D_b = D_b, Lambda_b = Lambda_b,
                    laD.re = laD.re))
    }, envir)
}


cmp <- function(rfm) {
    las <- calcMatrices(rfm)
    if (Matrix:::isDiagonal(rfm@pp$U_b)) {
        las2 <- calcMatricesDiagOnly(rfm)
        res <- all.equal(las, las2)
        if (!isTRUE(res)) {
            print(res)
            stop("Diagonal Test failed")
        }
    }
    res <- c(all.equal(las$D_e, rfm@pp$D_e),
             all.equal(las$D_b, rfm@pp$D_b),
             all.equal(las$Lambda_b, rfm@pp$Lambda_b),
             all.equal(las$A, rfm@pp$A),
             all.equal(las$B, rfm@pp$B()),
             all.equal(las$K, rfm@pp$K()),
             all.equal(las$L, rfm@pp$L),
             all.equal(las$J, rfm@pp$J()))
    if (!all(sapply(res, isTRUE))) {
        print(res)
        stop("Test failed")
    } else {
        cat("passed\n")
    }
}

test <- function(formula, data, ...) {
    rfm <- rlmer(formula, data, method="DASexp", ..., doFit=FALSE)
    rfm@pp$updateMatrices()
    theta0 <- theta(rfm)
    len <- length(theta0)
    ## test for original theta
    cat("Testing original theta... ")
    cmp(rfm)
    ## test for zero theta (all zero)
    cat("Testing all zero theta... ")
    setTheta(rfm, rep(0, len), fit.effects=FALSE)
    cmp(rfm)
    ## test one theta zero sequentially
    if (len > 1) {
        for (l in 1:len) {
            cat("Testing theta_", l, "=0... ", sep="")
            lth <- theta0
            lth[l] <- 0
            setTheta(rfm, lth, fit.effects=FALSE)
            cmp(rfm)
        }
    }
}

ttest <- function(formula, data) {
    test(formula, data)
    test(formula, data, rho.e = smoothPsi)
    test(formula, data, rho.b = smoothPsi)
    test(formula, data, rho.e = smoothPsi, rho.b = chgDefaults(smoothPsi, k=1, s=10))
}

cat("----- Dyestuff -----\n")
ttest(Yield ~ (1|Batch), Dyestuff)
cat("----- Sleepstudy -----\n")
ttest(Reaction ~ Days + (Days|Subject), sleepstudy)
cat("----- Sleepstudy2 -----\n")
sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
ttest(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)


### test update bug
zero.theta.DASexp <- robustlmm:::zero.theta.DASexp
rfm <- rlmer(Yield ~ (1|Batch), Dyestuff, wExp.e = 0, wExp.b = 0)
rho0 <- smoothPsi
rho <- chgDefaults(rho0, k = 1, s = rho0@tDefs[2])

rfm1 <- update(rfm, method="DASexp",
               method.effects = "IRWLS",
               rho.e = rho, rho.b = rho,
               wExp.e = 2, wExp.b = 2)
theta0 <- 3 
val1 <- zero.theta.DASexp(theta0, rfm1, 0)
val2 <- zero.theta.DASexp(theta0, rfm1, 0)
stopifnot(all.equal(val1, val2, 1e-6))

## test refClass copy problem
rfm2 <- update(rfm, method="DASexp",
               method.effects = "IRWLS",
               rho.e = rho, rho.b = rho,
               wExp.e = 2, wExp.b = 2)
val5 <- zero.theta.DASexp(theta0, rfm2, 0)
stopifnot(all.equal(val1, val5))

showPP <- function(field) {
    print(head(rfm@pp[[field]]))
    print(head(rfm1@pp[[field]]))
    print(head(rfm2@pp[[field]]))
}

aeq <- function(fields = names(getRefClass(class(rfm1@pp))$fields())) {
    res <- character(0)
    for (field in fields) {
        res <- c(res, all.equal(rfm1@pp[[field]], rfm2@pp[[field]]))
    }
    names(res) <- fields
    nt <- res != "TRUE"
    if (any(nt)) res[nt] else TRUE
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
