######################################################
## Functions to generate commonly required objects   ##
#######################################################

##' std.b: Return the spherical random effects or
##'   "Standardize" the Matrix matrix: \eqn{\Lambda^{-1} matrix / \sigma}{Lambda^-1 matrix / sigma}
##'
##' @title Standardized values
##' @param object rlmerMod object
##' @param sigma to use for standardization
##' @param matrix matrix to standardize
##' @param drop apply drop to result?
##' @param t transpose result
##' @rdname std
std.b <- function(object, sigma = object@pp$sigma, matrix, drop=TRUE, t=FALSE) 
    object@pp$stdB(sigma, matrix, drop, t)

##' std.e: Calculate the standardized residuals or
##'   "Standardize" the Matrix sigma: \eqn{R^{-1} matrix / \sigma}{R^-1 matrix / sigma}
##'
##' @rdname std
std.e <- function(object, sigma = object@pp$sigma, matrix, drop=TRUE) {
    if (missing(matrix)) return(object@resp$wtres / sigma)
    ## for the moment: just divide by sigma
    if (drop) matrix <- drop(matrix)
    matrix/sigma
}

##' sumRho.b: Calculate sum of rho of random effects
##'
##' @title sumRho
##' @param object rlmerMod object
##' @param sigma to use for standardization
##' @param lambda apply normalization by lambda
##' @param ... ignored
##' @rdname sumRho
sumRho.b <- function(object, sigma = object@pp$sigma, lambda = FALSE, ...) {
    ret <- sum(object@rho.b@rho(dist.b(object, sigma))) ## dist instead of std
    if (lambda) ret/object@rho.b@EDpsi() else ret
}

##' sumRho.e: Calculate sum of rho of residuals
##'
##' @rdname sumRho
sumRho.e <- function(object, sigma = object@pp$sigma, lambda = FALSE, ...) {
    ret <- sum(object@rho.e@rho(dist.e(object, sigma))) ## dist instead of std
    if (lambda) ret/object@rho.e@EDpsi() else ret
}

##' Calculate scaled Mahalanobis distances per group
##'
##' @title Distances per group
##' @param object rlmerMod object
##' @param sigma to use for standardization
##' @param k group indicator vector
##' @param b.s spherical random effects
##' @param ... ignored
.dk <- function(object, sigma, k = object@k, bs = u(object), ...) {
    ## aggregate and rescale by size of group
    sqrt(aggregate(bs^2, list(k), mean)$x) / sigma ## mean instead of sum
    
}

##' dist.b: Calculate the distance from 0 standardized by sigma.
##'   This is just value divided by sigma for uncorrelated
##'   observations. For correlated items, this is the Mahalanobis
##'   distance from 0.
##'
##' @title Calculate distance
##' @param object object to use
##' @param sigma scale for standardization
##' @rdname dist
dist.b <- function(object, sigma = object@pp$sigma, ...) {
    k <- object@k
    .dk(object, sigma, k, ...)[k]
}

##' dist.e: Calculate dist for residuals
##'   always assume they are uncorrelated
##'
##' @rdname dist
dist.e <- function(object, sigma = object@pp$sigma) {
    std.e(object, sigma) ## just the usual rescaled residuals
}

##' wgt.b: Calculate the robustness weights psi(d) / d,
##'   standardized by sigma. The robustness weights are calculated
##'   with d the Mahalanobis distance. Each group of correlated items
##'   then gets a constant weight. 
##'
##' @title Calculate robustness weights
##' @param object object to use
##' @param sigma scale for standardization
##' @param use.c.sigma whether use c.sigma tuning parameters or the default ones.
##' @rdname wgt
##' @export
wgt.b <- function(object, sigma = object@pp$sigma, use.c.sigma = FALSE)
    if (use.c.sigma) object@rho.sigma.b@wgt(dist.b(object, sigma)) else
       object@rho.b@wgt(dist.b(object, sigma))

##' wgt.e: robustness weights of residuals
##'
##' @rdname wgt
##' @export
wgt.e <- function(object, sigma = object@pp$sigma, use.c.sigma = FALSE)
    if (use.c.sigma) object@rho.sigma.e@wgt(dist.e(object, sigma)) else
       object@rho.e@wgt(dist.e(object, sigma))

### Calculate robustness weights * squared effect
### Return sensible result in the infinite case
## wgt(x)^wExp * y^2
## uses M-estimator of scale robustness weights for wExp == 0.
## Assume: x infinite <=> y infinite
.wgtxy2 <- function(rho, x, y, wExp)
    .wgtxy(rho, x, y*y, wExp)
## wgt(x)^wExp * y
.wgtxy <- function(rho, x, y, wExp) {
    switch(wExp + 1,
       { ## wExp == 0
           ret <- rho@rho(x)/(x*x)*2 * y
           ret[is.infinite(x) & is.infinite(y)] <- rho@rho(Inf)^2
       },
       { ## wExp == 1
           ret <- rho@wgt(x) * y
           ret[is.infinite(x) & is.infinite(y)] <- rho@psi(Inf) ## * Inf
       },
       { ## wExp == 2
           ret <- rho@wgt(x)^2 * y
           ret[is.infinite(x) & is.infinite(y)] <- rho@psi(Inf)^2
       })
    ret
}

##' wgtSum.b: Calculate \eqn{\sum(\psi(b^*/\sigma))}{sum(psi(b^*/sigma)}
##'
##' @title wgtSum(b)
##' @param object rlmerMod object
##' @param sigma sigma to be used
##' @param lambda whether to normalize with EDpsi()
##' @param wExp wExp parameter
##' @param use.c.sigma use rho.b or rho.sigma.b
##' @rdname wgtSum
wgtSum.b <- function(object, sigma, lambda = FALSE, wExp = 0, use.c.sigma = FALSE) {
    rho <- if (use.c.sigma) object@rho.sigma.b else object@rho.b
    ret <- sum(.wgtxy2(rho, dist.b(object, sigma), std.b(object, sigma), wExp))
    if (lambda) ret/rho@EDpsi() else ret
}

##' wgtSum.e: Calculate \eqn{\sum(\psi(r/\sigma))}{sum(psi(r/sigma)}
##'
##' @rdname wgtSum
wgtSum.e <- function(object, sigma, lambda = FALSE, wExp = 0, use.c.sigma = FALSE) {
    rho <- if (use.c.sigma) object@rho.sigma.e else object@rho.e
    ret <- sum(.wgtxy2(rho, dist.e(object, sigma), std.e(object, sigma), wExp))
    if (lambda) ret/rho@EDpsi() else ret
}

##' detCov.b: Calculate determinant of D / 2
##'
##' @title determinant of covariance matrices
##' @param object rlmerMod object
##' @param logarithm return logarithm of determinant
##' @param sigma which value of sigma to use.
##' @rdname detCov
detCov.b <- function(object, logarithm = TRUE, sigma = object@pp$sigma) {
    ## need to import diag from Matrix or call explicitly:
    ##ldiag <- Matrix::diag(object@EDpsi)
    ldiag <- diag(object@pp$U_b)
    stopifnot(all(ldiag >= 0))
    ## get number of fitted random effects
    q <- sum(!object@pp$zeroB)
    ## ignore the 0 terms:
    ## these correspond to dropped covariance parameters
    ldiag <- ldiag[ldiag > 0]
    ## return 0 if all are 0
    ldet <- if (q == 0) return(0) else 2*(q*log(sigma) + sum(log(ldiag)))
    if (logarithm) ldet else exp(ldet)
}

##' detCov.e: Calculate determinant of R / 2
##' @rdname detCov
detCov.e <- function(object, logarithm = TRUE, sigma = object@pp$sigma) {
    ## along mkDet of Matrix
    ## assume simple diagonal covariance matrix
    ldet <- 2* object@pp$n * log(sigma)
    if (logarithm) ldet else exp(ldet)
}

##' detCov.rest: Calculate determinant of
##' \eqn{Z^T V_e^{-1} Z + V_b}{Z^T V_e^{-1} Z + V_b}
##' ML case only
##' @rdname detCov
detCov.rest <- function(object, logarithm = TRUE, sigma = object@pp$sigma) {
    idx <- !object@pp$zeroB
    mat <- crossprod(object@pp$.U_eZ[, idx]) + crossprod(solve(object@pp$U_b[idx, idx]))
    ldet <- determinant(mat, logarithm=TRUE)$modulus[[1]] - 2 * sum(idx) * log(sigma)
    if (logarithm) ldet else exp(ldet)
}

##' @method robll rlmerMod
##' @rdname robll
##' @param add.q if TRUE, result is incremented by q
##' @S3method robll rlmerMod
robll.rlmerMod <- function(object, norm.only = FALSE, add.q = TRUE, ...) {
    ## calculate norm parts
    ret <- sumRho.e(object, object@pp$sigma, lambda=TRUE, ...) +
        sumRho.b(object, object@pp$sigma, lambda=TRUE, ...)
    if (norm.only) return(ret)
    ## add constants and determinants
    log(2*pi)*(object@pp$n +
               if (add.q) sum(!object@pp$zeroB) else 0) +
                   detCov.e(object) + detCov.b(object) + 2*ret
}

##' @method Q rlmerMod
##' @rdname Q
##' @S3method Q rlmerMod
Q.rlmerMod <- function(object, determinant = TRUE, log = determinant, ...) {
    idx <- !object@pp$zeroB
    q <- sum(idx)
    update.cache <- FALSE
    if (is.null(object@cache$Q.D.re) || any(idx != object@cache$Q.idx)) {
        update.cache <- TRUE
        object@cache$Q.idx <- idx
        ##cat("updating Q cache\n")
        
        Z <- t(object@pp$Zt)
        Z <- as(Z, "dgeMatrix")
        
        ## be careful about 0 variance components (without correlation)
        ## drop the corresponding columns in RicZ and Dic
        D.re <- object@pp$D_b
        if (! all(idx)) {
            Z <- Z[,idx]
            D.re <- D.re[idx,idx]
        }
        object@cache$Q.D.resp <- object@pp$D_e / object@rho.e@EDpsi()
        if (q > 0)
            object@cache$Q.0part <- crossprod(Z, object@cache$Q.D.resp %*% Z)
        object@cache$Q.D.re <- D.re / object@rho.b@EDpsi()
        object@cache$Q.Z <- Z
    } else  {
        ##cat("not updating...\n")
        D.re <- object@cache$Q.D.re
    }
    
    Dic <- std.b(object, 1, Diagonal(object@pp$q))
    Dic <- as(Dic, "dgeMatrix")
    
    if (! all(idx)) Dic <- Dic[idx,idx]

    if (q > 0)
        Q.0 <- object@cache$Q.0part + crossprod(Dic, D.re %*% Dic)

    if (determinant && !.isREML(object)) {
        ## ML case
        if (q > 0) {
            ret <- determinant(Q.0, log)$modulus
            if (log)
                return(ret - 2*q*log(object@pp$sigma))
            else
                return(ret/object@pp$sigma^(2*q))
        } else return(0) # no VC > 0
    }

    ## REML case (or determinant not requested)
    if (update.cache || is.null(object@cache$Q.1)) {
        X <- object@pp$X
        object@cache$Q.1 <- crossprod(X, object@cache$Q.D.resp %*% X)
        ## only need to calculate Q.3 if there are any nnz vc
        if (q > 0) {
            object@cache$Q.3 <- crossprod(X, object@cache$Q.D.resp %*% object@cache$Q.Z)
            object@cache$Q.partDet <- crossprod(object@cache$Q.3,
                                             solve(object@cache$Q.1, object@cache$Q.3))
        }
    }
    p <- nrow(object@cache$Q.1)

    if (determinant) {
        if (update.cache || is.null(object@cache$Q.detQ.1raw)) {
            object@cache$Q.detQ.1raw <- ret <- determinant(object@cache$Q.1)$modulus
        } else ret <- object@cache$Q.detQ.1raw
        ## only need to calculate rest if there are any nnz vc
        if (q > 0) {
            tmat <- Q.0 - object@cache$Q.partDet
            tval <- determinant(tmat)$modulus
            if (is.infinite(tval)) {
                ## use alternative calculation method in the near singular case
                tval <- abs(prod(diag(qr.R(qr(tmat)))))
                warning("Q, caught infinite value, using QR:", tval)
                tval <- log(tval)
            }
            ret <- ret + tval
        }
        
        ret <- ret - 2*(p+q)*log(object@pp$sigma)
        if (!log) {
            ret <- exp(ret) ## ignore sign, always positive???
            attr(ret, "logarithm") <- FALSE
        }
        ##cat("    Q:", ret, "\n")
        return(ret)
    }

    ## FIXME: matrix is of wrong size!!??
    ##        should it always be of full size??

    if (update.cache || is.null(object@cache$Q.raw)) {
        ## cannot use Matrix here: otherwise this fails
        ## af ... <- t(Q.3) (first execution only)
        Q <- matrix(0, p+q, p+q)
        Q[1:p, 1:p] <- as(object@cache$Q.1, "matrix")
        if (q > 0) {
            Q[1:p, p+1:q] <- Q.3 <- as(object@cache$Q.3, "matrix")
            Q[p+1:q, 1:p] <- t(Q.3)
        }
        object@cache$Q.raw <- Q
    } else Q <- object@cache$Q.raw
    
    if (q > 0) Q[p+1:q, p+1:q] <- as(Q.0, "matrix")
    
    as(Q/object@pp$sigma^2, "Matrix")
}

### Calculate the Laplace Approximation
### assume wExp = 0
## laplace.rlmerMod <- function(object, ...) {
##     ret <- robll(object, FALSE, add.q = FALSE)
##     ## cat("Laplace: ", ret,
##     ##     if (isREML(object)) log(2*pi)*object@pp$p) else 0,
##     ##     Q(object, TRUE, TRUE, wExp = wExp), "\n")
##     if (.isREML(object)) ret <- ret - log(2*pi)*object@pp$p
##     if (object@use.laplace)
##         ret <- ret + Q(object, TRUE, TRUE)
##     attributes(ret) <- NULL
##     ret
## }

##' Update the value of the deviance in the rlmerMod object
##'
##' @title Calculate the deviance
##' @param object rlmerMod object
##' @param ... ignored
updateDeviance <- function(object, ...) {
    if (object@use.laplace || !.isREML(object)) {
        ret <- robll(object, FALSE, add.q = FALSE)
        if (.isREML(object)) ret <- ret - log(2*pi)*object@pp$p
        ret <- ret + if (object@use.laplace) Q(object, TRUE, TRUE) else detCov.rest(object, TRUE)
    } else {
        ## REML non-laplace case
        stop("not implemented: deviance without using Laplace approx and REML = TRUE")
        ## FIXME: RZX depends on theta and is not updated
        ## idx <- !object@pp$zeroB
        ## V <- tcrossprod(object@pp$.U_eZ[,idx] %*% object@pp$U_b[idx,idx]) + Diagonal(object@pp$n)
        ## ret <- (object@pp$n-object@pp$p)*(log(2*pi) + 2*log(object@pp$sigma)) +
        ##     determinant(crossprod(object@pp$X) - crossprod(object@pp$RZX), TRUE)$modulus[[1]] +
        ##         determinant(V, TRUE)$modulus[[1]] + 2*robll(object, TRUE, add.q = FALSE)
    }
    attributes(ret) <- NULL
    object@pp$setDeviance(ret)
    invisible(ret)
}
##' @rdname updateDeviance
##' @method laplace rlmerMod
##' @S3method laplace rlmerMod
laplace.rlmerMod <- updateDeviance

##' @rdname gradient
##' @method gradient rlmerMod
##' @S3method gradient rlmerMod
gradient.rlmerMod <- function(object, ...) {
    rho.resp <- object@rho.e
    rho.re <- object@rho.b
    
    ## correct by rho@EDpsi() as in robll
    Ricpsir <- std.e(object, object@pp$sigma,
                     Matrix(rho.resp@psi(std.e(object, object@pp$sigma))/
                            rho.resp@EDpsi()))
    grad.u <- drop(crossprod(object@pp$U_b, (object@pp$Zt %*% Ricpsir))) - 
        ## rho.re@psi(object@pp$b.s/object@pp$sigma)/rho.re@EDpsi()/object@pp$sigma
        wgt.b(object)  * std.b(object, object@pp$sigma) / object@pp$sigma / rho.re@EDpsi()
    
    ## setting gradient to 0 for components corresponding to
    ## variance components equal to 0
    if (any(idx <- object@pp$zeroB))
        grad.u[idx] <- 0

    ## FIXME: this is probably wrong for nontrivial U_e.
    -1*c(drop(crossprod(object@pp$X, Ricpsir)), grad.u)
}

##' Find blocks of correlated random effects
##'
##' @title Find blocks in Lambda
##' @param obj reTrms-object
##' @return list of blocks
findBlocks <- function(obj) {
    obj <- as(obj, "merPredD")
    LambdaInd <- obj$Lambdat
    LambdaInd@x[] <- as.double(obj$Lind)
    LambdaInd <- t(LambdaInd)
    LambdaInd <- as(LambdaInd, "matrix") ## to avoid attempt to apply non function error
    bend <- unique(apply(LambdaInd != 0, 2, function(x) max(which(x))))
    nblocks <- length(bend)
    bstart <- c(1, bend[-nblocks]+1)
    bidx <- lapply(1:nblocks, function(i) seq.int(bstart[i],bend[i]))
    blocks <- lapply(bidx, function(idx) LambdaInd[idx,idx])
    bind <- match(blocks, ublocks <- unique(blocks))
    k <- unlist(lapply(1:nblocks, function(i) rep(i, length(bidx[[i]]))))
    bdim <- sapply(ublocks, NCOL)
    bidx <- lapply(1:length(ublocks), function (i) 
                   matrix(unlist(bidx[bind == i]),nrow = bdim[i]))
    q <- sapply(bidx, length)
    list(blocks = ublocks, ind = bind, idx = bidx, dim = bdim, q = q, k = k)
}

#######################################################
## Summary / printing methods                        ##
#######################################################

### Print method (calls print(lmerMod) adds some info about rho)
.printRlmerMod <- function(x, ...) {
    s <- lme4:::summary.merMod(x)
    out <- capture.output(print(s, ...))
    ## check whether the methods slot is there
    ## (for compatibility with older versions of rlmer)
    if (.hasSlot(x, "method"))
        out[1] <- sprintf("Linear mixed model fit by %s ['%s']",
                          x@method, class(x))
    cat(out, sep="\n")
    cat("\nRho functions used for fitting:\n")
    cat("   Residuals: ", robustbase:::.sprintPsiFunc(x@rho.e), "\n")
    cat("      w exp.: ", x@wExp.e, "\n")
    if (!isTRUE(all.equal(x@rho.e@tDefs, x@rho.sigma.e@tDefs)))
        cat("     c.sigma: ", paste(names(x@rho.sigma.e@tDefs), round(x@rho.sigma.e@tDefs, 3),
                                    sep=" = ", collapse=", "), "\n")
    
    cat(" Random Eff.: ", robustbase:::.sprintPsiFunc(x@rho.b), "\n")
    cat("      w exp.: ", x@wExp.b, "\n")
    if (!isTRUE(all.equal(x@rho.b@tDefs, x@rho.sigma.b@tDefs)))
        cat("     c.sigma: ", paste(names(x@rho.sigma.b@tDefs), round(x@rho.sigma.b@tDefs, 3),
                                    sep=" = ", collapse=", "), "\n")
    invisible(x)
}

##' @S3method print rlmerMod
print.rlmerMod <- .printRlmerMod
##' @exportMethod show
setMethod("show", "rlmerMod", function(object) .printRlmerMod(object))

##' @rdname getInfo
##' @method getInfo lmerMod
##' @S3method getInfo lmerMod
getInfo.lmerMod <- function(object, ...) {
    lsum <- summary(object)
    .namedVector <- function(mat) {
        if (is.vector(mat)) return(mat)
        names <- rownames(mat)
        ret <- drop(mat)
        names(ret) <- names
        ret
    }
    .getVC <- function(varcor) {
        vc <- lapply(varcor, function(grp) attr(grp, "stddev"))
        ret <- unlist(vc, use.names = FALSE)
        names(ret) <-
            unlist(lapply(1:length(vc), function(i)
                          paste(names(vc)[i], names(vc[[i]]))))
        ##.namedVector(ret)
        ret
    }
    .getCorr <- function(varcor) {
        ret <- lapply(1:length(varcor), function(i) {
            grp <- varcor[[i]]
            corr <- attr(grp, "correlation")
            names <- outer(paste(names(varcor)[i], colnames(corr)),
                           paste("x", rownames(corr)), paste)
            ret <- as.vector(corr[lower.tri(corr)])
            names(ret) <- as.vector(names[lower.tri(names)])
            ret
        })
        unlist(ret)
    }
    ret <- list(data = object@call$data,
                coef = .namedVector(lsum$coefficients[,1,drop=FALSE]),
                varcomp = .getVC(lsum$varcor),
                sigma = lsum$sigma)
    corrs <- .getCorr(lsum$varcor)
    if (length(corrs) > 0) ret$correlations <- corrs
    if (.isREML(object)) {
        ret$REML <- object@devcomp$cmp['REML']
    } else {
        ret$deviance <- object@devcomp$cmp['dev']
    }
    ret
}

##' @rdname getInfo
##' @method getInfo rlmerMod
##' @S3method getInfo rlmerMod
getInfo.rlmerMod <- function(object, ...) {
    linfo <- getInfo(as(object, "lmerMod"))
    linfo$rho.e <- robustbase:::.sprintPsiFunc(object@rho.e, TRUE)
    linfo$wExp.e <- object@wExp.e
    linfo$rho.b <- robustbase:::.sprintPsiFunc(object@rho.b, TRUE)
    linfo$wExp.b <- object@wExp.b
    linfo
}

##' Compare the fits of multiple lmerMod or rlmerMod objects
##'
##' @title Create a comparison chart for multiple fits
##' @param ... objects to compare
##' @param digits number of digits to show in print
##' @param dnames names of objects given as arguments (optional)
##' @return comparison table
##' @export
compare <- function(..., digits = 3, dnames = NULL) {
    linfos <- list(...)
    if (!missing(dnames) && !is.null(dnames)) names(linfos) <- dnames
    linfos <- lapply(linfos, getInfo)
    ## check if all methods work at least on the same dataset
    if (length(unique(sapply(linfos, function(x) x$data))) > 1)
        stop("Comparison has to be for objects fitted to the same dataset")
    ## local helper functions
    .NULLtoNA <- function(lst)
        lapply(lst, function(x) if (is.null(x)) NA else x)
    .getComp <- function(lst, comp) 
        sapply(lst, function(x) .NULLtoNA(x[comp]))
    .getComp2 <- function(lst, comp)
        sapply(lst, function(x) x[comp])
    .dropComp <- function(lst, comp)
        lapply(lst, function(x) { x[comp] <- NULL; x })
    .getNames <- function(lst) 
        unique(unlist(lapply(lst, names)))
    .combineComp <- function(lst) {
        names <- .getNames(lst)
        tmp <- .getComp2(lst, names)
        if (!is.matrix(tmp))
            tmp <- matrix(tmp, ncol = length(tmp))
        rownames(tmp) <- names
        format(tmp, digits = digits)
    }
    lnames <- .getNames(linfos)
    ## header
    call <- match.call()
    call$digits <- NULL
    call$names <- NULL
    split <- rep("", length(linfos))
    ret <- rbind(Coef=split, 
                 .combineComp(.getComp(linfos, "coef")))
    ret <- rbind(ret, NULL=split, VarComp=split, 
                 .combineComp(.getComp(linfos, "varcomp")))
    if (any(sapply(linfos, function(linfo) !is.null(linfo$correlations))))
        ret <- rbind(ret, NULL=split, Correlations=split,
                     .combineComp(.getComp(linfos, "correlations")))
    linfos <- .dropComp(linfos, c("data", "coef", "varcomp", "correlations"))
    ret <- rbind(ret, NULL=split)
    for (name in .getNames(linfos)) {
        ret <- rbind(ret, format(.getComp(linfos, name), digits = digits))
        rownames(ret)[nrow(ret)] <- name
    }
    ret <- gsub("\\s*(NA|NULL)", "", ret)
    colnames(ret) <- if (missing(dnames) && is.null(dnames)) as.character(call)[-1] else dnames
    rownames(ret)[rownames(ret) == "NULL"] <- ""
    ret   
}
 
##' @S3method update rlmerMod
update.rlmerMod <- function(object, formula., ..., evaluate = TRUE) {
    ## update call
    ## set old object as init object
    ## run it
    if (is.null(call <- object@call))
        stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.))
        call$formula <- update.formula(formula(object), formula.)
    ## set init to object, if not explicitly given (and no new data given)
    if (is.null(extras[["data"]])) {
        if (is.null(extras[["init"]])) {
            lcall <- sys.call(sys.parent())
            extras$init <- object
        }
        ## copy pp and resp (to really get a new object)
        extras$init@pp <- object@pp$copy()
        ## reset calledInit... fields to FALSE:
        fields <- grep("calledInit", names(getRefClass(class(extras$init@pp))$fields()), value=TRUE)
        Map(function(field) extras$init@pp$field(field, FALSE), fields)
        extras$init@resp <- object@resp$copy()
    } else {
        extras$init <- NULL
    }
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}
