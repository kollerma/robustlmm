zero.theta.DASexp <- function(par0, object, VERBOSE, ...) {
    ## take care of values outside boundary
    par <- par0
    par[idx <- par0 < object@lower] <- 0
    ## continue with sanitized parameters
    ## compare with updateTheta
    setTheta(object, par, fit.effects = TRUE, ...)
    skbs <- .sumKappa(object)
    
    ret <- rep(0, length(par))
    sigmae <- object@pp$sigma
    ## FIXME: this is not correct in the non-diagonal case
    db <- dist.b(object, sigmae)
    idxTheta <- 0
    for (block in seq_along(object@blocks)) {
        if (VERBOSE > 3) {
            cat("Block:", block, "\n")
            if (VERBOSE > 4) cat(sprintf("skbs[[%i]]: %.12f\n", block, skbs[[block]]))
        }
        ldim <- object@dim[block]
        lq <- object@q[block]
        lind <- object@ind[object@k] == block ## FIXME: == as.vector(object@idx[[blok]])
        wgtds <- object@rho.sigma.b[[block]]@wgt(db[lind])
        wgtb.s <- sqrt(wgtds) * object@pp$b.s[lind]
        idxTheta <- max(idxTheta) + 1:(ldim*(ldim+1)/2)
        wgtb.ss <- wgtb.s / sigmae
        if (all(abs(object@pp$b.s[lind] / sigmae) < 1e-7)) {
            if (VERBOSE > 3) cat("all(abs(b.ss)) < 1e-7 is TRUE, setting theta to 0\n")
            ret[idxTheta] <- 0
        } else {
            if (ldim > 1) {
                ## non-diag case
                tmp <- tcrossprod(Matrix(wgtb.ss, ldim)) - skbs[[block]]
                if (VERBOSE > 3) str(tmp)
                ret[idxTheta] <- tmp[lower.tri(tmp, TRUE)]
            } else {
                ## simple 1d case
                ret[idxTheta] <- sum(wgtb.ss^2) - drop(skbs[[block]])
            }
            if (VERBOSE > 3) ret[idxTheta]
        }
    }
    ## set the values, such that there is a root at zero
    ret[idx] <- -1*par0[idx]

    if (VERBOSE > 3) {
        cat(sprintf("theta:   %.12f\n", theta(object)))
        cat("coef:   ", object@pp$beta,"\n")
        cat("b.s:    ", b.s(object), "\n")
        cat("sigmae: ", object@pp$sigma, "\n")
        if (VERBOSE > 4) {
            cat("idx:    ", idx, "\n")
            cat("dk:     ", .dk(object, sigmae), "\n")
        }
    }
    
    if (VERBOSE > 2)
        cat("zero.theta(", par0, ") =", ret, "\n")
    ret
}
