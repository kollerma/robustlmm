fitEffects <- function(par, object, rel.tol = 1e-8, max.iter = 1000) {
  n <- object@pp$n
  p <- object@pp$p
  q <- object@pp$q
  idx <- c(rep(TRUE, p), !.zeroB(object))
  par[!idx] <- 0 ## set random effects corresponding to dropped vc to 0
  setU(object, par[p+1:q])
  setFixef(object, par[1:p])
  if (inherits(object@pp, "rlmerPredD")) {
    ## cat("idx: ", idx, "\n")
    ## cat("par:", par, "\n")
    LtZt <- .U_btZt.U_et(object)
    MAT1 <- Matrix(0, n+q, p+q)
    MAT1[1:n,1:p] <- ..U_eX(object)
    MAT1[1:n,p+1:q] <- t(LtZt)
    MAT1[n+1:q,p+1:q] <- sqrt(.Lambda_b(object))
    MAT2 <- Matrix(0, p+q, n)
    MAT2[1:p,] <- t(..U_eX(object))
    MAT2[p+1:q,] <- LtZt
    .U_ey <- solve(.U_e(object), object@resp$y)
    converged <- FALSE
    iter <- 1
    while(!converged && iter <= max.iter) {
      ## set parameters
      ## calculate weights
      sqW <- Diagonal(x=sqrt(c(w.e <- wgt.e(object), wgt.b(object))))
      ## build matrix (a factor, resp)
      MAT <- sqW %*% MAT1
      w.y <- w.e * .U_ey
      ## solve to get new parameters
      par1 <- drop(solve(crossprod(MAT), MAT2 %*% w.y))
      ## check convergence
      converged <- sum(abs(par - par1)) < rel.tol * sum(abs(par))
      ## cat("IRWLS Iter ", iter, ", par:", par1, ", conv:",
      ##     sprintf("%.12f < %.12f", sum(abs(par - par1)),
      ##             rel.tol * sum(abs(par))), "\n")
      ## update parameters
      setU(object, par1[p+1:q])
      setFixef(object, par1[1:p])
      par <- par1
      iter <- iter + 1
    }
  } else {
    FitEffects <- new(FitEffects, object@pp$ptr()$.pointer,
                      object@resp$ptr()$.pointer, rel.tol, max.iter);
    FitEffects$fit()
  }
  object
}
