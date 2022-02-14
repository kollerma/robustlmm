require(robustlmm)
require(rgl)

.d <- function(bs, s=length(bs)) {
    if (s == 1) return(bs)
    if (is.matrix(bs)) rowSums(bs*bs) else sum(bs*bs)
}
colorRampRgb <- function(vals, colors = c("blue", "red")) {
  vals <- vals
  mat <- colorRamp(colors)(vals)
  rgb(mat[, 1], mat[, 2], mat[, 3], (1 - vals) * 255, maxColorValue = 255)
}
plotArrow <- function(i, from, to) {
  p0 <- from[i, ]
  p1 <- p0 - to[, i]
  if (any(is.na(p1))) {
    return()
  }
  arrow3d(p0, p1, type = "lines")
}

Lkk <- matrix(c(0.555052111908150, -0.620130518640817, -0.620130518640816,  0.908619143206113), 2, 2)
Sk <- matrix(c(0.24507904474548758, 0.0, 0.00311004948666043, 0.24417103302921830), 2, 2)
skappa <- 1.57223325559661
s <- 2
gh <- robustlmm:::ghq(13, FALSE)
gh$nodes <- gh$nodes * sqrt(2)
gh$weights <- gh$weights / sqrt(pi)
ghZ <- as.matrix(expand.grid(gh$nodes, gh$nodes, gh$nodes, gh$nodes))
ghw <- apply(as.matrix(expand.grid(gh$weights, gh$weights, gh$weights, gh$weights)), 1, prod)

wgt <- chgDefaults(smoothPsi, k = 1.2, s = 10)@wgt
wgt.sigma <- chgDefaults(smoothPsi, k = 1.5, s = 10)@wgt
psi.sigma <- chgDefaults(smoothPsi, k = 1.5, s = 10)@psi
wgtDelta <- function(u) (psi.sigma(u) - psi.sigma(u-skappa))/s

lTbk <- matrix(c(0.6532058, 0.3874415, 0.3874415, 0.4323062), 2, 2)

optGH.wrong <- function(TkbLTri) {
  lTbk <- matrix(c(TkbLTri[1:2], TkbLTri[2:3]), 2, 2)
  qrlTbk <- try(qr(lTbk), silent=TRUE)
  if (is(qrlTbk, "try-error") || qrlTbk$rank < s) {
    return(c(NA, NA, NA))
  }
  btilde <- ghZ[,1:2] - wgt(.d(ghZ[,1:2],2)) * ghZ[, 1:2] %*% Lkk -
    ghZ[, 3:4] %*% Sk
  tmp1 <- colSums(qr.solve(qrlTbk, t(btilde))^2)
  tmp2 <- btilde[,1] * btilde[,2]
  a <- sum(wgtDelta(tmp1) * ghw)
  B <- matrix(colSums(wgt.sigma(tmp1) * ghw *
                        matrix(c(btilde[,1]*btilde[,1], tmp2, tmp2,
                                 btilde[,2]*btilde[,2]), length(ghw))),2)
  lTbk1 <- B/a
  ## return(lTbk1)
  (lTbk - lTbk1)[-3]
}

optGH <- function(TkbLTri) {
  lTbk <- matrix(c(TkbLTri[1:2], TkbLTri[2:3]), 2, 2)
  lLTbk <- try(chol(lTbk), silent = TRUE)
  if (is(lLTbk, "try-error")) {
    return(c(NA, NA, NA))
  }
  btilde <- ghZ[,1:2] - wgt(.d(ghZ[,1:2],2)) * ghZ[, 1:2] %*% Lkk -
    ghZ[, 3:4] %*% Sk
  tmp1 <- colSums(backsolve(lLTbk, t(btilde))^2)
  tmp2 <- btilde[,1] * btilde[,2]
  a <- sum(wgtDelta(tmp1) * ghw)
  B <- matrix(colSums(wgt.sigma(tmp1) * ghw *
                        matrix(c(btilde[,1]*btilde[,1], tmp2, tmp2,
                                 btilde[,2]*btilde[,2]), length(ghw))),2)
  lTbk1 <- B/a
  ## return(lTbk1)
  (lTbk - lTbk1)[-3]
}

optCub <- function(TkbLTri) {
  require(cubature)
  lTbk <- matrix(c(TkbLTri[1:2], TkbLTri[2:3]), 2, 2)
  lLTbk <- try(chol(lTbk), silent = TRUE)
  if (is(lLTbk, "try-error")) {
    return(c(NA, NA, NA))
  }
  funA <- function(u) {
    btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] - crossprod(Sk, u[3:4])
    wgtDelta(drop(crossprod(backsolve(lLTbk, btilde))))*prod(dnorm(u))
  }
  aCub <- robustlmm:::test_Hcubature_ninfInf2(funA, 4, 1, 5e-2)
  ## aCub <- robustlmm:::test_GaussianQuadratureNd_ninfInf(funA, 4, 1, 13)
  funB <- function(u) {
    btilde <- u[1:2] - wgt(.d(u[1:2],2)) * Lkk %*% u[1:2] -
      crossprod(Sk, u[3:4])
    wgt.sigma(drop(crossprod(backsolve(lLTbk, btilde))))*tcrossprod(btilde)*
      prod(dnorm(u))
  }
  BCub <- matrix(robustlmm:::test_Hcubature_ninfInf2(funB, 4, 4, 5e-2), 2)
  ## BCub <- matrix(robustlmm:::test_GaussianQuadratureNd_ninfInf(funB, 4, 4, 13), 2)
  lTbk1Cub <- BCub/aCub
  ## return(lTbk1)
  (lTbk - lTbk1Cub)[-3]
}

optGH(lTbk[lower.tri(lTbk, diag = TRUE)])
## optCub(lTbk[lower.tri(lTbk, diag = TRUE)])

len <- 10
testvec <- seq(0, 1, length.out = len + 2)[-c(1, len + 2)]
testall <- expand.grid(testvec, testvec, testvec)
testall <- as.matrix(testall[testall[, 1] != testall[, 2] | testall[, 2] != testall[, 3] |
                               testall[, 1] != testall[, 3], ])
testall <- rbind(lTbk[lower.tri(lTbk, diag = TRUE)], testall)
if (len > 10) {
    require(parallel)
    setDefaultCluster(makeForkCluster(4))
    vals <- parApply(NULL, testall, 1, optGH)
} else {
    vals <- apply(testall, 1, optGH)
}

open3d()
rgl.clear()
for (i in 1:NROW(testall)) plotArrow(i, testall, vals)

lastTestall <- testall
lastVals <- vals
rm("sols")
for (i in 1:30) {
    lastTestall <- lastTestall - t(lastVals)
    lastVals <- apply(lastTestall, 1, optGH)
    lastCsVals <- colSums(abs(lastVals))
    lastIdx <- !is.na(lastCsVals) & lastCsVals > 0.005
    solIdx <- !is.na(lastCsVals) & lastCsVals < 0.005
    if (any(solIdx)) {
       if (exists("sols")) {
         sols <- rbind(sols, lastTestall[solIdx, , drop = FALSE])
       } else {
         sols <- lastTestall[solIdx, , drop = FALSE]
       }
    }
    if (!any(lastIdx)) {
      break
    }
    lastTestall <- lastTestall[lastIdx, , drop = FALSE]
    lastVals <- lastVals[, lastIdx, drop = FALSE]
    for (i in 1:NROW(lastTestall)) plotArrow(i, lastTestall, lastVals)
}

sols <- sols[order(sols[, 1]), ]
usols <- sols[1, , drop = FALSE]
count <- 1
for (i in 2:NROW(sols)) {
   dist <- sqrt(colSums((t(usols) - sols[i, ])^2))
   if (all(dist > 0.01)) {
       usols <- rbind(usols, sols[i, , drop = FALSE])
       count <- c(count, 1)
   } else {
      count[dist <= 0.01] <- count[dist <= 0.01] + 1
   }
}

points3d(usols, size = 10, color = "red")

round(usols, 5)
count
## Var1    Var2    Var3
## 821 0.00000 0.00000 0.00000
## 954 0.61105 0.28992 0.32525
optCub(usols[2, ])
optGH.wrong(usols[2, ])
