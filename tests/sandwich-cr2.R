## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Bias-reduced (CR2 / CR3) robust cluster sandwich, its Bell-McCaffrey /
## approximate-Hotelling degrees of freedom, and the Mallows-weighted meat.
require(robustlmm)

msqrt <- function(M, pw) {
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    e$vectors %*% (e$values^pw * t(e$vectors))
}

## 1. With classical psi-functions the robust CR2 is the CR2 estimator of
##    generalised least squares (Bell and McCaffrey 2002), computed here in
##    observation space from the fitted marginal covariance, and its df are
##    the Satterthwaite df of that estimator under the working model.
set.seed(3)
J <- 10; m <- 6
g <- factor(rep(seq_len(J), each = m)); t <- rep(0:5, J)
w <- rnorm(J)[g]; x <- rnorm(J * m)
y <- 1 + t + rnorm(J, sd = 2)[g] + rnorm(J, sd = .5)[g] * t + rnorm(J * m)
d <- data.frame(g, t, w, x, y)
fitC <- rlmer(y ~ t + w + x + (t | g), d, rho.e = cPsi, rho.b = cPsi,
              rho.sigma.e = cPsi, rho.sigma.b = cPsi, method = "DASvar")
V2 <- vcov(fitC, type = "sandwich", correction = "CR2")
X <- getME(fitC, "X"); s <- sigma(fitC)
Zub <- as.matrix(t(fitC@pp$Zt) %*% fitC@pp$U_b)
e <- as.numeric(y - X %*% fixef(fitC))
idx <- split(seq_along(y), g)
Xt <- do.call(rbind, lapply(idx, function(ix)
    msqrt(s^2 * (diag(length(ix)) + tcrossprod(Zub[ix, , drop = FALSE])),
          -1/2) %*% X[ix, ]))
et <- unlist(lapply(idx, function(ix)
    msqrt(s^2 * (diag(length(ix)) + tcrossprod(Zub[ix, , drop = FALSE])),
          -1/2) %*% e[ix]))
M <- solve(crossprod(Xt)); H <- Xt %*% M %*% t(Xt)
Aadj <- lapply(idx, function(ix) msqrt(diag(length(ix)) - H[ix, ix], -1/2))
meat <- Reduce(`+`, lapply(seq_along(idx), function(j) {
    ix <- idx[[j]]
    tcrossprod(t(Xt[ix, ]) %*% Aadj[[j]] %*% et[ix])
}))
Vn <- M %*% meat %*% M
stopifnot(max(abs(as.matrix(V2) - Vn)) / max(abs(Vn)) < 1e-8,
          attr(V2, "correction") == "CR2", attr(V2, "n.clusters") == J)
cw <- as.numeric(colnames(X) == "w")
Q <- matrix(0, nrow(Xt), J)
for (j in seq_len(J)) {
    ix <- idx[[j]]
    Q[ix, j] <- t(cw %*% M %*% t(Xt[ix, ]) %*% Aadj[[j]])
}
Gm <- t(Q) %*% (diag(nrow(Xt)) - H) %*% Q
df_n <- sum(diag(Gm))^2 / sum(Gm^2)
df_p <- robustlmm:::.cr2_df(V2, rbind(cw))
stopifnot(abs(df_p - df_n) < 1e-6 * df_n)

## 2. G1 is unchanged by the refactoring: equals J/(J-1) times "none", and
##    CR2/CR3 never warn (the small-J warning concerns G1 / none).
V1 <- suppressWarnings(vcov(fitC, type = "sandwich"))
V0 <- suppressWarnings(vcov(fitC, type = "sandwich", correction = "none"))
stopifnot(max(abs(V1 - V0 * J / (J - 1))) < 1e-12 * max(abs(V1)))
r <- tryCatch(vcov(fitC, type = "sandwich", correction = "CR3"),
              warning = function(w) "warned")
stopifnot(!identical(r, "warned"))

## 3. Mallows design weights enter the sandwich meat: the "none" sandwich
##    equals sigma^2 A^-1 (sum_j s_j s_j') A^-1 with s_j = sum_i eta_i x_i psi_i,
##    and differs from the unweighted meat when eta < 1.
set.seed(11)
J <- 30; m <- 6
g <- factor(rep(seq_len(J), each = m))
x <- rnorm(J * m); x[sample.int(J * m, 10)] <- rnorm(10, 6, 1)
y <- 1 + 0.5 * x + rnorm(J)[g] + rnorm(J * m)
fitM <- rlmer(y ~ x + (1 | g), data.frame(g, x, y), design.weights = "mcd",
              method = "DASvar")
eta <- fitM@pp$eta
stopifnot(any(eta < 1))
parts <- robustlmm:::implicitIF(fitM)
A <- solve(solve(as.matrix(parts$Jbb)) +
           as.matrix(parts$Jbb_inv_Jbu %*% solve(parts$S, t(parts$Jbb_inv_Jbu))))
psi <- fitM@rho.e@psi(parts$r_std)
sj <- rowsum(parts$A_X * (eta * psi), g)
Vm <- parts$sigma^2 * solve(A, t(solve(A, t(crossprod(sj)))))
VM <- suppressWarnings(vcov(fitM, type = "sandwich", correction = "none"))
stopifnot(max(abs(VM - Vm)) < 1e-8 * max(abs(Vm)))
sj0 <- rowsum(parts$A_X * psi, g)
Vm0 <- parts$sigma^2 * solve(A, t(solve(A, t(crossprod(sj0)))))
stopifnot(Vm0[2, 2] > 1.2 * Vm[2, 2])
## caseweightIF carries the same weights
cwM <- caseweightIF(fitM)
i1 <- which.min(eta)
stopifnot(abs(cwM$Gb_cw[2, i1] - eta[i1] * psi[i1] * parts$A_X[i1, 2]) < 1e-10)

## 4. Crossed random effects: CR2 is not defined and says so; G1 works.
set.seed(5)
dC <- expand.grid(a = factor(1:8), b = factor(1:6))
dC$y <- rnorm(8)[dC$a] + rnorm(6)[dC$b] + rnorm(nrow(dC))
fitX <- rlmer(y ~ 1 + (1 | a) + (1 | b), dC, method = "DASvar")
e2 <- tryCatch(suppressWarnings(vcov(fitX, type = "sandwich", cluster = "a",
                                     correction = "CR2")),
               error = function(e) conditionMessage(e))
stopifnot(is.character(e2), grepl("nested within one cluster", e2))

## the fallback: anova() with the sandwich uses G1 with F(1, J - 1) and
## says so; the cluster warning names the crossed factor once
wx <- character()
aX <- withCallingHandlers(anova(fitX, vcov_type = "sandwich", cluster = "a"),
                          warning = function(w) {
                              wx <<- c(wx, conditionMessage(w))
                              invokeRestart("muffleWarning")
                          })
stopifnot(all(aX$DenDF == 7), length(wx) == 1L, grepl("b not nested", wx))

## 5. CR2 does not depend on how the fixed effects are parametrised:
##    shifting and rescaling covariates (x -> 1e4 x + 2000, w -> w + 2000,
##    a badly conditioned Jacobian) transforms the covariance by T^-1 and
##    leaves every contrast's df unchanged. (A relative eigenvalue cut-off
##    in the whitening root used to drop a direction here: df = Inf.)
fitR <- rlmer(y ~ t + w + x + (t | g), d, method = "DASvar")
d2 <- transform(d, x = 1e4 * x + 2000, w = w + 2000)
fitR2 <- suppressWarnings(rlmer(y ~ t + w + x + (t | g), d2,
                                method = "DASvar"))
Va <- vcov(fitR, type = "sandwich", correction = "CR2")
Vb <- vcov(fitR2, type = "sandwich", correction = "CR2")
Tm <- diag(4); Tm[1, 3:4] <- 2000; Tm[4, 4] <- 1e4
Ti <- solve(Tm)
Vt <- Ti %*% as.matrix(Va) %*% t(Ti)
stopifnot(max(abs(Vt - as.matrix(Vb)) / sqrt(outer(diag(Vt), diag(Vt)))) < 1e-5)
dfa <- sapply(1:4, function(i) robustlmm:::.cr2_df(Va, Ti[i, , drop = FALSE]))
dfb <- sapply(1:4, function(i) robustlmm:::.cr2_df(Vb, diag(4)[i, , drop = FALSE]))
stopifnot(all(is.finite(dfb)), max(abs(dfa - dfb) / dfa) < 1e-4)

## 6. Approximate-Hotelling df for k >= 2 contrasts equal the observation-
##    space computation of Pustejovsky and Tipton (2018) for the GLS CR2
##    (classical psi, as in 1): eta = k (k + 1) / sum_st Var(Z_st).
htz <- function(L) {
    k <- nrow(L)
    C <- msqrt(L %*% M %*% t(L), -1/2) %*% L
    Q <- lapply(seq_len(k), function(s) {
        Qs <- matrix(0, nrow(Xt), J)
        for (j in seq_len(J)) {
            ix <- idx[[j]]
            Qs[ix, j] <- t(C[s, , drop = FALSE] %*% M %*% t(Xt[ix, ]) %*%
                           Aadj[[j]])
        }
        Qs
    })
    IH <- diag(nrow(Xt)) - H
    G <- function(s, t) t(Q[[s]]) %*% IH %*% Q[[t]]
    tot <- 0
    for (s in seq_len(k)) for (t in seq_len(k))
        tot <- tot + sum(G(s, s) * G(t, t)) + sum(diag(G(s, t) %*% G(s, t)))
    k * (k + 1) / tot
}
J <- 10
for (L in list(rbind(c(0, 0, 1, 0), c(0, 0, 0, 1)),
               rbind(c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1)),
               rbind(c(0, 1, -1, 0), c(1, 0, 0, 2)))) {
    e1 <- htz(L); e2 <- robustlmm:::.cr2_df(V2, L)
    stopifnot(abs(e1 - e2) < 1e-6 * e1)
}

## 7. A factor level present in a single cluster: that direction has no
##    estimable sandwich variance. vcov warns, the level's df / interval /
##    test are NA, the other coefficients keep theirs.
d3 <- d
d3$f <- factor(ifelse(d3$g == "1", "B", "A"))
d3$h <- factor(ifelse(as.integer(d3$g) <= 5, "lo", "hi"))
fitD <- rlmer(y ~ t + h + f + (1 | g), d3, method = "DASvar")
wD <- character()
VD <- withCallingHandlers(vcov(fitD, type = "sandwich", correction = "CR2"),
                          warning = function(w) {
                              wD <<- c(wD, conditionMessage(w))
                              invokeRestart("muffleWarning")
                          })
stopifnot(length(wD) == 1L, grepl("single cluster", wD),
          ncol(attr(VD, "cr2")$degen) == 1L)
dfD <- sapply(1:4, function(i) robustlmm:::.cr2_df(VD, diag(4)[i, , drop = FALSE]))
stopifnot(is.na(dfD[4]), all(is.finite(dfD[1:3])))
ciD <- suppressWarnings(confint(fitD, vcov_type = "sandwich"))
stopifnot(all(is.na(ciD["fB", ])), all(is.finite(ciD[1:3, ])))
aD <- suppressWarnings(anova(fitD, vcov_type = "sandwich"))
stopifnot(is.na(aD["f", "Pr(>F)"]), all(is.finite(aD[1:3, "Pr(>F)"])))

## 8. Nested design: with cluster = the top level every random effect is
##    nested within a cluster, so CR2 applies without a warning, and
##    anova() / confint() accept the cluster.
set.seed(4)
dn <- expand.grid(i = 1:5, class = 1:3, school = 1:12)
dn$school <- factor(dn$school)
dn$class <- factor(paste(dn$school, dn$class))
dn$x <- rnorm(nrow(dn)); dn$z <- rnorm(12)[dn$school]
dn$y <- 1 + dn$x + 0.5 * dn$z + rnorm(12)[dn$school] +
    rnorm(36, sd = .7)[dn$class] + rnorm(nrow(dn))
fitN <- rlmer(y ~ x + z + (1 | school/class), dn, method = "DASvar")
wN <- character()
aN <- withCallingHandlers(anova(fitN, vcov_type = "sandwich",
                                cluster = "school"),
                          warning = function(w) {
                              wN <<- c(wN, conditionMessage(w))
                              invokeRestart("muffleWarning")
                          })
cN <- withCallingHandlers(confint(fitN, vcov_type = "sandwich",
                                  cluster = "school"),
                          warning = function(w) {
                              wN <<- c(wN, conditionMessage(w))
                              invokeRestart("muffleWarning")
                          })
stopifnot(length(wN) == 0L, attr(cN, "correction") == "CR2",
          all(is.finite(aN$DenDF)), all(aN$DenDF < 12),
          isTRUE(all.equal(attr(cN, "df"), aN$DenDF)))
## without cluster the model's two grouping factors are ambiguous
e3 <- tryCatch(anova(fitN, vcov_type = "sandwich"),
               error = function(e) conditionMessage(e))
stopifnot(is.character(e3), grepl("specify `cluster`", e3))

## 9. A 5-level between-cluster factor at J = 10: the approximate-Hotelling
##    df leave fewer than one denominator df, and anova() says so.
set.seed(9)
d5 <- d
d5$f5 <- factor(rep(letters[1:5], each = 2))[as.integer(d5$g)]
fit5 <- rlmer(y ~ t + f5 + (1 | g), d5, method = "DASvar")
w5 <- character()
a5 <- withCallingHandlers(anova(fit5, vcov_type = "sandwich"),
                          warning = function(w) {
                              w5 <<- c(w5, conditionMessage(w))
                              invokeRestart("muffleWarning")
                          })
stopifnot(any(grepl("^f5: ", w5)),
          is.na(a5["f5", "DenDF"]) || a5["f5", "DenDF"] < 1)

cat("sandwich-cr2.R: all checks passed\n")
