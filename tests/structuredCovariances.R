## Tests for lme4 >= 2.0-0 structured random-effect covariances:
##   diag() (Phase 1), cs() compound symmetric (Phase 2),
##   ar1() autoregressive (Phase 3).
##
## Structured covariances need lme4 >= 2.0-0; skip otherwise so the test
## stays a no-op on older lme4.
require(robustlmm)

if (packageVersion("lme4") < "2.0.0") {
    cat("lme4 < 2.0-0: structured covariances not available, skipping.\n")
} else {

set.seed(2)

## ---- diag(): Phase 1 still works -----------------------------------
rdiag <- rlmer(Reaction ~ Days + diag(Days | Subject), sleepstudy,
               method = "DASvar")
stopifnot(all(is.finite(theta(rdiag))),
          ## diagonal structure => two variance parameters, zero covariance
          length(theta(rdiag)) == 2L)

## ---- cs(): 2x2 is identical to the unstructured fit ----------------
## A 2x2 compound-symmetric covariance has a single correlation, so it
## coincides with the unstructured (us) covariance: the projection must
## be an exact identity and reproduce the us fit bit-for-bit.
for (method in c("DASvar", "DAStau")) {
    rcs <- rlmer(Reaction ~ Days + cs(Days | Subject), sleepstudy,
                 method = method)
    rus <- rlmer(Reaction ~ Days + (Days | Subject), sleepstudy,
                 method = method)
    stopifnot(max(abs(theta(rcs) - theta(rus))) < 1e-7,
              max(abs(fixef(rcs) - fixef(rus))) < 1e-7,
              abs(sigma(rcs) - sigma(rus)) < 1e-7)
}

## ---- cs(): 3x3 enforces the compound-symmetric structure -----------
nsubj <- 40
truesd <- c(2, 1.5, 1.2)
rho <- 0.4
R <- matrix(rho, 3, 3)
diag(R) <- 1
Lt <- t(chol(diag(truesd) %*% R %*% diag(truesd)))
dat <- expand.grid(f = factor(c("a", "b", "c")), rep = 1:5,
                   subj = factor(seq_len(nsubj)))
bvals <- t(Lt %*% matrix(rnorm(3 * nsubj), 3, nsubj))
dat$y <- 10 + bvals[cbind(as.integer(dat$subj), as.integer(dat$f))] +
    rnorm(nrow(dat), sd = 0.8)

rcs3 <- rlmer(y ~ 1 + cs(0 + f | subj), dat, method = "DASvar")
mcs3 <- lmer(y ~ 1 + cs(0 + f | subj), dat, REML = TRUE)

cm <- attr(VarCorr(rcs3)$subj, "correlation")
sd3 <- attr(VarCorr(rcs3)$subj, "stddev")
sdLmer <- attr(VarCorr(mcs3)$subj, "stddev")
## all off-diagonal correlations equal (compound symmetric)
stopifnot(isTRUE(all.equal(rep(cm[2, 1], 3), cm[upper.tri(cm)],
                           tolerance = 1e-6)),
          ## on clean data the robust fit tracks the classical cs fit;
          ## the marginal sds are preserved by the projection
          all(abs(sd3 - sdLmer) / sdLmer < 0.15),
          ## genuine convergence, no spurious false-convergence warning
          rcs3@optinfo$conv$opt == 0)

## the fitted covariance must reproduce lme4's parametrisation: VarCorr
## (which reads the reCovs 'par' slots) must match the theta-implied cov
vc <- VarCorr(rcs3)$subj
covFromTheta <- {
    th <- theta(rcs3)
    L <- matrix(0, 3, 3)
    L[lower.tri(L, diag = TRUE)] <- th
    sigma(rcs3)^2 * tcrossprod(L)
}
stopifnot(max(abs(as.matrix(vc) - covFromTheta)) < 1e-6)

## ---- ar1(): autoregressive lag structure --------------------------
## The homogeneous ar1 (the ar1() default) has a single variance and a
## lag-decaying correlation Cor(i, j) = rho^|i - j|. Simulate a 4-level
## within-subject factor and check the fit reproduces that structure.
nc <- 4
tsig <- 1.5
trho <- 0.6
Rar <- trho^abs(outer(seq_len(nc), seq_len(nc), "-"))
Lar <- t(chol(tsig^2 * Rar))
dar <- expand.grid(f = factor(letters[seq_len(nc)]), rep = 1:4,
                   subj = factor(seq_len(60)))
bar <- t(Lar %*% matrix(rnorm(nc * 60), nc, 60))
dar$y <- 5 + bar[cbind(as.integer(dar$subj), as.integer(dar$f))] +
    rnorm(nrow(dar), sd = 0.7)

rar <- rlmer(y ~ 1 + ar1(0 + f | subj), dar, method = "DASvar")
mar <- lmer(y ~ 1 + ar1(0 + f | subj), dar, REML = TRUE)

vcar <- VarCorr(rar)$subj
cmar <- unname(attr(vcar, "correlation"))
sdar <- as.numeric(attr(vcar, "stddev"))
stopifnot(## homogeneous: all marginal sds equal
          max(abs(sdar - sdar[1])) < 1e-8,
          ## autoregressive: row 1 correlations equal rho^lag
          max(abs(cmar[1, ] - cmar[1, 2]^(0:(nc - 1)))) < 1e-8,
          ## tracks the classical ar1 fit on clean data
          max(abs(sdar - as.numeric(attr(VarCorr(mar)$subj, "stddev"))) /
              as.numeric(attr(VarCorr(mar)$subj, "stddev"))) < 0.15,
          rar@optinfo$conv$opt == 0)

## ---- ar1() with a NEGATIVE rho: the ECME profile-rho step recovers it
## (the former log-correlation least-squares projection discarded negative
## off-diagonals and could only return rho >= 0). -------------------------
set.seed(101)
trhoN <- -0.5
RarN <- trhoN^abs(outer(seq_len(nc), seq_len(nc), "-"))
LarN <- t(chol(tsig^2 * RarN))
darN <- expand.grid(f = factor(letters[seq_len(nc)]), rep = 1:4,
                    subj = factor(seq_len(60)))
barN <- t(LarN %*% matrix(rnorm(nc * 60), nc, 60))
darN$y <- 5 + barN[cbind(as.integer(darN$subj), as.integer(darN$f))] +
    rnorm(nrow(darN), sd = 0.7)
rarN <- rlmer(y ~ 1 + ar1(0 + f | subj), darN, method = "DASvar")
rhoN <- unname(attr(VarCorr(rarN)$subj, "correlation"))[1, 2]
stopifnot(rhoN < 0,                       # sign recovered (old method could not)
          abs(rhoN - trhoN) < 0.15)       # close to truth

## ---- saturated structured block: cluster size = RE dimension -------
## A cs()/ar1() block of dimension nc needs MORE than nc observations per
## cluster; with exactly one observation per (subject, visit) the cluster
## size equals nc, the block is not identifiable, and the robust DAS fit
## would fail deep inside the projection with a cryptic error. The guard
## must instead stop up front with an informative message.
datSat <- expand.grid(f = factor(letters[seq_len(nc)]), rep = 1L,
                      subj = factor(seq_len(40)))
barS <- t(Lar %*% matrix(rnorm(nc * 40), nc, 40))
datSat$y <- 5 + barS[cbind(as.integer(datSat$subj), as.integer(datSat$f))] +
    rnorm(nrow(datSat), sd = 0.7)
stopifnot(min(tabulate(datSat$subj)) == nc)   # genuinely saturated
errSat <- tryCatch(
    rlmer(y ~ 1 + ar1(0 + f | subj), datSat, method = "DASvar"),
    error = function(e) conditionMessage(e))
stopifnot(is.character(errSat),
          grepl("not identifiable without", errSat, fixed = TRUE),
          grepl("more than nc observations per cluster", errSat, fixed = TRUE))

## ---- non-saturated companion still fits (cluster size > nc) ---------
datRep <- expand.grid(f = factor(letters[seq_len(nc)]), rep = 1:2,
                      subj = factor(seq_len(40)))
barR <- t(Lar %*% matrix(rnorm(nc * 40), nc, 40))
datRep$y <- 5 + barR[cbind(as.integer(datRep$subj), as.integer(datRep$f))] +
    rnorm(nrow(datRep), sd = 0.7)
stopifnot(min(tabulate(datRep$subj)) > nc)
rRep <- rlmer(y ~ 1 + ar1(0 + f | subj), datRep, method = "DASvar")
stopifnot(all(is.finite(theta(rRep))))

cat("structured covariance tests (diag, cs, ar1) passed.\n")

}
