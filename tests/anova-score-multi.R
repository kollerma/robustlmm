## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## anova(test = "score") for several added independent variance
## components: cone (chi-bar) statistic of the whitened robust score.
require(robustlmm)

cps <- robustlmm:::.cone_proj_sq

## 1. cone projection: orthant cases and a skewed cone against a
##    brute-force search over a grid of nonnegative coefficients
I2 <- diag(2)
stopifnot(isTRUE(all.equal(cps(I2, c(1, -2)), 1)),
          isTRUE(all.equal(cps(I2, c(-1, -1)), 0)),
          isTRUE(all.equal(cps(I2, c(2, 3)), 13)))
A <- matrix(c(1, 0.3, 0.2, 1), 2)
b <- c(0.7, -1.1)
grid <- expand.grid(l1 = seq(0, 3, by = 0.002), l2 = seq(0, 3, by = 0.002))
res <- colSums((b - A %*% t(as.matrix(grid)))^2)
stopifnot(abs(cps(A, b) - (sum(b^2) - min(res))) < 1e-4)

## 2. scope: separate terms and a diag() term give the same design;
##    correlated slopes are rejected
set.seed(5)
J <- 30; m <- 8
d <- data.frame(g = factor(rep(seq_len(J), each = m)),
                x = rnorm(J * m), z = rnorm(J * m))
d$y <- 1 + 0.5 * d$x + rnorm(J)[d$g] + 0.6 * rnorm(J)[d$g] * d$x +
    rnorm(J * m)
f0  <- rlmer(y ~ x + (1 | g), d, method = "DASvar")
f1  <- suppressWarnings(rlmer(y ~ x + (1 | g) + (0 + x | g) + (0 + z | g),
                              d, method = "DASvar"))
f1d <- suppressWarnings(rlmer(y ~ x + diag(1 + x + z | g), d,
                              method = "DASvar"))
f1c <- suppressWarnings(rlmer(y ~ x + (1 + x | g), d, method = "DASvar"))
sc <- robustlmm:::.anova_score_scope(f0, f1)
scd <- robustlmm:::.anova_score_scope(f0, f1d)
stopifnot(is.matrix(sc$z_add), ncol(sc$z_add) == 2L,
          isTRUE(all.equal(unname(sc$z_add), cbind(d$x, d$z))),
          isTRUE(all.equal(unname(scd$z_add), unname(sc$z_add))),
          identical(sc$label, c("(x | g)", "(z | g)")))
e <- tryCatch(robustlmm:::.anova_score_scope(f0, f1c), error = identity)
stopifnot(inherits(e, "error"),
          grepl("exactly one variance parameter", conditionMessage(e)))

## 3. the test: Df = 2, identical statistic for both formulations,
##    reproducible under a seed, p-value in (0, 1]
a1 <- suppressWarnings(anova(f0, f1, test = "score", nsim = 19, seed = 3))
a2 <- suppressWarnings(anova(f0, f1d, test = "score", nsim = 19, seed = 3))
a3 <- suppressWarnings(anova(f0, f1, test = "score", nsim = 19, seed = 3))
p1 <- a1[["Pr(>=Score)"]][2L]
stopifnot(a1$Df[2L] == 2L,
          isTRUE(all.equal(a1$Score[2L], a2$Score[2L])),
          identical(p1, a3[["Pr(>=Score)"]][2L]),
          p1 > 0, p1 <= 1,
          any(grepl("2 added variance components", attr(a1, "heading"))),
          ncol(attr(a1, "boot")$s_j) == 2L)

## 4. the statistic is nonnegative and zero when the score points away
##    from the cone
st <- robustlmm:::.anova_score_stat_multi(f0, sc$z_add)
stopifnot(st$T >= 0, all(dim(st$M) == c(2L, 2L)),
          all(eigen(st$M, only.values = TRUE)$values > 0))

## 5. values, not just signs: each column's per-cluster contributions are
##    the scalar statistic's, and with one column the cone statistic is
##    max(0, S)^2 of the scalar test
for (k in 1:2) {
    s1 <- robustlmm:::.anova_score_stat(f0, sc$z_add[, k])
    stopifnot(isTRUE(all.equal(unname(st$s_j[, k]), unname(s1$s_j))))
}
st1 <- robustlmm:::.anova_score_stat_multi(f0, sc$z_add[, 1, drop = FALSE])
s1 <- robustlmm:::.anova_score_stat(f0, sc$z_add[, 1])
stopifnot(isTRUE(all.equal(st1$T, max(0, s1$S)^2)))

## 6. the sensitivity matrix: M_kl = a (v_k'v_l)^2 + c sum_i v_ki^2 v_li^2 is
##    d/dlambda E (v_k' psi(r))^2 at lambda = 0 for r ~ N(0, I + lambda
##    v_l v_l'); check (a, c) against a Monte Carlo derivative with common
##    random numbers
cst <- robustlmm:::.score_multi_consts(f0@rho.e)
set.seed(1)
nn <- 4; vk <- rnorm(nn); vl <- rnorm(nn)
Zs <- matrix(rnorm(5e5 * nn), ncol = nn)
Ef <- function(lam, v) {
    e <- eigen(diag(nn) + lam * tcrossprod(v), symmetric = TRUE)
    R <- e$vectors %*% (sqrt(e$values) * t(e$vectors))
    mean((matrix(f0@rho.e@psi(Zs %*% R), ncol = nn) %*% vk)^2)
}
h <- 0.02
for (v in list(vl, vk)) {
    num <- (Ef(h, v) - Ef(-h, v)) / (2 * h)
    cf <- cst[["a"]] * sum(vk * v)^2 + cst[["c"]] * sum(vk^2 * v^2)
    stopifnot(abs(num - cf) < 0.01 * abs(cf) + 0.005)
}

cat("anova-score-multi.R: all checks passed\n")
