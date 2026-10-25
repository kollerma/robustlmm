## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## F3b: small-J sandwich inference with a bias-reduced (CR2) robust
## sandwich. The decision rule was fixed before the run.
##
## Candidates per term (w: between, k = 1; f: between factor, k = 2;
## x: within, k = 1):
##   G1_chisq   G1 sandwich, chi-square(k)                 (current anova)
##   G1_FJ1     G1 sandwich, F(k, J - 1)                   (F3 baseline)
##   CR2_FJ1    CR2 sandwich, F(k, J - 1)
##   CR3_FJ1    CR3 (one-step leave-one-cluster-out) sandwich, F(k, J - 1)
##   CR2_HTZ    CR2 sandwich, approximate Hotelling T^2 with estimated df
##              eta (Pustejovsky and Tipton 2018; k = 1: Bell-McCaffrey t)
##   def_satt   default vcov, Satterthwaite F (anova(ddf = "satterthwaite"))
## plus 95% CI coverage of the intercept and the t slope (T5's coefficients)
## for G1 + t(J-1) and CR2 + t(BM df).
## The CR2 pieces are computed by the prototype below, which reproduces
## vcov(fit, type = "sandwich") to 1e-14 for correction = "G1" and a GLS
## CR2 computed in observation space to 2e-11 (checked before the run).

args <- commandArgs(trailingOnly = TRUE)
R_rep <- if (length(args) >= 1) as.integer(args[1]) else 1000L
ncores <- if (length(args) >= 2) as.integer(args[2]) else 5L
outdir <- "f3b_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cr2_pieces <- function(fit) {
    parts <- robustlmm:::implicitIF(fit)
    p <- fit@pp$p; n <- fit@pp$n
    eta <- fit@pp$eta; if (length(eta) != n) eta <- rep(1, n)
    De <- eta * parts$D_e
    AX <- parts$A_X; AZ <- parts$A_ZUb
    Pen <- as.matrix(parts$Lambda_ratio %*% parts$D_b_act)
    Juu <- crossprod(AZ, De * AZ) + Pen
    Juu <- (Juu + t(Juu)) / 2
    cl <- as.integer(factor(robustlmm:::resolveCluster(fit, NULL, n)))
    J <- max(cl)
    ## u -> cluster via the grouping-factor levels (single factor here)
    Zt <- as(fit@pp$Zt, "TsparseMatrix")
    u_cl <- integer(nrow(Zt)); u_cl[Zt@i + 1L] <- cl[Zt@j + 1L]
    psi <- fit@rho.e@psi(parts$r_std)
    Aj <- vector("list", J); gj <- matrix(0, p, J)
    for (j in seq_len(J)) {
        ix <- which(cl == j); uj <- which(u_cl == j)
        Xj <- AX[ix, , drop = FALSE]; Zj <- AZ[ix, uj, drop = FALSE]
        dj <- De[ix]
        Bxu <- crossprod(Xj, dj * Zj)
        Aj[[j]] <- crossprod(Xj, dj * Xj) -
            Bxu %*% solve(Juu[uj, uj, drop = FALSE], t(Bxu))
        gj[, j] <- colSums(Xj * (eta[ix] * psi[ix]))
    }
    list(A = Reduce(`+`, Aj), Aj = Aj, gj = gj, sigma = parts$sigma,
         J = J, p = p)
}
mpow <- function(M, pw, tol = 1e-10) {
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    v <- ifelse(e$values > tol, e$values^pw, 0)
    e$vectors %*% (v * t(e$vectors))
}
cr_vcov <- function(pc) {
    S <- mpow(pc$A, -1/2); p <- pc$p; J <- pc$J
    IF2 <- IF3 <- IF1 <- matrix(0, p, J); Mj <- vector("list", J)
    for (j in seq_len(J)) {
        M <- S %*% pc$Aj[[j]] %*% S; M <- (M + t(M)) / 2; Mj[[j]] <- M
        Sg <- S %*% pc$gj[, j]
        IF1[, j] <- S %*% Sg
        IF2[, j] <- S %*% mpow(diag(p) - M, -1/2) %*% Sg
        IF3[, j] <- S %*% mpow(diag(p) - M, -1) %*% Sg
    }
    s2 <- pc$sigma^2
    list(G1 = s2 * tcrossprod(IF1) * J / (J - 1), CR2 = s2 * tcrossprod(IF2),
         CR3 = s2 * tcrossprod(IF3), S = S, Mj = Mj)
}
htz_eta <- function(cv, J, p, L) {
    L <- rbind(L); k <- nrow(L)
    Ainv <- cv$S %*% cv$S
    Lt <- mpow(L %*% Ainv %*% t(L), -1/2) %*% L
    ct <- cv$S %*% t(Lt)
    u <- array(0, c(p, J, k)); d <- array(0, c(J, k, k))
    for (j in seq_len(J)) {
        M <- cv$Mj[[j]]; a <- mpow(diag(p) - M, -1/2) %*% ct
        uu <- M %*% a; u[, j, ] <- uu
        d[j, , ] <- crossprod(a, uu)
    }
    G <- function(s, t) diag(d[, s, t], J) - crossprod(u[, , s], u[, , t])
    tot <- 0
    for (s in seq_len(k)) for (t in seq_len(k)) {
        Gst <- G(s, t)
        tot <- tot + sum(G(s, s) * G(t, t)) + sum(Gst * t(Gst))
    }
    k * (k + 1) / tot
}

one_rep <- function(seed, J, scen) {
    set.seed(seed)
    m <- 10L
    g <- factor(rep(seq_len(J), each = m))
    t <- rep(0:9, J)
    w <- rnorm(J)[g]
    f <- factor(rep_len(c("a", "b", "c"), J))[g]
    x <- rnorm(J * m)
    y <- 250 + 10 * t + rnorm(J, sd = 25)[g] + rnorm(J, sd = 6)[g] * t +
        rnorm(J * m, sd = 25)
    if (scen == "contam") {
        bad <- sample.int(J * m, round(0.05 * J * m))
        y[bad] <- y[bad] + 100
    }
    d <- data.frame(g = g, t = t, w = w, f = f, x = x, y = y)
    nout <- 6 * 3 + 3 + 4 + 2
    fit <- tryCatch(suppressMessages(suppressWarnings(
        rlmer(y ~ t + w + f + x + (t | g), d, method = "DASvar"))),
        error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, nout))
    b <- fixef(fit); nm <- names(b); p <- length(b)
    pc <- tryCatch(cr2_pieces(fit), error = function(e) NULL)
    if (is.null(pc)) return(rep(NA_real_, nout))
    cv <- cr_vcov(pc)
    cols <- list(w = "w", f = c("fb", "fc"), x = "x")
    pv <- NULL
    for (tt in names(cols)) {
        cc <- match(cols[[tt]], nm); k <- length(cc)
        L <- matrix(0, k, p); L[cbind(seq_len(k), cc)] <- 1
        Wf <- function(V) as.numeric(crossprod(b[cc], solve(V[cc, cc, drop = FALSE], b[cc])))
        W1 <- Wf(cv$G1); W2 <- Wf(cv$CR2); W3 <- Wf(cv$CR3)
        eta <- htz_eta(cv, J, p, L)
        pv <- c(pv,
                pchisq(W1, k, lower.tail = FALSE),
                pf(W1 / k, k, J - 1, lower.tail = FALSE),
                pf(W2 / k, k, J - 1, lower.tail = FALSE),
                pf(W3 / k, k, J - 1, lower.tail = FALSE),
                if (eta > k - 1) pf(W2 * (eta - k + 1) / (eta * k), k, eta - k + 1,
                                    lower.tail = FALSE) else NA_real_,
                eta)
    }
    ## default vcov + Satterthwaite (per-term F)
    aS <- tryCatch(suppressWarnings(anova(fit, ddf = "satterthwaite")),
                   error = function(e) NULL)
    ps <- if (is.null(aS)) rep(NA_real_, 3) else aS[c("w", "f", "x"), "Pr(>F)"]
    ## coverage of intercept and t: G1 + t(J-1), CR2 + t(BM df)
    truth <- c(250, 10)
    cov_G1 <- as.numeric(abs(b[1:2] - truth) <= qt(0.975, J - 1) * sqrt(diag(cv$G1))[1:2])
    bm <- vapply(1:2, function(i) { L <- matrix(0, 1, p); L[1, i] <- 1
        htz_eta(cv, J, p, L) }, numeric(1))
    cov_CR2 <- as.numeric(abs(b[1:2] - truth) <= qt(0.975, pmax(bm, 1)) * sqrt(diag(cv$CR2))[1:2])
    c(pv, ps, cov_G1, cov_CR2, bm[2], bm[1])[seq_len(nout)]
}

cells <- rbind(expand.grid(J = c(8L, 12L, 18L, 30L), scen = "clean",
                           stringsAsFactors = FALSE),
               data.frame(J = 18L, scen = "contam"))
cn <- c(paste(rep(c("w", "f", "x"), each = 6),
              c("G1_chisq", "G1_FJ1", "CR2_FJ1", "CR3_FJ1", "CR2_HTZ", "eta"),
              sep = "_"),
        paste0("def_satt_", c("w", "f", "x")),
        "covG1_int", "covG1_t", "covCR2_int", "covCR2_t", "bm_t", "bm_int")
if (sys.nframe() == 0L && !isTRUE(getOption("f3b.source_only"))) {
    library(parallel)
    cl <- makeCluster(ncores)
    clusterCall(cl, function() {
        suppressMessages(library(robustlmm))
        NULL
    })
    clusterExport(cl, c("one_rep", "cr2_pieces", "cr_vcov", "mpow", "htz_eta"))
    for (i in seq_len(nrow(cells))) {
        J <- cells$J[i]; scen <- cells$scen[i]
        f <- file.path(outdir, sprintf("cell_J%d_%s.rds", J, scen))
        if (file.exists(f)) next
        t0 <- proc.time()
        res <- parLapplyLB(cl, 950000L + 10000L * i + seq_len(R_rep),
                           function(s, J, scen) one_rep(s, J, scen),
                           J = J, scen = scen)
        res <- do.call(rbind, res); colnames(res) <- cn
        el <- proc.time() - t0
        saveRDS(list(J = J, scen = scen, res = res, wall = el[["elapsed"]],
                     ncores = ncores), f)
        cat(sprintf("cell J=%d %s done in %.0f s\n", J, scen, el[["elapsed"]]))
    }
    stopCluster(cl)

    summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
        o <- readRDS(file.path(outdir, sprintf("cell_J%d_%s.rds", cells$J[i],
                                               cells$scen[i])))
        r <- o$res
        do.call(rbind, lapply(c("w", "f", "x"), function(tt) {
            rej <- function(v) mean(v < 0.05, na.rm = TRUE)
            data.frame(J = o$J, scen = o$scen, term = tt,
                       n_ok = sum(is.finite(r[, paste0(tt, "_G1_chisq")])),
                       G1_chisq = rej(r[, paste0(tt, "_G1_chisq")]),
                       G1_FJ1 = rej(r[, paste0(tt, "_G1_FJ1")]),
                       CR2_FJ1 = rej(r[, paste0(tt, "_CR2_FJ1")]),
                       CR3_FJ1 = rej(r[, paste0(tt, "_CR3_FJ1")]),
                       CR2_HTZ = rej(r[, paste0(tt, "_CR2_HTZ")]),
                       def_satt = rej(r[, paste0("def_satt_", tt)]),
                       eta_med = stats::median(r[, paste0(tt, "_eta")], na.rm = TRUE),
                       covG1_int = mean(r[, "covG1_int"], na.rm = TRUE),
                       covG1_t = mean(r[, "covG1_t"], na.rm = TRUE),
                       covCR2_int = mean(r[, "covCR2_int"], na.rm = TRUE),
                       covCR2_t = mean(r[, "covCR2_t"], na.rm = TRUE),
                       cpu_min = o$wall * o$ncores / 60)
        }))
    }))
    print(summ, digits = 3)
    saveRDS(summ, file.path(outdir, "summary.rds"))
}
