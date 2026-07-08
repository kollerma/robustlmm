## WS16 step 1 (PLAN-WS16-satterthwaite.md §2, §3.1): the IF-based
## covariance of the variance parameters (sigma-hat, theta-hat) -- the
## "A" matrix of the Satterthwaite df formula
##     df(c) = 2 v(c)^2 / (g' A g).
##
## A is the cluster sandwich of the (sigma, theta) rows of the joint
## M-estimation linearization: with s_j the per-cluster contribution
## to the full score (F_beta, F_u, F_sigma, F_theta) at the fit
## (columns of .scoreByCluster, sum_j s_j = 0), the first-order
## expansion par-hat - par ~= -Jpar^{-1} sum_j s_j gives
##     Cov(par-hat) ~= sum_j (Jpar^{-1} s_j)(Jpar^{-1} s_j)'
## (G1 small-sample factor J/(J-1) optional), of which .vcov_varpar
## returns the (sigma, theta) block, transformed from the internal
## log-sigma coordinate.
##
## Why not simpler routes (validated 2026-06-11/12, see
## IF-thread1/thread1_validation/ws16_satterthwaite_anchor.R):
## - the model-delta IF Cov(y) IF' overstates Var(sigma-hat) by
##   exactly 2 on clean Gaussian data (quadratic-functional effect:
##   the linearized local-shift delta drops the factor 2 of
##   Var(quadratic form)) -- df comes out exactly half of lmerTest's;
## - caseweight = local-shift x raw residual is a linear-score-only
##   identity; applied to the quadratic sigma/theta scores it inflates
##   Var(sigma-hat) ~12x.
## The score sandwich E[s s'] carries both effects correctly.
##
## Grouping structure (WS17): the cluster sandwich needs a partition of
## {observations} U {random effects} into INDEPENDENT clusters. Such a
## partition exists iff some grouping factor G contains every random-
## effect factor by nesting (A nested in B iff each level of A meets a
## single level of B). Supported:
##   - a single grouping factor (G = it);
##   - nested designs, e.g. (1|school/class): G = the coarsest factor
##     (school); each cluster owns its observations AND all the random
##     effects nested in it (school-RE plus the class:school-REs).
## Crossed designs (e.g. (1|subj)+(1|item)) have no such G -- the score
## contributions are not independent under any single grouping -- and
## still stop() (WS17 Phase 2). Mallows design weights ARE supported
## (WS11 step 5): the e-side score rows carry eta, matching the
## estimator's eta-weighted equations.

## The top-level independent cluster factor G, or NULL if the design is
## crossed. For a single factor it is that factor; for a nested design it
## is the coarsest factor that all others nest into.
.top_cluster <- function(fit) {
    fl <- fit@flist
    if (length(fl) == 1L)
        return(list(g = fl[[1L]], name = names(fl)[1L]))
    nested_in <- function(A, B) all(rowSums(table(A, B) > 0L) == 1L)
    ## try coarsest (fewest levels) first: G must contain all others
    for (b in order(vapply(fl, nlevels, integer(1)))) {
        B <- fl[[b]]
        if (all(vapply(seq_along(fl), function(a)
                       a == b || nested_in(fl[[a]], B), logical(1))))
            return(list(g = B, name = names(fl)[b]))
    }
    NULL
}

## Per-cluster decomposition of the full RSE score at the fitted
## parameters. Returns a (p + q + 1 + L) x J matrix S whose columns
## sum to .scoreVec(par0, fit) (~= 0 at the fit) -- that identity is
## pinned in tests/vcov-varpar.R.
.scoreByCluster <- function(fit) {
    stopifnot(is(fit, "rlmerMod"))
    tc <- .top_cluster(fit)
    if (is.null(tc))
        stop(".scoreByCluster needs a single or nested grouping ",
             "structure; this fit has crossed grouping factors (",
             paste(names(fit@flist), collapse = ", "),
             ") with no common top-level cluster.")
    ## WS11 step 5: the per-cluster score carries the Mallows design
    ## weights eta on its e-side rows (F_beta, F_u e-part, F_sigma),
    ## matching .scoreVec; the b-side penalty and F_theta are not
    ## eta-weighted (eta enters theta via tau_b). eta == 1 reproduces
    ## the unweighted score.
    eta <- fit@resp$eta
    g  <- tc$g                 # top-level independent cluster factor
    J  <- nlevels(g)
    gi <- as.integer(g)        # observation -> top cluster

    pp     <- fit@pp
    p      <- pp$p; q <- pp$q; n <- pp$n
    if (length(eta) != n) eta <- rep.int(1, n)
    L      <- length(getME(fit, "theta"))
    rho_e  <- fit@rho.e
    rho_b  <- fit@rho.b
    rho_sb <- fit@rho.sigma.b
    sigma  <- .sigma(fit)
    u_std  <- pp$b.s / sigma

    ## per-observation pieces, exactly as in .scoreVec (at the fit,
    ## r_std equals resp$wtres / sigma; recompute from components so
    ## the identity test is meaningful)
    U_e   <- pp$U_e
    X     <- pp$X
    Z     <- t(pp$Zt)
    r     <- as.numeric(fit@resp$y - fit@resp$offset -
                        X %*% .fixef(fit) - Z %*% (pp$U_b %*% pp$b.s))
    r_std <- as.numeric(U_e %*% r) / sigma
    psi_e <- rho_e@psi(r_std)

    U_eX <- as.matrix(U_e %*% X)             # n x p
    ZL   <- as.matrix(U_e %*% Z %*% pp$U_b)  # n x q

    ## u-index -> top cluster. Each random-effect coefficient (row of Zt)
    ## loads only on observations within a single top cluster (guaranteed
    ## by nesting), so map it to the cluster of any loading observation.
    ## This generalises the single-factor case (where it reduces to the
    ## level index) to nested designs at any depth, across all idx blocks.
    Zt_t    <- as(pp$Zt, "TsparseMatrix")
    u_level <- integer(q)
    u_level[Zt_t@i + 1L] <- gi[Zt_t@j + 1L]   # duplicates share a cluster

    S <- matrix(0, p + q + 1L + L, J)

    ## F_beta rows: obs i contributes eta[i] * U_eX[i, ] * psi_e[i]
    contrib_beta <- t(U_eX * (psi_e * eta))  # p x n
    for (j in seq_len(p))
        S[j, ] <- S[j, ] + rowsum(contrib_beta[j, ], gi)[, 1L]

    ## F_u rows: obs part ZL[i, ] * psi_e[i]; penalty part
    ## -(lambda_e / lambda_b_q) psi_b(u_std) belongs to the u-index's
    ## own cluster
    lambda_e   <- rho_e@EDpsi()
    lambda_b_q <- pp$D_b@x
    psi_b      <- .applyBlockPsi_b(fit, u_std, rho_b)
    contrib_u  <- t(ZL * (psi_e * eta))      # q x n  (e-side eta-weighted)
    for (j in seq_len(q)) {
        S[p + j, ] <- S[p + j, ] + rowsum(contrib_u[j, ], gi)[, 1L]
        S[p + j, u_level[j]] <- S[p + j, u_level[j]] -
            (lambda_e / lambda_b_q[j]) * psi_b[j]
    }

    ## F_sigma row: per-observation summands of the DAS scale equation
    tau_e     <- pp$tau_e()
    arg_sigma <- r_std / tau_e
    w_sigma   <- fit@rho.sigma.e@wgt(arg_sigma)
    v_sigma   <- eta * tau_e^2 * w_sigma * (arg_sigma^2 - pp$kappa_e)
    S[p + q + 1L, ] <- rowsum(v_sigma, gi)[, 1L]

    ## F_theta rows: per-RE-block summands, mirroring .scoreTheta_diag
    ## / .scoreTheta_block
    block_sizes <- vapply(fit@idx, nrow, integer(1))
    kappa_per_block <- pp$kappa_b
    if (all(block_sizes == 1L)) {
        tau_b  <- .getTauB_diag(fit, pp, rho_b, rho_sb)
        arg_th <- u_std / tau_b
        for (l in seq_along(rho_sb)) {
            is_in <- which(fit@ind == l)
            w_th  <- rho_sb[[l]]@wgt(arg_th[is_in])
            terms <- w_th * u_std[is_in]^2 -
                kappa_per_block[[l]] * w_th * tau_b[is_in]^2
            acc <- numeric(J)
            tt  <- rowsum(terms, u_level[is_in])
            acc[as.integer(rownames(tt))] <- tt[, 1L]
            S[p + q + 1L + l, ] <- acc
        }
    } else {
        Tb_full <- .getTbFull_block(fit, pp)
        row_off <- 0L
        for (type in seq_along(fit@idx)) {
            idx_mat <- fit@idx[[type]]
            s <- nrow(idx_mat); K <- ncol(idx_mat)
            psi_fn  <- rho_sb[[type]]@psi
            wgt_fn  <- rho_sb[[type]]@wgt
            kappa_t <- kappa_per_block[[type]]
            pat  <- fit@blocks[[type]] != 0
            Lind <- fit@blocks[[type]][pat]
            ord  <- order(Lind)
            n_th <- length(Lind)
            for (k in seq_len(K)) {
                bk_idx  <- as.vector(idx_mat[, k])
                cl_k    <- u_level[bk_idx[1L]]   # this level's top cluster
                b_block <- u_std[bk_idx]
                Tk      <- as.matrix(Tb_full[bk_idx, bk_idx,
                                             drop = FALSE])
                Lk      <- t(chol(Tk))
                T_b     <- forwardsolve(Lk, b_block)
                d_k     <- sum(T_b * T_b)
                if (s > 1L) {
                    w_eta   <- wgt_fn(d_k)
                    w_delta <- (psi_fn(d_k) -
                                psi_fn(d_k - s * kappa_t)) / s
                } else {
                    w_eta   <- wgt_fn(d_k)
                    w_delta <- w_eta * kappa_t
                }
                diff_k  <- (w_eta * tcrossprod(b_block) -
                            w_delta * Tk) / K
                entries <- diff_k[pat][ord]
                S[p + q + 1L + row_off + seq_len(n_th), cl_k] <-
                    S[p + q + 1L + row_off + seq_len(n_th), cl_k] +
                    entries
            }
            row_off <- row_off + n_th
        }
    }

    rownames(S) <- c(paste0("beta", seq_len(p)),
                     paste0("u", seq_len(q)),
                     "log_sigma",
                     paste0("theta", seq_len(L)))
    colnames(S) <- levels(g)
    S
}

## WS17 Phase 2 (crossed designs): a PER-OBSERVATION decomposition of the
## full RSE score. Returns a (p + q + 1 + L) x n matrix G whose columns
## sum to .scoreVec(par0, fit) (~= 0 at the fit). The e-side rows
## (F_beta, F_u e-part, F_sigma) are naturally per observation; the
## b-side penalty (F_u) and F_theta are per random effect, and are
## distributed equally across the observations that load on each random
## effect (via the Zt incidence) so that G is genuinely observation-
## level. Aggregating G by any clustering and forming the sandwich gives
## that clustering's cluster-robust meat; the crossed covariance combines
## several clusterings by Cameron-Gelbach-Miller inclusion-exclusion.
## Diagonal V_b only (block-diagonal / crossed random slopes: deferred).
.scoreByObs <- function(fit) {
    stopifnot(is(fit, "rlmerMod"))
    block_sizes <- vapply(fit@idx, nrow, integer(1))
    if (!all(block_sizes == 1L))
        stop(".scoreByObs (crossed designs) supports diagonal random ",
             "effects only; crossed random slopes are not yet covered.")
    eta <- fit@resp$eta
    pp  <- fit@pp
    p   <- pp$p; q <- pp$q; n <- pp$n
    if (length(eta) != n) eta <- rep.int(1, n)
    L   <- length(getME(fit, "theta"))
    rho_e  <- fit@rho.e; rho_b <- fit@rho.b; rho_sb <- fit@rho.sigma.b
    sigma  <- .sigma(fit)
    u_std  <- pp$b.s / sigma

    U_e <- pp$U_e; X <- pp$X; Z <- t(pp$Zt)
    r     <- as.numeric(fit@resp$y - fit@resp$offset -
                        X %*% .fixef(fit) - Z %*% (pp$U_b %*% pp$b.s))
    r_std <- as.numeric(U_e %*% r) / sigma
    psi_e <- rho_e@psi(r_std)
    U_eX  <- as.matrix(U_e %*% X)              # n x p
    ZL    <- as.matrix(U_e %*% Z %*% pp$U_b)   # n x q

    ## incidence / distribution weights: W[r, i] = 1/n_r if obs i loads
    ## on random effect r (n_r = #obs loading on r), 0 otherwise.
    Zt_t <- as(pp$Zt, "TsparseMatrix")
    ri <- Zt_t@i + 1L; ci <- Zt_t@j + 1L
    n_r <- tabulate(ri, nbins = q)
    W  <- sparseMatrix(i = ri, j = ci, x = 1 / n_r[ri], dims = c(q, n))

    G <- matrix(0, p + q + 1L + L, n)
    ## F_beta and F_u e-side rows (per observation)
    G[seq_len(p), ]       <- t(U_eX * (psi_e * eta))
    G[p + seq_len(q), ]   <- t(ZL   * (psi_e * eta))
    ## F_u penalty (per RE) distributed across loading observations
    lambda_e   <- rho_e@EDpsi()
    lambda_b_q <- pp$D_b@x
    psi_b      <- .applyBlockPsi_b(fit, u_std, rho_b)
    pen        <- -(lambda_e / lambda_b_q) * psi_b      # length q
    G[p + seq_len(q), ] <- G[p + seq_len(q), ] +
        as.matrix(Diagonal(x = pen) %*% W)
    ## F_sigma row (per observation)
    tau_e   <- pp$tau_e()
    arg_s   <- r_std / tau_e
    w_s     <- fit@rho.sigma.e@wgt(arg_s)
    G[p + q + 1L, ] <- eta * tau_e^2 * w_s * (arg_s^2 - pp$kappa_e)
    ## F_theta rows (per RE) distributed across loading observations
    tau_b  <- .getTauB_diag(fit, pp, rho_b, rho_sb)
    arg_th <- u_std / tau_b
    kappa_per_block <- pp$kappa_b
    for (l in seq_along(rho_sb)) {
        is_in <- which(fit@ind == l)
        w_th  <- rho_sb[[l]]@wgt(arg_th[is_in])
        vvec  <- numeric(q)
        vvec[is_in] <- w_th * u_std[is_in]^2 -
            kappa_per_block[[l]] * w_th * tau_b[is_in]^2
        G[p + q + 1L + l, ] <- as.numeric(crossprod(vvec, W))
    }
    rownames(G) <- c(paste0("beta", seq_len(p)), paste0("u", seq_len(q)),
                     "log_sigma", paste0("theta", seq_len(L)))
    G
}

## solve(-Jpar) %*% M with a pseudo-inverse fallback at a boundary.
.applyJparInv <- function(IF, M) {
    tryCatch(-solve(IF$Jpar, M), error = function(e) {
        sv  <- svd(IF$Jpar)
        tol <- max(dim(IF$Jpar)) * .Machine$double.eps * sv$d[1L]
        pos <- sv$d > tol
        dinv <- numeric(length(sv$d)); dinv[pos] <- 1 / sv$d[pos]
        -(sv$v %*% (dinv * t(sv$u))) %*% M
    })
}

## Cameron-Gelbach-Miller multiway cluster-robust covariance of the full
## parameter vector for a CROSSED design. With per-observation influence
## columns Z = -Jpar^{-1} G (G from .scoreByObs), the covariance is the
## inclusion-exclusion sum over non-empty subsets T of the grouping
## factors: A = sum_T (-1)^{|T|+1} A_{C(T)}, where C(T) clusters by the
## interaction of the factors in T and A_C aggregates the influence
## columns by C before the outer product. For two crossed factors this is
## the familiar A_subj + A_item - A_{subj:item}.
.vcov_varpar_crossed_Apar <- function(fit, IF, correction) {
    G  <- .scoreByObs(fit)
    ZI <- .applyJparInv(IF, G)              # d x n per-observation influence
    fl <- fit@flist
    K  <- length(fl)
    d  <- nrow(ZI)
    A_par <- matrix(0, d, d)
    for (m in seq_len(K)) {                  # subset size
        combs <- utils::combn(K, m)
        sgn   <- (-1)^(m + 1L)
        for (cc in seq_len(ncol(combs))) {
            sub <- combs[, cc]
            fac <- if (length(sub) == 1L) fl[[sub]]
                   else do.call(interaction, c(fl[sub], list(drop = TRUE)))
            ci  <- as.integer(fac); nc <- nlevels(fac)
            Xi_c <- rowsum(t(ZI), ci)        # nc x d (aggregate columns)
            A_c  <- crossprod(Xi_c)          # sum_c xi_c xi_c'
            if (correction == "G1" && nc > 1L) A_c <- A_c * nc / (nc - 1)
            A_par <- A_par + sgn * A_c
        }
    }
    list(A_par = A_par, J = min(vapply(fl, nlevels, integer(1))))
}

## IF-based covariance of (sigma-hat, theta-hat): the (sigma, theta)
## block of the joint cluster sandwich (see file header). Pass a
## pre-computed implicitIF_full(fit) via IF to reuse its Jacobian.
## Single-factor / nested designs use the one-way cluster sandwich
## (.scoreByCluster); crossed designs use the Cameron-Gelbach-Miller
## multiway covariance (.vcov_varpar_crossed_Apar).
.vcov_varpar <- function(fit, IF = NULL, correction = c("G1", "none")) {
    correction <- match.arg(correction)
    if (is.null(IF)) IF <- implicitIF_full(fit)
    ## WS14: variance components estimated at their lower bound (~0) sit on
    ## the boundary of the parameter space. Conditional on those being 0,
    ## the remaining (sigma, interior-theta) covariance is regular, and the
    ## Satterthwaite df computed on it is exactly the df of the reduced
    ## model (the model with those components dropped). We therefore zero
    ## the boundary components' rows/cols of A explicitly and, when the
    ## remaining block is non-degenerate ("reducible"), still return a
    ## usable covariance; only a genuinely singular fit (sigma or an
    ## interior theta non-identifiable) is unreliable.
    th0   <- getME(fit, "theta")
    lower <- getME(fit, "lower")
    bcomp <- which(lower == 0 & th0 < 1e-8)          # boundary theta indices
    singular <- isTRUE(IF$solve_info$singular)
    boundary <- singular || length(bcomp) > 0L
    crossed <- is.null(.top_cluster(fit))
    if (crossed) {
        cr    <- .vcov_varpar_crossed_Apar(fit, IF, correction)
        A_par <- cr$A_par
        J     <- cr$J
    } else {
        S  <- .scoreByCluster(fit)
        J  <- ncol(S)
        Xi <- .applyJparInv(IF, S)
        A_par <- tcrossprod(Xi)
        if (correction == "G1" && J > 1L) A_par <- A_par * J / (J - 1)
    }

    p <- fit@pp$p; q <- fit@pp$q
    L <- length(getME(fit, "theta"))
    sel <- c(p + q + 1L, p + q + 1L + seq_len(L))
    A <- A_par[sel, sel, drop = FALSE]
    ## log-sigma -> sigma
    sig0 <- .sigma(fit)
    A[1L, ] <- A[1L, ] * sig0
    A[, 1L] <- A[, 1L] * sig0
    A <- (A + t(A)) / 2
    ## The Cameron-Gelbach-Miller multiway estimator (crossed designs) is
    ## not guaranteed positive semi-definite in finite samples; project to
    ## the nearest PSD matrix by truncating negative eigenvalues (the
    ## standard CGM fix), so the variances and Satterthwaite df stay well
    ## defined. The one-way sandwich (single/nested) is already PSD.
    psd_corrected <- FALSE
    if (crossed) {
        ee <- eigen(A, symmetric = TRUE)
        if (any(ee$values < 0)) {
            A <- ee$vectors %*% (pmax(ee$values, 0) * t(ee$vectors))
            A <- (A + t(A)) / 2
            psd_corrected <- TRUE
        }
    }
    ## WS14: explicitly zero the boundary components' rows/cols (do not
    ## rely on the pseudo-inverse having zeroed them to machine precision),
    ## so the df denominator g' A g excludes them cleanly. The covariance
    ## is then "reducible" -- i.e. equals the reduced model's and yields a
    ## valid conditional df -- iff the remaining (sigma + interior-theta)
    ## block is positive definite. sigma is row/col 1; theta_k is row/col
    ## 1 + k.
    if (length(bcomp) > 0L) {
        brows <- 1L + bcomp
        A[brows, ] <- 0
        A[, brows] <- 0
    }
    ## "Reducible" iff, after dropping the boundary components, the
    ## remaining (sigma + interior-theta) block is positive definite. This
    ## is the right test on its own: the boundary always makes Jpar
    ## rank-deficient (so IF$solve_info$singular is uninformative here),
    ## but ANY *additional* non-identifiability -- a non-identifiable sigma
    ## or interior theta -- leaves a (near-)zero eigenvalue in this reduced
    ## block, which fails the test.
    keep   <- setdiff(seq_len(1L + L), 1L + bcomp)    # sigma + interior theta
    Asub   <- A[keep, keep, drop = FALSE]
    emin   <- if (nrow(Asub) > 0L)
        min(eigen(Asub, symmetric = TRUE, only.values = TRUE)$values)
        else Inf
    reducible <- emin > 1e-10 * max(1, max(abs(diag(Asub))))
    if (boundary && !reducible)
        warning("variance-component boundary (singular fit): the ",
                "(sigma, theta) covariance is not identifiable at ",
                "first order; returned values are unreliable. See the ",
                "boundary discussion in the IF-thread1 papers.")
    nm <- c("sigma", paste0("theta", seq_len(L)))
    dimnames(A) <- list(nm, nm)
    attr(A, "n.clusters")          <- J
    attr(A, "boundary")            <- boundary
    attr(A, "boundary.components") <- bcomp
    attr(A, "reducible")           <- reducible
    attr(A, "multiway")            <- crossed
    attr(A, "psd.corrected")       <- psd_corrected
    A
}

## WS16 (df default): a deterministic, machine-independent measure of the
## work implicitIF_full(fit) costs, used by summary(df = "auto") to decide
## whether to compute the Satterthwaite df by default. It is a pure
## function of the problem dimensions (n, p, q, the number of variance
## parameters L and the df method) -- NOT a wall-clock estimate and NOT a
## trial run -- so the same fit takes the same auto-decision on every
## machine. (A runtime probe would be more accurate but itself runs the
## costly tau machinery once, ~0.2-1.3 s on large fits, defeating the
## "summary stays fast" goal exactly when we mean to skip.)
##
## "Work" counts the dominant operation of implicitIF_full: numerically
## differentiating the RSE score, one O(n) score evaluation per Jacobian
## column.
##  - DASvar evaluates all d = p + q + 1 + L columns directly, each O(n),
##    so its work is d * n.
##  - DAStau reuses the converged tau on the p + q + 1 non-theta columns
##    (WS22) and re-runs the (tau | s) alternation only on the L theta
##    columns; profiling (Rprof on the reference machine) shows one such
##    calcTau column costs about 40 plain columns -- a ratio of the
##    algorithms (8.0e-3 vs 2.0e-4 s/column), not an absolute speed -- so
##    its work is 40 * L * n.
## Both methods cross the default cutoff (5000) at the same problem size
## the old one-second wall-clock budget did on the reference machine:
## about n = 170 for a single-intercept DASvar fit and n = 125 for DAStau.
## The cutoff is options(robustlmm.summary.df.max = <work units>), default
## 5000; it is dimensionless and reproducible, so it needs no per-machine
## calibration.
.satterthwaite_df_workload <- function(fit) {
    pp <- fit@pp
    n  <- as.double(pp$n); p <- pp$p; q <- pp$q
    L  <- length(getME(fit, "theta"))
    if (identical(fit@method, "DASvar"))
        (p + q + 1L + L) * n               # d direct FD columns, each O(n)
    else                                   # DAStau (and any iterative tau)
        40 * max(L, 1L) * n                # L calcTau columns, ~40x heavier
}

## Should summary() compute the Satterthwaite df by default for this fit?
## Yes if the IF is already cached (then it is essentially free, WS16
## caching), or if the deterministic workload is within the cutoff. The
## workload is returned via attr(., "work") so the caller can phrase the
## skip note.
.satterthwaite_df_default_ok <- function(fit) {
    if (!is.null(.IFfullCached(fit)))
        return(structure(TRUE, work = 0, cached = TRUE))
    cutoff <- getOption("robustlmm.summary.df.max", 5000)
    work   <- .satterthwaite_df_workload(fit)
    structure(work <= cutoff, work = work, cached = FALSE)
}

## Satterthwaite degrees of freedom for fixed-effects contrasts
## (PLAN-WS16 §2-3, step 1 second half):
##     df(c) = 2 v(c)^2 / (g' A g),
##     v(c)  = c' V(sigma, theta) c,   V = sigma^2 unsc(theta),
##     g     = d v / d (sigma, theta),
##     A     = .vcov_varpar(fit)  (IF-based cluster sandwich).
## sigma enters v multiplicatively, so the sigma-coordinate of g is
## analytic (2 v / sigma); the theta coordinates are central FD over
## pp$setTheta + unsc() recomputes (setTheta invalidates the unsc/M
## caches; D_b / Lambda_b are rho-only, theta-independent). The fit's
## pp is restored on exit.
##
## Lmat: contrasts as rows (default: one per fixed effect). Only the
## default vcov is supported (the sandwich version needs vcov_sandwich
## recomputed at perturbed theta -- score-dependent, deferred).
.satterthwaite_df <- function(fit, Lmat = NULL, IF = NULL,
                              h_rel = 1e-5) {
    stopifnot(is(fit, "rlmerMod"))
    p <- fit@pp$p
    if (is.null(Lmat)) {
        Lmat <- diag(p)
        rownames(Lmat) <- names(.fixef(fit))
    }
    Lmat <- rbind(Lmat)
    stopifnot(ncol(Lmat) == p)

    A    <- .vcov_varpar(fit, IF)
    pp   <- fit@pp
    sig0 <- .sigma(fit)
    th0  <- getME(fit, "theta")
    L    <- length(th0)

    on.exit(pp$setTheta(th0))
    v_at <- function(theta) {
        pp$setTheta(theta)
        V <- sig0^2 * pp$unsc()
        rowSums((Lmat %*% V) * Lmat)   # c_j' V c_j per row
    }

    ## WS14: theta components at the lower bound (~0) are held at the
    ## boundary -- their A row/col is zeroed, so they are excluded from the
    ## df. Skip the central finite difference for them: it would step theta
    ## below 0 (an infeasible negative variance) into a possibly invalid
    ## unsc() evaluation. Their gradient is irrelevant (multiplied by 0).
    bcomp <- attr(A, "boundary.components")
    if (is.null(bcomp)) bcomp <- integer(0)

    v0 <- v_at(th0)
    g  <- matrix(0, 1L + L, nrow(Lmat))
    g[1L, ] <- 2 * v0 / sig0
    for (l in setdiff(seq_len(L), bcomp)) {
        h  <- h_rel * max(abs(th0[l]), 0.1)
        tp <- th0; tp[l] <- tp[l] + h
        tm <- th0; tm[l] <- tm[l] - h
        g[1L + l, ] <- (v_at(tp) - v_at(tm)) / (2 * h)
    }

    df <- vapply(seq_len(nrow(Lmat)), function(j) {
        denom <- drop(crossprod(g[, j], A %*% g[, j]))
        if (denom <= 0) return(Inf)
        2 * v0[j]^2 / denom
    }, numeric(1))
    names(df) <- rownames(Lmat)
    attr(df, "v")          <- v0
    attr(df, "se")         <- sqrt(v0)
    attr(df, "boundary")   <- attr(A, "boundary")
    attr(df, "reducible")  <- attr(A, "reducible")
    attr(df, "n.boundary") <- length(bcomp)
    df
}

## Denominator degrees of freedom for a multi-row Wald contrast via the
## Satterthwaite F approximation (the multivariate generalisation used
## by lmerTest::contestMD). For a contrast matrix L (rows = the q
## restrictions of a hypothesis L beta = 0), decompose the contrast
## covariance L V L' into its rank positive eigen-directions, take the
## per-direction Satterthwaite df nu_m from .satterthwaite_df, and
## combine:
##     E   = sum_{nu_m > 2} nu_m / (nu_m - 2),
##     df2 = 2 E / (E - q')     (when E > q').
## Returns df1 = q' (numerator df = rank of the contrast) and df2 (the
## F denominator df; NA when the combination is undefined, e.g. every
## direction has nu_m <= 2). Default vcov only (matches .satterthwaite_df).
.satterthwaite_Fdenom <- function(fit, Lmat, IF = NULL) {
    Lmat <- rbind(Lmat)
    V  <- as.matrix(vcov(fit, type = "default"))
    VL <- tcrossprod(Lmat %*% V, Lmat)
    eig <- eigen(VL, symmetric = TRUE)
    tol <- sqrt(.Machine$double.eps) * max(c(eig$values, 0))
    pos <- eig$values > tol
    q   <- sum(pos)
    if (q == 0L) return(list(df1 = 0L, df2 = NA_real_))
    ## direction contrasts in coefficient space: row m = v_m' L
    dirs <- crossprod(eig$vectors[, pos, drop = FALSE], Lmat)  # q x p
    nu <- as.numeric(.satterthwaite_df(fit, Lmat = dirs, IF = IF))
    nu <- nu[is.finite(nu) & nu > 2]
    if (length(nu) == 0L) return(list(df1 = q, df2 = NA_real_))
    E   <- sum(nu / (nu - 2))
    df2 <- if (E > q) 2 * E / (E - q) else NA_real_
    list(df1 = q, df2 = df2)
}
