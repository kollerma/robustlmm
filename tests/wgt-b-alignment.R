## Disabled on the CRAN release branch to keep the overall check time
## under the 10-minute CRAN limit; runs in full on master and in CI.
quit()

## Alignment of the random-effect robustness weights with the spherical
## random effects b_s, and the labeling that is built on top of it.
##
## Background (bugfix): lme4 >= 2.0-0 expands a heterogeneous diag() term
## into one 1x1 block type PER COORDINATE (findBlocks), so object@idx
## holds interleaved b-indices (type 1 = b 1,3,5,..., type 2 = b 2,4,6,...).
## wgt.b() used to concatenate the weights per block type, silently
## permuting the vector out of b order for such fits; every consumer that
## indexes by b (fitEffects, getME(., "w_b"), qq plots,
## .re_min_weight_by_group and hence the anova bootstrap guard /
## null = "robust" trimming) then used the wrong weights or labels.
## wgt.b() now assigns each weight back through its b-index, so
## getME(., "w_b_vector") is aligned with getME(., "b_s").
##
## Checks:
## 1. Heterogeneous diag() fit: the block types really are interleaved
##    (the shape that used to trigger the bug), and getME(., "w_b_vector")
##    matches an independent reconstruction (weight of b_s[j] computed
##    directly from dist.b and the block type owning j).
## 2. The b layout is level-major within each term (verified against the
##    Z column names, which lme4 labels with the grouping levels), and
##    .re_min_weight_by_group's labeled minima match an independent
##    per-level reconstruction.
## 3. getME(., "w_b") carries the correct level names on the rows and the
##    correct per-(level, coordinate) values; same for "w_sigma_b_vector"
##    alignment.
## 4. Simple (1 | g) fit: same checks (regression guard for the
##    single-block-type shape that was already correct).

suppressMessages(require(robustlmm))

## weight of b_s[j], computed independently of wgt.b()'s loop order:
## look up the block type owning j via object@idx and apply its wgt
indepWgtB <- function(fit, use.rho.sigma = FALSE) {
    db  <- robustlmm:::dist.b(fit, center = use.rho.sigma)
    rho <- robustlmm:::rho.b(fit, if (use.rho.sigma) "sigma" else "default")
    w   <- rep(NA_real_, length(db))
    for (bt in seq_along(fit@blocks)) {
        bind <- as.vector(fit@idx[[bt]])
        w[bind] <- rho[[bt]]@wgt(db[bind])
    }
    stopifnot(!anyNA(w))
    w
}

## independent per-level minima using the level-major Gp layout, after
## that layout has itself been verified against the Z column names
checkFit <- function(fit) {
    w_pkg <- getME(fit, "w_b_vector")
    w_ind <- indepWgtB(fit)
    stopifnot(all.equal(as.numeric(w_pkg), w_ind, tolerance = 1e-12))

    cnms <- fit@cnms; fl <- fit@flist; asg <- attr(fl, "assign")
    Gp <- fit@Gp
    zn <- dimnames(robustlmm:::getZ(fit))[[2]]
    lab <- character(length(w_ind))
    for (i in seq_along(cnms)) {
        nc <- length(cnms[[i]]); levs <- levels(fl[[asg[i]]])
        ## level-major layout: Z columns of term i are named
        ## rep(levels, each = nc)
        stopifnot(identical(zn[Gp[i] + seq_len(nc * length(levs))],
                            rep(levs, each = nc)))
        lab[Gp[i] + seq_len(nc * length(levs))] <-
            paste0(names(fl)[asg[i]], ": ", rep(levs, each = nc))
    }
    mw <- robustlmm:::.re_min_weight_by_group(fit)
    stopifnot(isTRUE(attr(mw, "labeled")))
    ref <- vapply(split(w_ind, factor(lab, levels = unique(lab))),
                  min, numeric(1))
    stopifnot(identical(names(mw), names(ref)),
              all.equal(as.numeric(mw), as.numeric(ref),
                        tolerance = 1e-12))

    ## getME(., "w_b"): per-factor data frames, rows = levels,
    ## columns = coefficients, values from the b-aligned vector
    wb <- getME(fit, "w_b")
    for (id in names(wb)) {
        df <- wb[[id]]
        terms <- which(names(cnms) == id)
        stopifnot(identical(colnames(df),
                            as.vector(unlist(cnms[terms]))),
                  identical(rownames(df), levels(fl[[id]])))
        col <- 0L
        for (i in terms) {
            nc <- length(cnms[[i]]); nl <- length(levels(fl[[asg[i]]]))
            m <- matrix(w_ind[Gp[i] + seq_len(nc * nl)], nrow = nc)
            for (k in seq_len(nc)) {
                col <- col + 1L
                stopifnot(all.equal(as.numeric(df[[col]]), m[k, ],
                                    tolerance = 1e-12))
            }
        }
    }

    ## sigma weights use the same alignment
    stopifnot(all.equal(as.numeric(getME(fit, "w_sigma_b_vector")),
                        indepWgtB(fit, use.rho.sigma = TRUE),
                        tolerance = 1e-12))
    invisible(TRUE)
}

data(sleepstudy, package = "lme4")
set.seed(17)
## contaminate two subjects so the weights are spread out and a
## permutation cannot cancel out
sleepstudy$y2 <- sleepstudy$Reaction +
    ifelse(sleepstudy$Subject %in% c("308", "309"), 120, 0)

## ---- 1.-3. heterogeneous diag() term ---------------------------------
fitDiag <- rlmer(y2 ~ Days + diag(Days | Subject), sleepstudy,
                 method = "DASvar")
## the shape that used to break: several interleaved 1x1 block types
stopifnot(length(fitDiag@blocks) == 2L,
          identical(as.vector(fitDiag@idx[[1L]]),
                    as.integer(seq(1L, 35L, by = 2L))),
          identical(as.vector(fitDiag@idx[[2L]]),
                    as.integer(seq(2L, 36L, by = 2L))))
checkFit(fitDiag)
## the fit really downweights the contaminated subjects and their labels
## point at them
mwDiag <- robustlmm:::.re_min_weight_by_group(fitDiag)
stopifnot(all(c("Subject: 308", "Subject: 309") %in%
              names(mwDiag)[mwDiag < 1]))
cat("wgt-b-alignment: diag() term ok\n")

## ---- 4. simple (1 | g) term ------------------------------------------
fitInt <- rlmer(y2 ~ Days + (1 | Subject), sleepstudy, method = "DASvar")
stopifnot(length(fitInt@blocks) == 1L)
checkFit(fitInt)
cat("wgt-b-alignment: (1 | g) term ok\n")

cat("wgt-b-alignment: all tests passed\n")
