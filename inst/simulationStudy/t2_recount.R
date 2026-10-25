## Shipped with robustlmm: a script of the studies for the small-sample and
## variance-component methods in vignette("simulationStudies"). Run it with
## Rscript from a working directory of your choice; results are written to a
## subdirectory there. The numbers in the vignette come from runs with the
## development version of robustlmm 3.5.1.
##
## Recompute T2's warp-speed sizes from the cell files of
## t2_threshold_sim.R, both with the raw statistic (the rule fixed before
## the run) and with the tie rule that was adopted after the boundary cells
## had finished
## (D <= 0.01 treated as 0, observed and bootstrap alike). Size at 5% =
## mean(D_obs > q_0.95(D*)) with ties: mean(D_obs' > q) where D' = D, or 0.
## Run from the directory in which t2_threshold_sim.R was run.
tie <- function(D) ifelse(D <= 0.01, 0, D)
rows <- list()
for (f in sort(list.files("t2_results", "^cell_.*[.]rds$",
                          full.names = TRUE))) {
    o <- readRDS(f); r <- o$res
    r <- r[is.finite(r[, "D_obs"]) & is.finite(r[, "D_raw"]) &
           is.finite(r[, "D_thr"]), , drop = FALSE]
    size <- function(Dobs, Db) mean(Dobs > stats::quantile(Db, 0.95))
    rows[[length(rows) + 1L]] <- data.frame(
        cell = sub("cell_(.*)[.]rds", "\\1", basename(f)), n = nrow(r),
        raw = size(r[, "D_obs"], r[, "D_raw"]),
        thr = size(r[, "D_obs"], r[, "D_thr"]),
        raw_tied = size(tie(r[, "D_obs"]), tie(r[, "D_raw"])),
        thr_tied = size(tie(r[, "D_obs"]), tie(r[, "D_thr"])),
        q95_raw_tied = unname(stats::quantile(tie(r[, "D_raw"]), 0.95)),
        q95_thr_tied = unname(stats::quantile(tie(r[, "D_thr"]), 0.95)),
        atom_obs = mean(abs(r[, "D_obs"]) <= 0.01),
        neg_far = sum(c(r[, "D_obs"], r[, "D_raw"], r[, "D_thr"]) < -0.01),
        min_D = min(c(r[, "D_obs"], r[, "D_raw"], r[, "D_thr"])),
        fired = mean(r[, "trig"] > 0))
}
tab <- do.call(rbind, rows)
print(tab, digits = 3, row.names = FALSE)
## D10: thresholded size in [0.035, 0.065] in every cell, and no farther
## from 0.05 than raw where the nuisance is on the boundary (sx = 0)
d10 <- function(thr, raw, cell)
    all(thr >= 0.035 & thr <= 0.065) &&
        all((abs(thr - 0.05) <= abs(raw - 0.05))[grepl("sx0.0", cell)])
cat(sprintf("D10 on the pre-registered statistic: %s; with the tie rule: %s\n",
            if (d10(tab$thr, tab$raw, tab$cell)) "met" else "NOT met",
            if (d10(tab$thr_tied, tab$raw_tied, tab$cell)) "met" else "NOT met"))
saveRDS(tab, "t2_results/recount.rds")
