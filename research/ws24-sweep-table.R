## Assemble the per-version sweep results into one table + rds.
##   Rscript ws24-sweep-table.R <outdir>
outdir <- commandArgs(trailingOnly = TRUE)[1]
ord <- c("3.1-2","3.2-0","3.2-2","3.2-3","3.2-5","3.3-1","3.3-2","3.3-3",
         "3.4-2","3.4-3","3.4-5","3.5.0-1","dev-prefix","dev-fixed")
rows <- lapply(ord, function(v) {
    f <- file.path(outdir, paste0("sweep-", v, ".rds"))
    if (!file.exists(f)) return(NULL)
    r <- readRDS(f)
    err <- function(x) if (is.na(x["sigma"])) NA_real_ else x["sigma"]/r$lmer["sigma"] - 1
    data.frame(version = v, robustlmm = r$robustlmm,
               lmer_sigma = unname(r$lmer["sigma"]),
               DASvar_sigma = unname(r$DASvar["sigma"]),
               DASvar_relerr = unname(err(r$DASvar)),
               DAStau_sigma = unname(r$DAStau["sigma"]),
               DAStau_relerr = unname(err(r$DAStau)),
               DASvar_se = unname(r$DASvar["se"]), lmer_se = unname(r$lmer["se"]),
               note = if (any(grepl("^ERROR", c(r$DASvar_msgs, r$DAStau_msgs))))
                          sub("ERROR: ", "", grep("^ERROR", r$DASvar_msgs, value = TRUE)[1])
                      else "", stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab$stamp_R <- readRDS(file.path(outdir, "sweep-3.4-5.rds"))$Rversion
tab$stamp_lme4 <- readRDS(file.path(outdir, "sweep-3.4-5.rds"))$lme4
saveRDS(tab, file.path(outdir, "ws24-version-sweep.rds"))
cat("Check dataset: J = 50 x 5 (n = 250), weights log-uniform in [0.255, 3.99],\n",
    "model y ~ x + (1 | g), rho.e = rho.b = cPsi -- at which the RSE reduces\n",
    "to REML, so sigma-hat must equal lmer's exactly.\n\n", sep = "")
p <- tab[, c("version","robustlmm","DASvar_sigma","DASvar_relerr",
             "DAStau_sigma","DAStau_relerr","note")]
names(p) <- c("version","pkg","DASvar sigma","rel.err","DAStau sigma","rel.err","note")
print(p, digits = 6, row.names = FALSE)
cat(sprintf("\nlmer reference sigma: %.6f (identical in every version)\n",
            unique(round(tab$lmer_sigma, 10))))
