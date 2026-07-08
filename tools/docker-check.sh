#!/bin/bash
# Run the CI-equivalent R CMD check inside a rocker container.
# Usage (from the package root, tarball already built with
#   R CMD build --no-manual --compact-vignettes=both .):
#   docker run --rm -v "$PWD":/pkg -w /tmp rocker/r-ver:<tag> /pkg/tools/docker-check.sh
# Mirrors .github/workflows/check-standard.yaml:
#   R CMD check --no-manual --as-cran --no-stop-on-test-error
set -uo pipefail

## rocker/drd ships R-devel as RD (default R is the release build)
RBIN=$(command -v RD || command -v R)
RSCRIPT=$(command -v RDscript || command -v Rscript)

echo "== $($RBIN --version | head -1) on $(uname -m) =="

apt-get update -qq > /dev/null 2>&1 || true
## qpdf: --as-cran PDF checks; cmake: nloptr source build; lapack/blas/
## gfortran: linking on images without them (drd); texlive: the Sweave
## vignette rebuild needs pdflatex (CI used tinytex + ae adjustbox
## collectbox)
apt-get install -y -qq qpdf cmake liblapack-dev libblas-dev gfortran \
    texlive-latex-base texlive-latex-recommended texlive-latex-extra \
    texlive-fonts-recommended > /dev/null 2>&1 || true

$RSCRIPT -e '
## current packages, not the image-pinned P3M snapshot: DESCRIPTION
## needs lme4 >= 2.0-1, newer than the oldrel image snapshot carries.
## P3M latest serves Ubuntu arm64 binaries where available; cloud CRAN
## is the source fallback (cmake installed above for nloptr).
options(repos = c(P3M = "https://p3m.dev/cran/__linux__/noble/latest",
                  CRAN = "https://cloud.r-project.org"))
ncpus <- max(1L, parallel::detectCores() - 1L)
deps <- c("lme4", "Matrix", "robustbase", "lattice", "nlme", "xtable",
          "Rcpp", "fastGHQuad", "numDeriv", "rlang", "reformulas",
          "ggplot2", "reshape2", "microbenchmark", "emmeans",
          "estimability", "lqmm", "MASS", "lemon", "RColorBrewer",
          "skewt", "fs", "dplyr", "ggh4x", "testthat", "knitr",
          "robustvarComp", "confintROB")
have <- function() {
    ip <- installed.packages()
    ok <- deps %in% rownames(ip)
    ok[deps == "lme4"] <- ok[deps == "lme4"] &&
        "lme4" %in% rownames(ip) &&
        package_version(ip["lme4", "Version"]) >= "2.0-1"
    deps[!ok]
}
miss <- have()
if (length(miss)) install.packages(miss, Ncpus = ncpus, quiet = TRUE)
miss <- have()
optional <- c("lqmm", "skewt", "ggh4x", "robustvarComp", "confintROB",
              "microbenchmark")
hard <- setdiff(miss, optional)
if (length(hard)) stop("missing hard deps: ", paste(hard, collapse = ", "))
if (length(miss)) message("missing optional: ", paste(miss, collapse = ", "))
'
[ $? -ne 0 ] && exit 1

cd /tmp
cp /pkg/robustlmm_*.tar.gz /tmp/
NOT_CRAN=true _R_CHECK_CRAN_INCOMING_=false _R_CHECK_CRAN_INCOMING_REMOTE_=false \
_R_CHECK_FORCE_SUGGESTS_=false $RBIN CMD check --no-manual --as-cran \
    --no-stop-on-test-error /tmp/robustlmm_*.tar.gz
status=$?
echo "== check exit: $status =="
dest="/pkg/docker-check-results/$($RBIN --version | head -1 | tr -c '[:alnum:].\n' '_')"
mkdir -p "$dest"
cp /tmp/robustlmm.Rcheck/00check.log /tmp/robustlmm.Rcheck/00install.out \
   "$dest/" 2>/dev/null || true
tail -20 /tmp/robustlmm.Rcheck/00check.log
exit $status
