## Rebuild the pre-built vignettes/simulationStudies.pdf from
## vignettes/simulationStudies.Rnw. The PDF is shipped as is (engine
## R.rsp::asis, see simulationStudies.pdf.asis), so that R CMD check
## does not re-run the vignette; re-run this script after changing the
## Rnw or the results it reads:
##   Rscript vignettes/make-simulationStudies.R
##
## Build-ignored (see .Rbuildignore), as is the Rnw.

oldwd <- setwd("vignettes")
utils::Sweave("simulationStudies.Rnw")
tools::texi2pdf("simulationStudies.tex", clean = TRUE)
tools::compactPDF("simulationStudies.pdf", gs_quality = "ebook")
unlink(c("simulationStudies.tex", list.files(pattern = "^simulationStudies-")))
setwd(oldwd)
