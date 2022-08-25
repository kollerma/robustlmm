% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processFile.R
\name{mergeProcessedFits}
\alias{mergeProcessedFits}
\title{Merge Processed Fits}
\usage{
mergeProcessedFits(processedFitList)
}
\arguments{
\item{processedFitList}{list of processed fits as produced by
\code{\link{processFit}}.}
}
\value{
similar list as returned by \code{\link{processFit}} just with
  matrix entries instead of vectors.
}
\description{
Combine list of processed fits into one list in matrix form.
}
\examples{
  z1 <- rep(1, 10)
  z2 <- 0:9
  K <- list()
  K[[1]] <- tcrossprod(z1, z1) ## Int
  K[[2]] <- tcrossprod(z1, z2) + tcrossprod(z2, z1) ## Int:Days
  K[[3]] <- tcrossprod(z2, z2)
  names(K) <- c("Subject.Intercept.", "Subject.Days.Intercept.",
                "Subject.Days")
  p <- length(z2)
  n <- length(unique(sleepstudy$Subject))
  groups <- cbind(rep(1:p, each = n), rep(1:n, p))

  preparedDataset <-
      prepareMixedEffectDataset(Reaction ~ Days + (Days|Subject),
                                sleepstudy,
                                lmeFormula = Reaction ~ Days,
                                heavyLmeRandom = ~ Days,
                                heavyLmeGroups = ~ Subject,
                                lqmmRandom = ~ Days,
                                lqmmGroup = "Subject",
                                lqmmCovariance = "pdSymm",
                                groups = groups,
                                varcov = K,
                                lower = c(0, -Inf, 0))
  set.seed(1)
  datasets <- generateMixedEffectDatasets(2, preparedDataset)

  fits <- fitDatasets_lmer(datasets)
  if (require(heavy)) {
      fits <- c(fits, fitDatasets_heavyLme(datasets))
  }
  if (require(lqmm)) {
      fits <- c(fits, fitDatasets_lqmm(datasets))
  }
  if (require(robustvarComp)) {
      fits <- c(fits, fitDatasets_varComprob_compositeTau(datasets))
  }
  processedFits <- lapply(fits, processFit, all = TRUE)
  merged <- mergeProcessedFits(processedFits)
  str(merged)
}