xtraR <- system.file("xtraR", package="robustlmm")
unitTests <- list.files(xtraR, pattern = "unitTests-", full.names = TRUE)
for (unitTest in unitTests) {
  cat(paste0("Running ", unitTest, "...\n"))
  source(unitTest)
}
