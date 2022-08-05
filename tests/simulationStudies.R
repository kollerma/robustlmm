require(robustlmm)

## verify all scripts listed in argument of viewCopyOfSimulationStudy actually
## exist
sourceDirectory <-
    system.file("simulationStudy", package = "robustlmm")
allScripts <- list.files(sourceDirectory, pattern = ".R$")
actualScripts <- eval(formals(viewCopyOfSimulationStudy)$study)
stopifnot(all(actualScripts %in% allScripts))
