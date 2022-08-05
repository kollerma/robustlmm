##' This is a convenience function to make it simple to access the simulation
##' study script files that are shipped with robustlmm.
##'
##' The function creates a copy of the script file that can be safely edited
##' without changing the original file.
##' @title Access Simulation Study Code
##' @param study Name of the script file, partial matching is supported via
##'   \code{\link{match.arg}}.
##' @param destinationPath optional path to directory in which the copy of the
##'   script should be created. By default the current working directory is
##'   used.
##' @param overwrite logical; should existing destination files be overwritten?
##' @examples
##' \dontrun{
##'   viewCopyOfSimulationStudy("sensitivityCurves")
##' }
##' @export
viewCopyOfSimulationStudy <-
    function(study = c(
        "sensitivityCurves.R",
        "consistencyAndEfficiencyDiagonal.R",
        "consistencyAndEfficiencyBlockDiagonal.R",
        "breakdown.R",
        "convergence.R",
        "robustnessDiagonal.R",
        "robustnessBlockDiagonal.R"
    ),
    destinationPath = getwd(),
    overwrite = FALSE) {
        study <- match.arg(study)
        original <- getOriginalFile(study)
        destination <- file.path(destinationPath, study)
        file.copy(original, destination, overwrite = overwrite)
        utils::file.edit(destination)
        return(destination)
    }

getOriginalFile <- function(study) {
    return(system.file(file.path("simulationStudy", study),
                       package = "robustlmm"))
}
