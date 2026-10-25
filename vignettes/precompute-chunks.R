## Precompute the output of the slowest code chunks of
## vignettes/robustlmm-3.5.0.Rnw, to keep the CRAN check time down. In
## the vignette these chunks have eval=FALSE, echo=FALSE, and a hidden
## chunk writes the code and the stored output, as Sweave would, via
## showChunk(), which stops if the code of the chunk has changed since
## this script was run.
##
## The script runs all evaluated chunks of the vignette in order (so
## objects and the random number stream are as in a live build) and
## also evaluates the chunks listed in `cached`, recording their output.
##
## Re-run after changing a cached chunk or the code it calls:
##   Rscript vignettes/precompute-chunks.R
##
## Build-ignored (see .Rbuildignore).

rnw <- "vignettes/robustlmm-3.5.0.Rnw"
out <- "vignettes/robustlmm-3.5.0-chunks.rds"
cached <- c("confint-df", "cooks-obs", "ransac-init", "ransac-consensus",
            "structured-fits")

## the chunks of a Sweave file: label, options and code
parseChunks <- function(file) {
    src <- readLines(file)
    starts <- grep("^<<.*>>=[[:space:]]*$", src)
    lapply(starts, function(s) {
        e <- s + match(TRUE, grepl("^@", src[-seq_len(s)]))
        opts <- strsplit(sub("^<<(.*)>>=.*$", "\\1", src[s]), ",")[[1]]
        opts <- trimws(opts)
        named <- grepl("=", opts)
        o <- setNames(sub("^[^=]*=", "", opts[named]),
                      trimws(sub("=.*$", "", opts[named])))
        list(label = if (any(!named)) opts[!named][1] else "",
             opts = o,
             code = if (e > s + 1) src[(s + 1):(e - 1)] else character())
    })
}

isFalse <- function(x) !is.null(x) && !is.na(x) && toupper(x) %in% c("FALSE", "F")
isTrue <- function(x) !is.null(x) && !is.na(x) && toupper(x) %in% c("TRUE", "T")

## Sweave's evaluation: auto-print visible values; returns, for each
## expression, its source lines and its output
runChunk <- function(code) {
    exprs <- parse(text = code, keep.source = TRUE)
    lapply(seq_along(exprs), function(i) {
        o <- utils::capture.output({
            v <- withVisible(eval(exprs[[i]], globalenv()))
            if (v$visible) print(v$value)
        })
        list(src = as.character(attr(exprs, "srcref")[[i]]), output = o)
    })
}

chunks <- parseChunks(rnw)
stopifnot(all(cached %in% vapply(chunks, `[[`, "", "label")))

oldwd <- setwd(dirname(rnw))
res <- list()
for (ch in chunks) {
    isCached <- ch$label %in% cached
    if (isFalse(ch$opts["eval"]) && !isCached) next
    ## the hidden chunks that print the stored output
    if (ch$label %in% paste0(cached, "-output")) next
    if (isTrue(ch$opts["fig"])) pdf(NULL)
    o <- runChunk(ch$code)
    if (isTrue(ch$opts["fig"])) dev.off()
    if (isCached)
        res[[ch$label]] <- list(code = ch$code, steps = o)
}
setwd(oldwd)

attr(res, "robustlmm") <- as.character(packageVersion("robustlmm"))
saveRDS(res, out)
