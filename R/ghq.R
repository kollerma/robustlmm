##' @importFrom fastGHQuad gaussHermiteData
ghq <- function(n = 1, modify = TRUE) {
    rule <- fastGHQuad::gaussHermiteData(n)
    return(list(nodes=rule$x,
                weights= if (modify) rule$w*exp(rule$x^2) else rule$w))
}
