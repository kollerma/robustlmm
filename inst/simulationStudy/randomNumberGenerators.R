########################################################
## Random number generators used in various scripts   ##
########################################################

if (exists("centerAndScale")) {
    ..backup..centerAndScale.. <- centerAndScale
}
if (exists(".Random.seed")) {
    ..backup..Random.seed.. <- .Random.seed
}

centerAndScale <- function(rand, doCenter = TRUE) {
    set.seed(1)
    ta <- rand(100000)
    tsum <- MASS::hubers(ta, k = 1.345)
    scale <- tsum$s
    center <- if (doCenter) {
        tsum$mu
    } else {
        0
    }
    if (doCenter) {
        return(function(n)
            (rand(n) - center) / scale)
    } else {
        return(function(n)
            rand(n) / scale)
    }
}

rcnorm <-
    function (n,
              mean = 0,
              sd = 1,
              epsilon = 0.1,
              meanc = mean,
              sdc = sqrt(10) * sd) {
        e <- rnorm(n, mean, sd)
        nc <- floor(epsilon * n)
        idx <- sample(1:n, nc)
        e[idx] <- rnorm(nc, meanc, sdc)
        e
    }

srnorm <- rnorm
srt3 <- centerAndScale(function(n)
    rt(n, df = 3), doCenter = FALSE)
srskt3 <-
    centerAndScale(function(n)
        skewt::rskt(n, df = 3, gamma = 2))
srcnorm <-
    centerAndScale(function(n)
        rcnorm(
            n,
            epsilon = 0.1,
            meanc = 4,
            sdc = 1
        ))

rm(centerAndScale)

if (exists("..backup..centerAndScale..")) {
    centerAndScale <- ..backup..centerAndScale..
    rm(..backup..centerAndScale..)
}
if (exists("..backup..Random.seed..")) {
    .Random.seed <- ..backup..Random.seed..
    rm(..backup..Random.seed..)
}
