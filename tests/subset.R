## Disabled on the CRAN release branch (historical github-transition set:
## slow / MC-random / platform-fragile tests kept out of the CRAN check;
## they run in full on master and in CI). See feedback in project notes.
quit()

## testing subset
require(robustlmm)

fit1 <- rlmer(Yield ~ (1 | Batch), Dyestuff, subset=Batch != "F")
fit2 <- rlmer(Yield ~ (1 | Batch), subset(Dyestuff, Batch != "F"))

stopifnot(all.equal(getInfo(fit1)[-1], getInfo(fit2)[-1]))
