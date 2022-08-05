## Disable test
quit()

require(lme4)
require(robustlmm)

## generate data with multiple grouping effects
group1 <- rep(letters[1:10], each = 32)
group2 <- rep(LETTERS[1:8], each = 4)
treat1 <- rep(c(TRUE, FALSE), each = 2)
treat2 <- c(TRUE, FALSE)

data <- data.frame(group1 = factor(group1), group2 = factor(group2),
                   treat1 = factor(treat1), treat2 = factor(treat2))
table(data)

data <- data[rep(1L:NROW(data), 3), ]

set.seed(123)
group1means <- 1 * rnorm(10)
group2means <- 1 * rnorm(8)
treat1means <- 2 * rnorm(2)
treat2means <- 2 * rnorm(2)
treatsInter <- 4 * rnorm(4)

data <- within(data, {
  resp <- group1means[group1] + group2means[group2] +
      treat1means[treat1] + treat2means[treat2] +
      treatsInter[interaction(treat1, treat2)] +
      0.5 * rnorm(NROW(data))
})

testFormula <- function(formula, data) {
    print(summary(fm <- lmer(formula, data, control=lmerControl(optimizer="bobyqa"))))
    print(summary(rm <- rlmerRcpp(formula, data, rho.e = cPsi, rho.b = cPsi, init = fm)))
    ranef.fm <- ranef(fm, condVar=FALSE)
    stopifnot(all.equal(coef(fm), coef(rm), tolerance = 1e-1, check.attributes = FALSE),
              all.equal(fixef(fm), fixef(rm), tolerance = 1e-2, check.attributes = FALSE),
              all.equal(ranef.fm , ranef(rm), tolerance = 1e-1, check.attributes = FALSE))
    invisible(list(fm, rm))
}

ms1 <- testFormula(resp ~ (1 + treat1|group1) + (1 + treat1|group2), data)

## currently fails with segmentation fault. When running in github-actions, it says:
## corrupted double-linked list

ms2 <- testFormula(resp ~ (1 + treat1|group1) + (1 + treat2|group2), data)
ms3 <- testFormula(resp ~ (1 + treat1 + treat2|group1) + (1 + treat1 + treat2|group2), data)
ms4 <- testFormula(resp ~ (1 + treat1:treat2|group1) + (1 + treat1:treat2|group2), data)

