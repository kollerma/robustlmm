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
group1means <- 4 * rnorm(10)
group2means <- 4 * rnorm(8)
treat1means <- 2 * rnorm(18)
treat2means <- 2 * rnorm(18)
treatsInter <- 2 * rnorm(80)
errors <- 0.5 * rnorm(NROW(data))

data <- within(data, {
  resp1 <- group1means[group1] + errors
  resp2 <- resp1 + group2means[group2]
  resp3 <- resp2 + as.logical(treat1) * (1 + treat1means[group1])
  resp4 <- resp3 + as.logical(treat2) * (1 + treat2means[10 + as.integer(group2)])
  resp5 <- resp4 + as.logical(treat2) * (1 + treat2means[group1])
  resp6 <- resp5 + as.logical(treat1) * (1 + treat1means[10 + as.integer(group2)])
  resp7 <- resp6 + as.logical(treat1) * as.logical(treat2) * (1 + treatsInter[interaction(group1, group2)])
})

testFormula <- function(formula, data) {
    print(summary(fm <- lmer(formula, data, control=lmerControl(optimizer="bobyqa"))))
    print(summary(rm <- rlmer(formula, data, rho.e = cPsi, rho.b = cPsi, init = lmerNoFit)))
    ranef.fm <- ranef(fm, condVar=FALSE)
    stopifnot(all.equal(coef(fm), coef(rm), tolerance = 1e-1, check.attributes = FALSE),
              all.equal(fixef(fm), fixef(rm), tolerance = 1e-2, check.attributes = FALSE),
              all.equal(ranef.fm , ranef(rm), tolerance = 1e-1, check.attributes = FALSE))
    invisible(list(fm, rm))
}

testFormula(resp1 ~ (1|group1), data)
testFormula(resp2 ~ (1|group1) + (1|group2), data)
testFormula(resp3 ~ treat1 + (1 + treat1|group1) + (1|group2), data)
testFormula(resp4 ~ treat1 + treat2 + (1 + treat1|group1) + (1 + treat2|group2), data)
testFormula(resp5 ~ treat1 + treat2 + (1 + treat1 + treat2|group1) + (1 + treat2|group2), data)
testFormula(resp6 ~ treat1 + treat2 + (1 + treat1 + treat2|group1) + (1 + treat1 + treat2|group2), data)
testFormula(resp7 ~ treat1*treat2 + (1 + treat1*treat2|group1) + (1 + treat1*treat2|group2), data)
