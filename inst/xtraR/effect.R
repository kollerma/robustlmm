## Make effect package work for rlmerMod objects
require(effects)
require(robustlmm)

rlmerMod.to.glm <- function(mod) {
  family <- gaussian()
  link <- family$link
  family <- family$family
  cl <- mod@call
  cl$control <- glm.control(epsilon = 1)
  .m <- match(c("formula", "family", "data", "weights", "subset",
                "na.action", "offset", "model", "contrasts"), names(cl),
              0L)
  cl <- cl[c(1L, .m)]
  cl[[1L]] <- as.name("glm")
  cl$formula <- effects:::fixmod(as.formula(cl$formula))
  mod2 <- eval(cl)
  mod2$coefficients <- lme4::fixef(mod)
  mod2$vcov <- as.matrix(vcov(mod))
  mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
  mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
  mod2$weights <- as.vector(with(mod2, prior.weights * (family$mu.eta(linear.predictors)^2/family$variance(fitted.values))))
  mod2$residuals <- with(mod2, prior.weights * (y - fitted.values)/weights)
  class(mod2) <- c("fakeglm", class(mod2))
  mod2
}

Effect.rlmerMod <- function(focal.predictors, mod, ...) {
  result <- Effect(focal.predictors, rlmerMod.to.glm(mod),
                   ...)
  result$formula <- as.formula(formula(mod))
  result
}

fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake, REML = FALSE)
plot(Effect(c("recipe", "temperature"), fm1))

## REML = FALSE is not possible.
rfm1 <- rlmer(angle ~ recipe * temperature + (1|recipe:replicate), cake)
plot(Effect(c("recipe", "temperature"), rfm1))


