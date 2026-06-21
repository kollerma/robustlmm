########################################################
## Contamination strategies for the MC breakdown study  #
########################################################
##
## Each strategy is a function (data, eps, ...) -> contaminated_data.
## eps is interpreted as a fraction of observations (or, for RE-level
## strategies on balanced designs, equivalently the fraction of groups).
##
## Observation-level strategies replace y[i] with abs(y[i]) * magnitude
## (matching the existing breakdown.R: magnitude = 1e6). RE-level
## strategies add a shift of magnitude_sigma * sigma_e to all observations
## of the selected groups.

##' Sort row indices so that group 1 is filled first, then group 2, etc.
.indexWithin <- function(data, groupname) {
    order(data[[groupname]], seq_len(nrow(data)))
}

##' Sort row indices so that one observation from each group is taken in
##' turn (round 1: first obs of each group; round 2: second obs of each
##' group; ...). Useful for spreading contamination across groups.
.indexRoundRobin <- function(data, groupname) {
    n <- nrow(data)
    g <- as.factor(data[[groupname]])
    by_group <- split(seq_len(n), g)
    max_per_group <- max(lengths(by_group))
    o <- unlist(lapply(seq_len(max_per_group), function(j) {
        sapply(by_group, function(idx)
            if (j <= length(idx)) idx[j] else NA_integer_)
    }))
    o[!is.na(o)]
}

##' Apply an additive y-shift of magnitude_sigma * sigma_e to a subset
##' of rows.
.applyShift <- function(data, contam, yname, magnitude_sigma, sigma_e) {
    if (length(contam) == 0L) return(data)
    if (is.null(sigma_e)) sigma_e <- stats::sd(data[[yname]])
    data[[yname]][contam] <- data[[yname]][contam] + magnitude_sigma * sigma_e
    data
}

##' Saturate one group at a time (within-group first).
strategy_within <- function(data, eps, yname, groupname,
                             magnitude_sigma = 10, sigma_e = NULL) {
    n <- nrow(data)
    n_contam <- ceiling(eps * n)
    contam <- if (n_contam > 0)
        .indexWithin(data, groupname)[seq_len(n_contam)]
    else integer(0)
    .applyShift(data, contam, yname, magnitude_sigma, sigma_e)
}

##' Round-robin across groups: one obs per group, then next obs per group, ...
strategy_round_robin <- function(data, eps, yname, groupname,
                                  magnitude_sigma = 10, sigma_e = NULL) {
    n <- nrow(data)
    n_contam <- ceiling(eps * n)
    contam <- if (n_contam > 0)
        .indexRoundRobin(data, groupname)[seq_len(n_contam)]
    else integer(0)
    .applyShift(data, contam, yname, magnitude_sigma, sigma_e)
}

##' Random observations.
strategy_random_obs <- function(data, eps, yname,
                                 magnitude_sigma = 10, sigma_e = NULL,
                                 seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    n <- nrow(data)
    n_contam <- ceiling(eps * n)
    contam <- if (n_contam > 0) sample.int(n, n_contam) else integer(0)
    .applyShift(data, contam, yname, magnitude_sigma, sigma_e)
}

##' High-leverage observations (top of a continuous covariate).
strategy_high_leverage <- function(data, eps, yname, leverage_var,
                                    magnitude_sigma = 10, sigma_e = NULL) {
    n <- nrow(data)
    n_contam <- ceiling(eps * n)
    contam <- if (n_contam > 0)
        order(data[[leverage_var]], decreasing = TRUE)[seq_len(n_contam)]
    else integer(0)
    .applyShift(data, contam, yname, magnitude_sigma, sigma_e)
}

##' RE-level contamination: shift all observations of a fraction eps of
##' the groups by magnitude_sigma * sigma_e.
strategy_re_shift <- function(data, eps, yname, groupname,
                               magnitude_sigma = 4, sigma_e = NULL,
                               seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(sigma_e)) sigma_e <- stats::sd(data[[yname]])
    g <- as.factor(data[[groupname]])
    levs <- levels(g)
    n_groups <- length(levs)
    n_contam_groups <- ceiling(eps * n_groups)
    if (n_contam_groups == 0) return(data)
    contam_levs <- sample(levs, n_contam_groups)
    rows <- as.character(g) %in% contam_levs
    shift <- magnitude_sigma * sigma_e
    data[[yname]][rows] <- data[[yname]][rows] + shift
    data
}

##' RE-level slope contamination (sleepstudy-specific):
##' for a fraction eps of the subjects, shift each obs by
##' magnitude_sigma * sigma_e * (covariate_value - mean(covariate)) so the
##' effect is concentrated at high covariate values.
strategy_re_shift_slope <- function(data, eps, yname, groupname,
                                     covariate, magnitude_sigma = 4,
                                     sigma_e = NULL, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(sigma_e)) sigma_e <- stats::sd(data[[yname]])
    g <- as.factor(data[[groupname]])
    levs <- levels(g)
    n_contam_groups <- ceiling(eps * length(levs))
    if (n_contam_groups == 0) return(data)
    contam_levs <- sample(levs, n_contam_groups)
    rows <- as.character(g) %in% contam_levs
    x <- data[[covariate]]
    x_centered <- x - mean(x)
    data[[yname]][rows] <- data[[yname]][rows] +
        magnitude_sigma * sigma_e * x_centered[rows]
    data
}

##' Per-design strategy registry.
##'
##' Returns a list mapping strategy name to a function (data, eps, seed)
##' -> contaminated_data. The contamination magnitude is fixed at
##' construction time: \code{magnitude_sigma} (default 10) is the y-shift
##' in units of \code{sigma_e} for observation-level strategies, and the
##' RE-shift in units of \code{sigma_e} for RE-level strategies.
##' \code{sigma_e} should be the true residual scale of the design
##' (e.g. \code{sigma(lmer(...))} on the uncontaminated data).
strategiesForDesign <- function(design_name, sigma_e, magnitude_sigma = 10) {
    m  <- magnitude_sigma
    se <- sigma_e
    switch(design_name,
        dyestuff = list(
            within       = function(data, eps, seed)
                strategy_within(data, eps, "Yield", "Batch",
                                 magnitude_sigma = m, sigma_e = se),
            round_robin  = function(data, eps, seed)
                strategy_round_robin(data, eps, "Yield", "Batch",
                                      magnitude_sigma = m, sigma_e = se),
            random_obs   = function(data, eps, seed)
                strategy_random_obs(data, eps, "Yield",
                                     magnitude_sigma = m, sigma_e = se,
                                     seed = seed),
            re_shift     = function(data, eps, seed)
                strategy_re_shift(data, eps, "Yield", "Batch",
                                   magnitude_sigma = m, sigma_e = se,
                                   seed = seed)
        ),
        penicillin = list(
            within_plate       = function(data, eps, seed)
                strategy_within(data, eps, "diameter", "plate",
                                 magnitude_sigma = m, sigma_e = se),
            within_sample      = function(data, eps, seed)
                strategy_within(data, eps, "diameter", "sample",
                                 magnitude_sigma = m, sigma_e = se),
            round_robin_plate  = function(data, eps, seed)
                strategy_round_robin(data, eps, "diameter", "plate",
                                      magnitude_sigma = m, sigma_e = se),
            round_robin_sample = function(data, eps, seed)
                strategy_round_robin(data, eps, "diameter", "sample",
                                      magnitude_sigma = m, sigma_e = se),
            random_obs         = function(data, eps, seed)
                strategy_random_obs(data, eps, "diameter",
                                     magnitude_sigma = m, sigma_e = se,
                                     seed = seed),
            re_shift_plate     = function(data, eps, seed)
                strategy_re_shift(data, eps, "diameter", "plate",
                                   magnitude_sigma = m, sigma_e = se,
                                   seed = seed),
            re_shift_sample    = function(data, eps, seed)
                strategy_re_shift(data, eps, "diameter", "sample",
                                   magnitude_sigma = m, sigma_e = se,
                                   seed = seed)
        ),
        sleepstudy = list(
            within_subject    = function(data, eps, seed)
                strategy_within(data, eps, "Reaction", "Subject",
                                 magnitude_sigma = m, sigma_e = se),
            round_robin       = function(data, eps, seed)
                strategy_round_robin(data, eps, "Reaction", "Subject",
                                      magnitude_sigma = m, sigma_e = se),
            high_leverage     = function(data, eps, seed)
                strategy_high_leverage(data, eps, "Reaction", "Days",
                                        magnitude_sigma = m, sigma_e = se),
            random_obs        = function(data, eps, seed)
                strategy_random_obs(data, eps, "Reaction",
                                     magnitude_sigma = m, sigma_e = se,
                                     seed = seed),
            re_shift          = function(data, eps, seed)
                strategy_re_shift(data, eps, "Reaction", "Subject",
                                   magnitude_sigma = m, sigma_e = se,
                                   seed = seed),
            re_shift_slope    = function(data, eps, seed)
                strategy_re_shift_slope(data, eps, "Reaction", "Subject",
                                         "Days",
                                         magnitude_sigma = m, sigma_e = se,
                                         seed = seed)
        ),
        stop("Unknown design: ", design_name))
}
