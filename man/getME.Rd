\name{getME}
\alias{getME}
\alias{theta}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\usage{
  getME(object, name = c("X", "Z", "Zt", "u", "b.s", "b",
    "Gp", "Lambda", "Lambdat", "U_b", "Lind", "flist",
    "beta", "theta", "n_rtrms", "devcomp", "offset",
    "lower", "rho_e", "rho_b", "rho_sigma_e",
    "rho_sigma_b", "M", "w_e", "w_b", "w_b_vector",
    "w_sigma_e", "w_sigma_b", "w_sigma_b_vector"))

  theta(object)
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
  \code{"\linkS4class{rlmerMod}"}, i.e. typically the
  result of \code{\link{rlmer}()}.}

  \item{name}{a character string specifying the name of the
  \dQuote{component}.  Possible values are:\cr \describe{
  \item{X}{fixed-effects model matrix}
  \item{Z}{random-effects model matrix} \item{Zt}{transpose
  of random-effects model matrix} \item{u}{conditional mode
  of the \dQuote{spherical} random effects variable}
  \item{b.s}{synonym for \dQuote{u}} \item{b}{onditional
  mode of the random effects variable} \item{Gp}{groups
  pointer vector.  A pointer to the beginning of each group
  of random effects corresponding to the random-effects
  terms.} \item{Lambda}{relative covariance factor of the
  random effects.} \item{U_b}{synonym for \dQuote{Lambda}}
  \item{Lambdat}{transpose of the relative covariance
  factor of the random effects.} \item{Lind}{index vector
  for inserting elements of \eqn{\theta}{theta} into the
  nonzeros of \eqn{\Lambda}{Lambda}} \item{flist}{a list of
  the grouping variables (factors) involved in the random
  effect terms} \item{beta}{fixed-effects parameter
  estimates (identical to the result of
  \code{\link{fixef}}, but without names)}
  \item{theta}{random-effects parameter estimates: these
  are parameterized as the relative Cholesky factors of
  each random effect term} \item{n_rtrms}{number of
  random-effects terms} \item{devcomp}{a list consisting of
  a named numeric vector, \dQuote{cmp}, and a named integer
  vector, \dQuote{dims}, describing the fitted model}
  \item{offset}{model offset} \item{lower}{lower bounds on
  model parameters (random effects parameters only)}
  \item{rho_e}{rho function used for the residuals}
  \item{rho_b}{list of rho functions used for the random
  effects} \item{rho_sigma_e}{rho function used for the
  residuals when estimating sigma} \item{rho_sigma_b}{list
  of rho functions used for the random effects when
  estimating the covariance parameters} \item{M}{list of
  matrices, blocks of the Henderson's equations and the
  matrices used for computing the linear approximations of
  the estimates of beta and spherical random effects.}
  \item{w_e}{robustness weights associated with the
  observations} \item{w_b}{robustness weights associated
  with the spherical random effects, returned in the same
  format as \code{\link{ranef}()}}
  \item{w_b_vector}{robustness weights associated with the
  spherical random effects, returned as one long vector}
  \item{w_sigma_e}{robustness weights associated with the
  observations when estimating sigma}
  \item{w_sigma_b}{robustness weights associated with the
  spherical random effects when estimating the covariance
  parameters, returned in the same format as
  \code{\link{ranef}()}} \item{w_sigma_b_vector}{robustness
  weights associated with the spherical random effects when
  estimating the covariance parameters, returned as one
  long vector} }}
}
\value{
  Unspecified, as very much depending on the
  \code{\link{name}}.
}
\description{
  Extract (or \dQuote{get}) \dQuote{components} -- in a
  generalized sense -- from a fitted mixed-effects model,
  i.e. from an object of class
  \code{"\linkS4class{rlmerMod}"} or
  \code{"\linkS4class{merMod}"}.

  The function \code{theta} is short for \code{getME(,
  "theta")}.
}
\details{
  The goal is to provide \dQuote{everything a user may
  want} from a fitted \code{"rlmerMod"} object \emph{as
  far} as it is not available by methods, such as
  \code{\link{fixef}}, \code{\link{ranef}},
  \code{\link{vcov}}, etc.
}
\examples{
## shows many methods you should consider *before* using getME():
methods(class = "rlmerMod")

## doFit = FALSE to speed up example
(fm1 <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
              method="DASvar", doFit=FALSE))
Z <- getME(fm1, "Z")
stopifnot(is(Z, "CsparseMatrix"),
          c(180,36) == dim(Z),
	  all.equal(fixef(fm1), getME(fm1, "beta"),
		    check.attr=FALSE, tol = 0))

## All that can be accessed [potentially ..]:
(nmME <- eval(formals(getME)$name))
\dontshow{
## internal consistency check ensuring that all work:
## "try(.)" because some are not yet implemented:
str(parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
                    simplify=FALSE))
}\% dont..
stopifnot(all.equal(theta(fm1), getME(fm1, "theta")))
}
\seealso{
  \code{\link{getCall}()}, More standard methods for
  rlmerMod objects, such as \code{\link{ranef}},
  \code{\link{fixef}}, \code{\link{vcov}}, etc.: see
  \code{methods(class="rlmerMod")}
}
\keyword{utilities}

