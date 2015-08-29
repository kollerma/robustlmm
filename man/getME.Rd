\name{getME}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\alias{getME.rlmerMod}
\alias{getME}
\alias{theta}
\usage{
\S3method{getME}{rlmerMod}(object,
      name = c("X", "Z", "Zt", "Ztlist", "y", "mu",
               "u", "b.s", "b", "Gp", "Tp", "Lambda",
               "Lambdat","A", "U_b", "Lind", "sigma",
               "flist", "beta", "theta", "n_rtrms",
               "n_rfacs", "cnms", "devcomp", "offset",
               "lower", "rho_e", "rho_b", "rho_sigma_e",
               "rho_sigma_b", "M", "w_e", "w_b",
               "w_b_vector", "w_sigma_e", "w_sigma_b",
               "w_sigma_b_vector", "is_REML"), \dots)

  theta(object)
}
\arguments{
  \item{object}{a fitted mixed-effects model of class
    \code{"\linkS4class{rlmerMod}"}, i.e. typically the
    result of \code{\link{rlmer}()}.}

  \item{name}{a character string specifying the name of the
    \dQuote{component}.  Possible values are:\cr \describe{
  \item{X}{fixed-effects model matrix} \item{Z}{random-effects model
  matrix} \item{Zt}{transpose of random-effects model matrix}
  \item{Ztlist}{list of components of the transpose of the
  random-effects model matrix, separated by individual variance
  component} \item{y}{response vector} \item{mu}{conditional mean of
  the response} \item{u}{conditional mode of the \dQuote{spherical}
  random effects variable} \item{b.s}{synonym for \dQuote{u}}
  \item{b}{onditional mode of the random effects variable}
  \item{Gp}{groups pointer vector.  A pointer to the beginning of each
  group of random effects corresponding to the random-effects terms.}
  \item{Tp}{theta pointer vector.  A pointer to the beginning of the
  theta sub-vectors corresponding to the random-effects terms,
  beginning with 0 and including a final element giving the total
  number of random effects} \item{Lambda}{relative covariance factor
  of the random effects.}  \item{U_b}{synonym for \dQuote{Lambda}}
  \item{Lambdat}{transpose of the relative covariance factor of the
  random effects.}  \item{Lind}{index vector for inserting elements of
  \eqn{\theta}{theta} into the nonzeros of \eqn{\Lambda}{Lambda}}
  \item{A}{Scaled sparse model matrix (class
  \code{"\link[Matrix:dgCMatrix-class]{dgCMatrix}"}) for the unit,
  orthogonal random effects, \eqn{U}, equal to \code{getME(.,"Zt")
  \%*\% getME(.,"Lambdat")}} \item{sigma}{residual standard error}
  \item{flist}{a list of the grouping variables (factors) involved in
  the random effect terms} \item{beta}{fixed-effects parameter
  estimates (identical to the result of \code{\link{fixef}}, but
  without names)} \item{theta}{random-effects parameter estimates:
  these are parameterized as the relative Cholesky factors of each
  random effect term} \item{n_rtrms}{number of random-effects terms}
  \item{n_rfacs}{number of distinct random-effects grouping factors}
  \item{cnms}{ the "component names", a 'list'.} \item{devcomp}{a
  list consisting of a named numeric vector, \dQuote{cmp}, and a named
  integer vector, \dQuote{dims}, describing the fitted model}
  \item{offset}{model offset} \item{lower}{lower bounds on model
  parameters (random effects parameters only)} \item{rho_e}{rho
  function used for the residuals} \item{rho_b}{list of rho functions
  used for the random effects} \item{rho_sigma_e}{rho function used
  for the residuals when estimating sigma} \item{rho_sigma_b}{list of
  rho functions used for the random effects when estimating the
  covariance parameters} \item{M}{list of matrices, blocks of the
  Henderson's equations and the matrices used for computing the linear
  approximations of the estimates of beta and spherical random
  effects.}  \item{w_e}{robustness weights associated with the
  observations} \item{w_b}{robustness weights associated with the
  spherical random effects, returned in the same format as
  \code{\link{ranef}()}} \item{w_b_vector}{robustness weights
  associated with the spherical random effects, returned as one long
  vector} \item{w_sigma_e}{robustness weights associated with the
  observations when estimating sigma} \item{w_sigma_b}{robustness
  weights associated with the spherical random effects when estimating
  the covariance parameters, returned in the same format as
  \code{\link{ranef}()}} \item{w_sigma_b_vector}{robustness weights
  associated with the spherical random effects when estimating the
  covariance parameters, returned as one long vector}
  \item{is_REML}{returns \code{TRUE} for rlmerMod-objects (for
    compatibility with lme4)}
}% describe
}% name
  \item{...}{potentially further arguments passed to and from methods;
    none here at the moment.}
}
\value{ Unspecified, as very much
  depending on the \code{\link{name}}.
}
\description{ Extract (or
  \dQuote{get}) \dQuote{components} -- in a generalized sense -- from
  a fitted mixed-effects model, i.e. from an object of class
  \code{"\linkS4class{rlmerMod}"} or \code{"\linkS4class{merMod}"}.

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
		    check.attributes=FALSE, tolerance = 0))

## All that can be accessed [potentially ..]:
(nmME <- eval(formals(robustlmm:::getME.rlmerMod)$name))
\dontshow{
## internal consistency check ensuring that all work:
## "try(.)" because some are not yet implemented:
str(parts <- sapply(nmME, function(nm) try(getME(fm1, nm)),
                    simplify=FALSE))
}\% dont..
stopifnot(all.equal(theta(fm1), getME(fm1, "theta")))
}
\seealso{
  \code{\link{getCall}()}; more standard methods for
  rlmerMod objects, such as \code{\link{ranef}},
  \code{\link{fixef}}, \code{\link{vcov}}, etc.: see
  \code{methods(class="rlmerMod")}
}
\keyword{utilities}

