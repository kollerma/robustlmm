\name{getME}
\alias{getME.rlmerMod}
\alias{getME}
\alias{theta}
\title{Extract or Get Generalize Components from a Fitted Mixed Effects Model}
\alias{getME.rlmerMod}
\alias{getME}
\alias{theta}
\usage{
\S3method{getME}{rlmerMod}(
  object,
  name = c("X", "Z", "Zt", "Ztlist", "mmList", "y", "mu", "u", "b.s", "b", "Gp", "Tp",
    "Lambda", "Lambdat", "Tlist", "A", "U_b", "Lind", "sigma", "flist", "fixef", "beta",
    "theta", "ST", "is_REML", "n_rtrms", "n_rfacs", "N", "n", "p", "q", "p_i", "l_i",
    "q_i", "k", "m_i", "m", "cnms", "devcomp", "offset", "lower", "rho_e", "rho_b",
    "rho_sigma_e", "rho_sigma_b", "M", "w_e", "w_b", "w_b_vector", "w_sigma_e",
    "w_sigma_b", "w_sigma_b_vector"),
  ...
)

theta(object)
}
\arguments{
\item{object}{a fitted mixed-effects model of class
\code{\linkS4class{rlmerMod}}, i.e. typically the result of
\code{\link{rlmer}()}.}

\item{name}{a character string specifying the name of the
\dQuote{component}.  Possible values are:\cr
\describe{
    \item{\code{"X"}:}{fixed-effects model matrix}
    \item{\code{"Z"}:}{random-effects model matrix}
    \item{\code{"Zt"}:}{transpose of random-effects model matrix}
    \item{\code{"Ztlist"}:}{list of components of the transpose of the random-effects model matrix,
             separated by individual variance component}
    \item{\code{"mmList"}:}{list of raw model matrices associated with random effects terms}
    \item{\code{"y"}:}{response vector}
    \item{\code{"mu"}:}{conditional mean of the response}
    \item{\code{"u"}:}{conditional mode of the \dQuote{spherical} random effects variable}
    \item{\code{"b.s"}:}{synonym for \dQuote{u}}
    \item{\code{"b"}:}{conditional mode of the random effects variable}
    \item{\code{"Gp"}:}{groups pointer vector.  A pointer to the beginning of each group
              of random effects corresponding to the random-effects terms.}
    \item{\code{"Tp"}:}{theta pointer vector.  A pointer to the beginning
              of the theta sub-vectors corresponding to the
              random-effects terms, beginning with 0 and including
              a final element giving the total number of random effects}
    \item{\code{"Lambda"}:}{relative covariance factor of the random effects.}
    \item{\code{"U_b"}:}{synonym for \dQuote{Lambda}}
    \item{\code{"Lambdat"}:}{transpose of the relative covariance factor of the random effects.}
    \item{\code{"Lind"}:}{index vector for inserting elements of \eqn{\theta}{theta} into the
                nonzeros of \eqn{\Lambda}{Lambda}}
    \item{\code{"A"}:}{Scaled sparse model matrix (class
     \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}) for
     the unit, orthogonal random effects, \eqn{U},
      equal to \code{getME(.,"Zt") \%*\% getME(.,"Lambdat")}}
    \item{\code{"sigma"}:}{residual standard error}
    \item{\code{"flist"}:}{a list of the grouping variables (factors) involved in the random effect terms}
    \item{\code{"fixef"}:}{fixed-effects parameter estimates}
    \item{\code{"beta"}:}{fixed-effects parameter estimates (identical to the result of \code{\link{fixef}}, but without names)}
    \item{\code{"theta"}:}{random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term}
    \item{\code{"ST"}:}{A list of S and T factors in the TSST' Cholesky
              factorization of the relative variance matrices of the random
              effects associated with each random-effects term.  The unit lower
              triangular matrix, \eqn{T}, and the diagonal matrix, \eqn{S}, for
              each term are stored as a single matrix with diagonal elements
              from \eqn{S} and off-diagonal elements from \eqn{T}.}
    \item{\code{"is_REML"}:}{returns \code{TRUE} for rlmerMod-objects (for compatibility with lme4)}
    \item{\code{"n_rtrms"}:}{number of random-effects terms}
    \item{\code{"n_rfacs"}:}{number of distinct random-effects grouping factors}
    \item{\code{"N"}:}{number of rows of \code{X}}
    \item{\code{"n"}:}{length of the response vector, \code{y}}
    \item{\code{"p"}:}{number of columns of the fixed effects model matrix, \code{X}}
    \item{\code{"q"}:}{number of columns of the random effects model matrix, \code{Z}}
    \item{\code{"p_i"}:}{numbers of columns of the raw model matrices, \code{mmList}}
    \item{\code{"l_i"}:}{numbers of levels of the grouping factors}
    \item{\code{"q_i"}:}{numbers of columns of the term-wise model matrices, \code{ZtList}}
    \item{\code{"k"}:}{number of random effects terms}
    \item{\code{"m_i"}:}{numbers of covariance parameters in each term}
    \item{\code{"m"}:}{total number of covariance parameters, i.e., the
                       same as \code{dim@nth} below.}
    \item{\code{"cnms"}:}{the \dQuote{component names}, a \sQuote{list}.}
    \item{\code{"devcomp"}:}{a list consisting of a named numeric vector,
        \code{cmp}, and a named integer vector, \code{dims}, describing
        the fitted model.  The elements of \code{cmp} are:\cr
        \describe{
            \item{ldL2}{always NA, for consistency with lme4 output}
            \item{ldRX2}{always NA, for consistency with lme4 output}
            \item{wrss}{always NA, for consistency with lme4 output}
            \item{ussq}{always NA, for consistency with lme4 output}
            \item{pwrss}{always NA, for consistency with lme4 output}
            \item{drsum}{always NA, for consistency with lme4 output}
            \item{REML}{always NA, for consistency with lme4 output}
            \item{dev}{always NA, for consistency with lme4 output}
            \item{sigmaML}{ML estimate of residual standard deviation}
            \item{sigmaREML}{REML estimate of residual standard deviation}
        } The elements of \code{dims} are:\cr
        \describe{
            \item{N}{number of rows of \code{X}}
            \item{n}{length of \code{y}}
            \item{p}{number of columns of \code{X}}
            \item{nmp}{\code{n-p}}
            \item{nth}{length of \code{theta}}
            \item{q}{number of columns of \code{Z}}
            \item{nAGQ}{see \code{\link{glmer}}}
            \item{compDev}{see \code{\link{glmerControl}}}
            \item{useSc}{\code{TRUE} if model has a scale parameter}
            \item{reTrms}{number of random effects terms}
            \item{REML}{\code{0} indicates the model was fitted by maximum
                likelihood, any other positive integer indicates fitting by
                restricted maximum likelihood}
            \item{GLMM}{\code{TRUE} if a GLMM}
            \item{NLMM}{\code{TRUE} if an NLMM}
        }
    }
    \item{\code{"offset"}:}{model offset}
    \item{\code{"lower"}:}{lower bounds on random-effects model
        parameters (i.e, "theta" parameters). In order to constrain
        random effects covariance matrices to be semi-positive-definite,
        this vector is equal to 0 for elements of
        the \code{theta} vector corresponding to diagonal elements of
        the Cholesky factor, \code{-Inf}
        otherwise. (\code{getME(.,"lower")==0} can be used as a test to
        identify diagonal elements, as in \code{isSingular}.)
    }
    \item{\code{"rho_e"}:}{rho function used for the residuals}
    \item{\code{"rho_b"}:}{list of rho functions used for the random effects}
    \item{\code{"rho_sigma_e"}:}{rho function used for the residuals when estimating sigma}
    \item{\code{"rho_sigma_b"}:}{list of rho functions used for the random effects when estimating the covariance parameters}
    \item{\code{"M"}:}{list of matrices, blocks of the Henderson's equations and the matrices used for computing the linear approximations of the estimates of beta and spherical random effects.}
    \item{\code{"w_e"}:}{robustness weights associated with the observations}
    \item{\code{"w_b"}:}{robustness weights associated with the spherical random effects, returned in the same format as \code{\link{ranef}()}}
    \item{\code{"w_b_vector"}:}{robustness weights associated with the spherical random effects, returned as one long vector}
    \item{\code{"w_sigma_e"}:}{robustness weights associated with the observations when estimating sigma}
    \item{\code{"w_sigma_b"}:}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned in the same format as \code{\link{ranef}()}}
    \item{\code{"w_sigma_b_vector"}:}{robustness weights associated with the spherical random effects when estimating the covariance parameters, returned as one long vector}

     %% -- keep at the very end:
     \item{\code{"ALL"}:}{get all of the above as a \code{\link{list}}.}
}}

\item{...}{potentially further arguments; not here.}
}
\value{
Unspecified, as very much depending on the \code{\link{name}}.
}
\description{
Extract (or \dQuote{get}) \dQuote{components} -- in a generalized
sense -- from a fitted mixed-effects model, i.e. from an object
of class \code{\linkS4class{rlmerMod}} or \code{\linkS4class{merMod}}.
}
\details{
The function \code{theta} is short for \code{getME(, "theta")}.

The goal is to provide \dQuote{everything a user may want} from a fitted
\code{rlmerMod} object \emph{as far} as it is not available by methods, such
as \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{vcov}}, etc.
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
	  all.equal(fixef(fm1), b1 <- getME(fm1, "beta"),
		    check.attributes=FALSE, tolerance = 0))

## A way to get *all* getME()s :
## internal consistency check ensuring that all work:
parts <- getME(fm1, "ALL")
str(parts, max=2)
stopifnot(identical(Z,  parts $ Z),
          identical(b1, parts $ beta))
stopifnot(all.equal(theta(fm1), getME(fm1, "theta")))
}
\seealso{
\code{\link{getCall}()};
more standard methods for rlmerMod objects, such as \code{\link{ranef}},
\code{\link{fixef}}, \code{\link{vcov}}, etc.:
see \code{methods(class="rlmerMod")}
}
\keyword{utilities}
