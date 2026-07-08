##' Global options consulted by \pkg{robustlmm}
##'
##' \pkg{robustlmm} reads a small number of global options, set with
##' \code{\link[base]{options}()} and queried with
##' \code{\link[base]{getOption}()}. With one experimental exception
##' (the Monte-Carlo DAS-tau calibration below), none of them change
##' the fitted estimates; they only tune optional diagnostics and the
##' default behaviour of \code{summary} (see
##' \code{\link{rlmerMod-class}}).
##'
##' @section Degrees-of-freedom options:
##'   These control whether \code{summary(object)} computes the robust
##'   Satterthwaite degrees of freedom and \code{Pr(>|t|)} column by
##'   default (\code{df = "auto"}); see the \dQuote{Coefficient-table
##'   degrees of freedom} section of \code{\link{rlmerMod-class}}.
##' \describe{
##'   \item{\code{robustlmm.summary.df.max}}{Numeric, default \code{5000}.
##'     The size cutoff for computing the Satterthwaite df under
##'     \code{df = "auto"}. The cost of the df is a deterministic
##'     dimension-only \dQuote{workload} \eqn{W} -- one \eqn{O(n)} score
##'     evaluation per parameter-Jacobian column: \eqn{W = (p + q + 1 + L)
##'     n} for method \code{"DASvar"} and \eqn{W = 40\,L\,n} for
##'     \code{"DAStau"}, where \eqn{n} is the number of observations,
##'     \eqn{p} the number of fixed effects, \eqn{q} the number of random
##'     effects and \eqn{L} the number of variance parameters. If \eqn{W}
##'     exceeds this cutoff and no influence function is cached on the fit,
##'     \code{summary} falls back to the plain \code{Estimate} /
##'     \code{Std. Error} / \code{t value} table and prints a note. The
##'     rule is dimensionless, so the same fit behaves identically on every
##'     machine. The default \code{5000} computes the df by default up to
##'     about \eqn{n = 170} for a single random intercept fit with
##'     \code{"DASvar"} (\eqn{n = 125} for \code{"DAStau"}). Set it higher
##'     to show the df on larger fits, or to \code{0} to always skip the
##'     automatic computation (you can still request it with
##'     \code{summary(object, df = "satterthwaite")}).}
##' }
##'
##' @section Monte-Carlo DAS-tau calibration (EXPERIMENTAL):
##'   These options enable and tune an \emph{experimental} Monte-Carlo
##'   calibration of the DAS-tau fixed point in \code{\link{rlmer}}.
##'   Without it, \code{method = "DAStau"} computes the consistency
##'   factors for non-diagonal random-effect blocks by Gauss-Hermite
##'   quadrature, which is limited to blocks of dimension \eqn{\le 2};
##'   fits containing a larger block fall back to \code{method =
##'   "DASvar"} with a warning. The Monte-Carlo path lifts this
##'   restriction: it computes the same self-consistent fixed point by
##'   plain Monte-Carlo integration, which works for any block
##'   dimension, including structured (\code{cs}/\code{ar1}) and
##'   unstructured blocks of dimension \eqn{> 2}. It is
##'   simulation-validated -- in the companion ar1 simulation study it
##'   removes the small calibration residual that the \code{DASvar}
##'   approximation leaves in the fitted correlation -- but it is not
##'   backed by a finite-sample theorem, and it has been validated on
##'   clean Gaussian data only. For blocks of dimension \eqn{\le 2} the
##'   classical quadrature path remains the default and the Monte-Carlo
##'   path offers no improvement there.
##' \describe{
##'   \item{\code{robustlmm.dastau.mc}}{Logical, default \code{FALSE}.
##'     Master switch. When \code{TRUE}, \code{method = "DAStau"} uses
##'     the Monte-Carlo calibration for all non-diagonal random-effect
##'     blocks of dimension \eqn{> 2} (instead of falling back to
##'     \code{"DASvar"} for the whole fit). The Monte-Carlo sample is
##'     drawn once per fit (common random numbers), deterministically
##'     seeded and moment-matched to exact zero mean and identity
##'     second moment, so repeated fits are identical and the caller's
##'     \code{.Random.seed} is left untouched. \code{rlmer} emits a
##'     message when the experimental path is active.}
##'   \item{\code{robustlmm.dastau.mc.all}}{Logical, default
##'     \code{FALSE}. Research switch: also use the Monte-Carlo path
##'     for blocks of dimension \eqn{2}, replacing the Gauss-Hermite
##'     quadrature. Intended only for comparing the two calibration
##'     paths; it offers no improvement over the quadrature.}
##'   \item{\code{robustlmm.dasmc.nsim}}{Integer, default \code{1e5}.
##'     Number of Monte-Carlo draws. Larger values reduce the
##'     (deterministic, seed-dependent) residual calibration error at
##'     linear cost in time and memory.}
##'   \item{\code{robustlmm.dasmc.seed}}{Integer, default
##'     \code{20260703}. Seed for the common-random-numbers draw. Fits
##'     are deterministic given this option; change it (e.g. per
##'     replicate in a simulation) to decorrelate the residual
##'     Monte-Carlo calibration error across fits. The global RNG
##'     state is saved and restored around the draw.}
##' }
##'
##' @section Developer options:
##' \describe{
##'   \item{\code{robustlmm.check_rhs_optimisation}}{Logical, default
##'     \code{FALSE}. When \code{TRUE}, \code{\link{rlmer}} cross-checks
##'     the vectorised right-hand-side computation in the block-diagonal
##'     \eqn{\theta} update against an explicit per-block loop and stops
##'     on any discrepancy. Intended for development and debugging only;
##'     it adds redundant work and is not needed in normal use.}
##' }
##'
##' @examples
##' ## show the df on larger fits (raise the size cutoff)
##' \dontrun{
##' options(robustlmm.summary.df.max = 20000)
##' }
##'
##' @name robustlmm-options
##' @seealso \code{\link{rlmerMod-class}}, \code{\link{rlmer}}
##' @keywords misc
NULL
