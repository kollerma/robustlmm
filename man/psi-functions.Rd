\name{psi-functions}
\alias{cPsi}
\alias{huberPsi}
\alias{smoothPsi}
\title{Classical, Huber and smoother Huber psi- or rho-functions}
\description{
  \eqn{\psi}{Psi}-functions are used by \code{\link{rlmer}}
  in the estimating equations and to compute robustness
  weights. Change tuning parameters using
  \code{\link{chgDefaults}} and convert to squared
  robustness weights using the \code{\link{psi2propII}}
  function.
}
\details{
  The \bold{\dQuote{classical} \eqn{\psi}{psi}-function
  \code{cPsi}} can be used to get a non-robust, i.e.,
  classical, fit. The \code{psi} slot equals the identity
  function, and the \code{rho} slot equals quadratic
  function. Accordingly, the robustness weights will always
  be 1 when using \code{cPsi}.

  The \bold{Huber \eqn{\psi}{psi}-function \code{huberPsi}}
  is identical to the one in the package \code{robustbase}.
  The \code{psi} slot equals the identity function within
  \eqn{\pm k}{+-k} (where \eqn{k}{k} is the tuning
  parameter). Outside this interval it is equal to \eqn{\pm
  k}{+-k}. The \code{rho} slot equals the quadratic
  function within \eqn{\pm k}{+-k} and a linear function
  outside.

  The \bold{smoothed Huber \eqn{\psi}{psi}-function} is
  very similar to the regular Huber
  \eqn{\psi}{psi}-function. Instead of a sharp bend like
  the Huber function, the smoothe Huber function bends
  smoothly. The first tuning contant, k, can be compared to
  the tuning constant of the original Huber function. The
  second tuning constant, s, determines the smoothness of
  the bend.
}
\examples{
plot(cPsi)
plot(huberPsi)
plot(smoothPsi)
curve(cPsi@psi(x), -3, 3)
curve(smoothPsi@psi(x, 1.345, 10), -3, 3, add=TRUE, col="red")
curve(huberPsi@psi(x, 1.345), -3, 3, add=TRUE, col="blue")
}
\seealso{
  \code{\link{chgDefaults}} and \code{\link{psi2propII}}
  for changing tuning parameters;
  \code{\link{psi_func_cached-class}} and
  \code{\link{psi_func-class}} for a more detailed
  description of the slots; \code{\link{psiFuncCached}} for
  a constructor function to create custom
  \eqn{\psi}{psi}-functions.
}
