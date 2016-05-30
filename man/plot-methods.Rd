\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,Rcpp_SmoothPsi-method}
\alias{plot,Rcpp_HuberPsi-method}
\alias{plot,Rcpp_PsiFunction-method}
\alias{plot,Rcpp_PsiFunction to Prop II PsiFunction wrapper-method}
\title{Plot an Object of the "Psi Function" Class}
\description{
  The \code{\link{plot}} method objects of class
  \code{PsiFunction} simply visualizes the
  \eqn{\rho()}, \eqn{\psi()}, and weight functions and their
  derivatives.
}
\usage{
\S4method{plot}{Rcpp_SmoothPsi}(x, y,
     which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
     main = "full", % shortMain = FALSE,
     col = c("black", "red3", "blue3", "dark green", "light green"),
     leg.loc = "right", ...)
\S4method{plot}{Rcpp_HuberPsi}(x, y,
     which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
     main = "full", % shortMain = FALSE,
     col = c("black", "red3", "blue3", "dark green", "light green"),
     leg.loc = "right", ...)
\S4method{plot}{Rcpp_PsiFunction}(x, y,
     which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
     main = "full", % shortMain = FALSE,
     col = c("black", "red3", "blue3", "dark green", "light green"),
     leg.loc = "right", ...)
\S4method{plot}{`Rcpp_PsiFunction to Prop II PsiFunction wrapper`}(x, y,
     which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
     main = "full", % shortMain = FALSE,
     col = c("black", "red3", "blue3", "dark green", "light green"),
     leg.loc = "right", ...)
}
\arguments{
  \item{x}{instance of class \code{PsiFunction} to be plotted}
  \item{y}{(optional) vector of abscissa values (to plot object at).}
  \item{which}{\code{\link{character}} vector of slots to be included in
    plot; by default, all of the slots are included}
  \item{main}{string or logical indicating the kind of plot title;
  either \code{"full"}, \code{"short"} or \code{FALSE} which chooses a
  full, a short or no main title at all.}
  \item{col}{colors to be used for the different slots}
  \item{leg.loc}{legend placement, see also \code{x} argument
    of \code{\link{legend}}}
  \item{...}{passed to \code{\link{matplot}}}
}
\note{
  If you want to specify your own title, use \code{main=FALSE}, and a
  subsequent \code{\link{title}(...)} call.
}
\seealso{
  \code{\link{psi-functions}}.
}
\examples{
plot(huberPsiRcpp)
plot(huberPsiRcpp, which=c("psi", "Dpsi", "wgt"),
     main="short", leg = "topleft")

plot(smoothPsi)
## Plotting aspect ratio = 1:1 :
plot(smoothPsi, asp=1, main="short",
     which = c("psi", "Dpsi", "wgt", "Dwgt"))
}
\keyword{methods}

