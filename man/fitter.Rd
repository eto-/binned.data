\name{fitter}
\alias{fitter}
\title{Advanced fitter for binned data}
\description{
}
\usage{
fitter <- function (b, pdf, start, fixed=list(),
                    statistics=c("Neyman", "Pearson", "Poissonian", "Multinomial"),
                    method=c("NLS", "MIGRAD", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), maxiter=-1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{The binned data set to fit.}
  \item{pdf}{The parametric function to fit over data.}
  \item{start}{A list with the named parameter of PDF.}
  \item{fixed}{An optional list of named parameter of PDF fixed during the minimization.}
  \item{statistics}{The statistical definition of the chi square under minimization, see \sQuote{Details}.}
  \item{method}{The minimization algorithm, see \sQuote{Details}.} 
  \item{maxiter}{Max number of iteration (available only in nls).}
}
\details{
  \code{fitter} minimizes the following function:
  \deqn{\sum_i{chi2(pdf(bins_i, c(start, fixed)), counts_i)},}{Sum chi2(pdf(bins[i], c(start, fixed)), counts[i]),}
  where the chi2 functions are defined for bin counting according to \code{"Neyman"}, \code{"Pearson"} or as
  likelihood ratio (\code{"Poissonian"}, \code{"Multinomial"}). See \code{\link{chi2.pearson}} documentation.

  \code{fitter} provides several minimization algorithms: \code{"NLS"} internally uses 
  a Levenberg-Marquardt nonlinear least-squares, as provided by the \code{"minpack.lm"} package.
  \code{"MIGRAD"} internally uses a variable metric minimizer, as provided by the
  \code{"migrad"} package. The other algorithms are provided by \code{optim}.
}
\value{
  \code{fitter} return a list with the following objects:
  \item{r.c}{the coefficient matrix returned by \code{coef}.}
  \item{chi2}{the chi square value returned by the minimization.}
  \item{df}{the number of degree of freedom.}
  \item{f}{the object returned by the minimizer.}
}

\seealso{chi2.neyman}
\examples{
  b <- binned(rnorm(20, 100,10))
  f <- fitter(d, function (mids, p) p$c, start=list(c=5), statistics="N")
  f$r.c
}
