\name{chi2.pearson}
\alias{chi2.pearson}
\alias{chi2.neyman}
\alias{chi2.saturared.poisson}
\alias{chi2.saturared.multinomial}
\title{chi-square models for binned data analysis}
\description{
  TBW
}
\usage{
chi2.pearson(model, count)
chi2.neyman(model, count)
chi2.saturared.poisson(model, count)
chi2.saturared.multinomial(model, count)
}
\arguments{
  \item{model}{the numeric vector with the theoretical counts}
  \item{count}{the numeric vector with the experimental counts}
}
\details{
  \code{chi2.pearson} and \code{chi2.neyman} return the Pearson/Neyman chi-square for
  given experimental data and model.

  \code{chi2.saturared.poisson} and \code{chi2.saturared.multinomial} 
  return the asymptotic saturated Pearson chi-square for the given 
  experimental data and model as described by Baker and Cousins (1983).
}
\value{ The numeric value of chi-square is returned }
\references{
  S. Baker and R. D. Cousins (1984) Clarification of the use of 
  chi-square and likelihood functions in fits to histograms. 
  \emph{Nuclear Instrument and Methods}, \bold{221}, 437-442.
}
