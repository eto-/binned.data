\name{binned}
\alias{binned}
\alias{as.binned}
\alias{as.binned.list}
\alias{as.binned.histogram}
\alias{as.binned.data.frame}
\title{Binned data}
\description{Convert existing data to a binned S3 class}
\usage{
binned(counts, bins, breaks=NULL)

\S3method{as.binned}{histogram}(h)

\S3method{as.binned}{data.frame}(d, breaks=NULL)

\S3method{as.binned}{list}(l, breaks=NULL)
}
\arguments{
  \item{counts}{the numerical vector of bin contents}
  \item{bins}{the numerical vector of bin centers}
  \item{breaks}{an optional numerical vector defining the egdes of each bin.}
  \item{h}{an object of type \code{"histogram"} usually, a result of a
        call to \code{\link{hist}}.}
  \item{d}{a data.frame containing at least one colum, see \sQuote{Details}.}
  \item{l}{a list of vectors, see \sQuote{Details}.}
}
\details{
  TBW
}
\value{
  \code{binned} and \code{as.binned} return a S3 object of type
  \code{"binned"} which is a data.frame with 2 columns: \code{"counts"} 
  and \code{"bins"} representing the binned data and the bin center.
  Optionally the attribute \code{"breaks"} may be present.
}
\examples{
  x <- -50:50
  b <- binned(sapply(x, function(i) rpois(1, 1000 * dnorm(i, 0, 10))), x)
  summary(b)

  b1 <- as.binned(hist(rnorm(1000, 0, 10), breaks=-50:50, plot=F))
  summary(b1)
}
