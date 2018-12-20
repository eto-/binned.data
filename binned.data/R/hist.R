# histograms
hist <- function (x=NULL, breaks=NULL, bin.width = 0, right=F, plot=T, ...) {
  xname <- paste(deparse(substitute(x), 500), collapse="\n")

  if (!length(x) && !length(breaks)) stop("data or breaks should be specified")
  if (length(breaks) > 0 && bin.width > 0) stop("breaks and bin.width are mutually exclusive")

  h_ <- NULL
  if (length(x) > 0) {
    if (bin.width > 0) h_ <- hist (x, breaks=seq(floor(min(x)/bin.width)*bin.width, ceiling(max(x)/bin.width)*bin.width, bin.width), right=right, plot=F, ...)
    else if (!length(breaks)) h_ <- hist.default (x, right=right, plot=F, ...)
    else if (length(breaks) == 1) h_ <- hist.default (x, breaks=breaks, right=right, plot=F, ...)
    else if (right) h_ <- hist.default (x[x > breaks[1] & x <= tail(breaks, 1)], breaks=breaks, include.lowest=F, right=right, plot=F, ...) 
    else h_ <- hist.default (x[x >= breaks[1] & x < tail(breaks, 1)], breaks=breaks, include.lowest=F, right=right, plot=F, ...)
  } else {
    h_ <- hist.default (breaks[1], breaks=breaks, plot=F, ...)
    h_$density[1] <- 0
    h_$counts[1] <- 0
  }

  h_$xname <- xname

#  h_$right <- right

  if(plot) plot(h_)

  invisible(h_)
}
