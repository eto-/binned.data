as.binned <- function(o, ...) UseMethod("as.binned",o)
as.binned.default <- function(o) stop(paste("unknown cast from", class(o)))
as.binned.binned <- function (o) o
as.binned.histogram <- function(h) binned(h$counts, h$mids, h$breaks)

as.binned.list <- function(l, breaks=NULL) {
  if (length(l) < 1) stop("empty list")
  as.binned(as.data.frame(l), breaks) 
}

as.binned.data.frame <- function(d, breaks=NULL) { 
  n <- c("counts", "bins")

  if (ncol(d) < 1) stop("empty data.frame")

  if (ncol(d) == 1 || !any(names(d) %in% n[1])) d$bins = seq_along(d[,1])
  id <- unlist(lapply(n, function(x) which(names(d) == x)))
  if (length(id) != 2) id <- 1:2
  x <- d[,id]
  colnames(x) <- n

  if (is.null(breaks)) {
    if (any(names(d) %in% "breaks")) {
      stop("unimplemented breaks from data.frame")
      breaks=d[,"breaks"]
    } else if (length(v <- unique(diff(x$bins))) == 1) 
      breaks <- c(x$bins - v/2, tail(x$bins, 1) + v/2)
  } else if (length(breaks) != length(x$bins) + 1) stop("inconsistent user supplied breaks")
  attr(x, "breaks") <- breaks

  class(x) <- c("binned", class(x))
  x 
}

binned <- function(counts, bins, breaks=NULL) {
  if (missing(counts)) return(as.binned(data.frame(numeric(0), numeric(0)), breaks))
  if (missing(bins)) bins <- seq_along(counts)
  as.binned(data.frame(counts, bins), breaks)
}

#.is.binned <- function(b) inherits(b, "binned")

summary.binned <- function (b) {
  attach(b)
  values <- c(length=length(counts))
  if (values[1]) {
    entries <- sum(as.numeric(counts))
    m <- sum (as.numeric(counts) * bins) / entries
    s <- sqrt (sum (as.numeric(counts) * (bins - m)^2) / ifelse(entries > 2, (entries - 1), entries))
    values <- c(values, mean=m, stddev=s, entries=entries, max=max(counts), mode=bins[which.max(counts)])
  }
  detach(b)
  class(values) <- c("summaryDefault", "table")
  values
}

to.histogram <- function(b) {
  if (is.null(attr(b, "breaks"))) stop("undefined breaks")

  xname <- paste(deparse(substitute(b), 500), collapse = "\n")

  h <- hist(numeric(0), breaks=attr(b, "breaks"));
  h$counts <- b$counts

  h$xname <- xname

  h
}

plot.binned <- function (b, xlab="bins", ylab="counts", ...) {
  if (is.null(attr(b, "breaks"))) plot(b$bins, b$counts, type="b", xlab=xlab, ylab=ylab, ...)
  else plot(to.histogram(b), main="", xlab=xlab, ylab=ylab, ...)
}

