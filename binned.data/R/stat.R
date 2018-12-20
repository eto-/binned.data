.effective.variance <- function (pdf, p, mids, mids.variance, counts.variance) {
  g <- grad(function(x) pdf(x, p), mids)

  counts.variance + mids.variance * g**2
}

chi2.uncertainties <- function(pdf, p, mids, counts, mids.se=c(), counts.se) {
  if (length(counts.se) == 1) {
    if (counts.se < 0) counts.se <- counts * - counts.se
    else counts.se = rep(counts.se, length(counts))
  }

  variance <- counts.se**2

  if (length(mids.se) > 0) {
    if (length(mids.se) == 1) mids.se = rep(mids.se, length(mids))
    
    variance <- .effective.variance (pdf, p, mids, mids.se**2, variance)
  }

  (pdf(mids, p) - counts)**2 / variance
}

chi2.neyman <- function (pdf, p, mids, counts, mids.se=c(), ...) chi2.uncertainties (pdf, p, mids, counts, mids.se, sqrt(counts))
chi2.pearson <- function (pdf, p, mids, counts, mids.se=c(), ...) chi2.uncertainties (pdf, p, mids, counts, mids.se, sqrt(pdf(mids, p)))
chi2.saturared.poisson <- function (pdf, p, mids, counts, ...) 2 * (pdf(mids, p) - counts) + chi2.saturared.multinomial(pdf, p, mids, counts)
chi2.saturared.multinomial <- function (pdf, p, mids, counts, ...) ifelse (counts > 0, 2 * counts * log (counts / pdf(mids, p)), 0)

# negative emg
dnemg <- function(x, mu=0, ...) demg(-x, -mu, ...)

# logexp distribution
dlogexp <- function(x, m, base=10) { k <- base**x; k / m * exp(-k / m) * log(base) }
rlogexp <- function(n, m, base=10) log(rexp(n, 1/m), base=base)

