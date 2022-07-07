# norm.test fitter
fitter <- function (d, pdf, start, fixed=list(), counts.se=c(), bins.se=c(),
		    statistics=c("gaussian", "neyman", "pearson", "poissonian", "multinomial", "uncertainties"), 
		    method=c("NLS", "MIGRAD", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), maxiter=-1, na.action) {

  statistics <- tolower(statistics)
  chi2.base <- switch(match.arg(statistics),
		      gaussian=, neyman= chi2.neyman,
		      pearson= chi2.pearson,
		      poissonian= chi2.saturared.poisson,
		      multinomial= chi2.saturared.multinomial,
		      uncertainties = chi2.uncertainties)

  method <- match.arg(method)

  f <- "fit"
  s <- "summary"
  df <- numeric()
  chi2 <- numeric()
  r.c <- numeric() 

  d <- as.binned(d)

  if (identical(chi2.base, chi2.neyman)) d <- subset(d, counts > 0)

  if (identical(chi2.base, chi2.uncertainties)) {
    if (length(counts.se) != 1 && length(counts.se) != nrow(d)) stop("counts uncertainty can be 1 element or an array long as the observations")
  } else if (length(counts.se)) stop(paste("counts uncertainty specified with", statistics, "model"))

  if (length(bins.se) > 0) {
    if (identical(chi2.base, chi2.saturared.poisson,) || identical(chi2.base, chi2.saturared.multinomial))
      stop("mids uncertainties are not implemented for poissonian and multinomial regressions")
    require(numDeriv)
  }

  if (missing(na.action)) na.action <- match.fun(unlist(options("na.action")))

  if (method == "NLS") {
    chi.fn <- function (p) { 
      if (length(fixed) > 0) p[names(fixed)] <- fixed; 
      r <- sqrt(chi2.base(pdf, p, d$bins, d$counts, bins.se, counts.se)) 
      na.action(r)
    }

    control = nls.control()
    if (maxiter > 0) control = nls.control(maxiter=maxiter)

    f <- nls.lm (par=start, fn=chi.fn, control=control)
    s <- summary(f)

    return(list(r.c=coef(s), chi2=deviance(f), df=df.residual(f), f=f, s=s))
  } else {
    require(bbmle)
    if (method == "MIGRAD") require(migrad)

    if(length(start) == 1 && length(fixed) == 0) fixed = list(fix.mle2.bug = 0) # for 1 parameter only mle2 forgets the name

    mll <- function(v) sum(na.action(chi2.base(pdf, as.list(v), d$bins, d$counts, bins.se, counts.se)))/2
    parnames(mll) <- c(names(start), names(fixed))

    f <- switch(method,
      MIGRAD= mle2(mll, start=start, vecpar=T, fixed=fixed, optimfun=migrad, optimizer="user", browse_obj=F),
      mle2(mll, start=start, vecpar=T, fixed=fixed, method=method))

    s <- bbmle::summary(f)

    return(list(r.c=s@coef, chi2=2 * f@min, df=nrow(d) - length(f@coef), f=f, s=s))
  }
}

# norm.test
norm.test <- function (x, core=c("gauss", "poisson", "logexp", "cauchy", "emg", "nemg", "exp", "user"), 
		       background=c("none", "flat", "linear", "exponential", "step", "gauss"), 
		       background.modifier=c("", "+", "!"), quiet=1, summary=F, plot=F, start=list(), fixed=list(), range=c(), user.pdf=NULL, ...) {

  core <- tolower(core)
  core <- match.arg(core)
  core.pdf <- norm.test.core.pdf(match.arg(core), user.pdf)

  background <- tolower(background)
  if (identical(eval(formals()$background), background)) background <- background[1];
  background <- match.arg(background, several.ok=T)
  background.pdf <- norm.test.background.pdf (background, match.arg(background.modifier))

  pdf <- function (mids, p) core.pdf(mids, p) + background.pdf(mids, p)

  d <- as.binned(x)
  if (length(range) == 2) d <- subset(d, bins > range[1] & bins < range[2])
  
  start <- norm.test.def.parameters(d, core, background, start, fixed)

  l <- fitter (d, pdf, start, fixed, ...)
  r.c <- l$r.c
  df <- l$df
  chi2 <- l$chi2
  f <- l$f
  s <- l$s
  statistics <- l$statistics

  if (summary) print(s)

  par.names <- c(names(start), names(fixed))
  list_val <- function (m, c) {
    v <- m[,c]
    names(v) <- rownames(m)
    as.list(v)
  }
  par.vals <- c(list_val(r.c, 1), fixed)
  l <- list(); l[names(fixed)] <- 0
  par.errs <- c(list_val(r.c, 2), l)

  bin.width = mean(diff(d$bins))
  pvalue <- 1 - pchisq(chi2, df)
  if (pvalue < quiet) {
    cat ("\n", "\tChi2 test on histogram", "\n", "\n", sep="")
    if (core != "user") {
      cat ("best values are mean = ", par.vals[["m"]], "+-", par.errs[["m"]], " sd = ", par.vals[["s"]], "+-", par.errs[["s"]], "\n", sep="")
      cat ("peak area = ", par.vals[["c"]] / bin.width, "+-", par.errs[["c"]] / bin.width, "\n", sep="")
    }
    cat ("chi2 = ", chi2, " df = ", df, " pvalue = ", pvalue, "\n", sep="")
  }

  r <- list(v=unlist(par.vals), e=unlist(par.errs), 
            peak.area = par.vals[["c"]] / bin.width, peak.area.se = par.errs[["c"]] / bin.width, 
	    chi2=chi2, df=df, pvalue=pvalue, 
	    statistics=statistics, background=background, core=core, 
	    pdf=pdf, background.pdf=background.pdf, core.pdf=core.pdf, fit=f, x=x, 
	    resid=d[,2] - pdf(d[,1], par.vals), 
	    range=range)
  class(r) <- "ntest"

  if (plot) plot.ntest(r)

  invisible(r)
}

predict.ntest <- function(object, newdata) {
  if (!inherits(object, "ntest"))
    warning("calling predict.ntest(<fake-ntest-object>) ...")
  if (missing(newdata) || is.null(newdata)) data <- object$x
  else data <- newdata
  data <- as.binned(data)

  return(t$pdf(data$bins, as.list(t$v)))
}

plot.ntest <- function (t, resid=F, components=F, ...) {
  plot.components <- function (x, t) {
    p <- as.list(t$v)
    bks <- t$background
    c.p <- tail(rainbow(length(bks) + 1), length(bks))
    lines(x, t$core.pdf(x, p), col="black", lwd=2)
    invisible(sapply(seq_along(bks), function(i) lines(x, norm.test.background.pdf(bks[i],"")(x, p), col=c.p[i], lwd=2)))
  }

  if (resid) plot (as.binned(t$x)$bins, t$resid, ...)
  else {
    if (inherits(t$x, "histogram")) {
      if(length(t$range) == 2) plot(t$x, xlim=t$range, ...)
      else plot(t$x, ...)
      lines(t$x$mids, t$pdf(t$x$mids, as.list(t$v)), col="red", lwd=2)
      if (components) plot.components(t$x$mids, t)
    } else {
      o <- as.binned(t$x)
      if (length(t$range) == 2) plot(o$bins, o$counts, "h", xlim=t$range, ...)
      else plot(o$bins, o$counts, "h", ...)
      lines(o$bins, t$pdf(o$bins, as.list(t$v)), col="red", lwd=2)
      if (components) plot.components(o$bins, t)
    }
  }
}


# norm.test ancillaries
norm.test.core.pdf <- function (pdf, f) {
  if (grepl("emg$", pdf)) require(emg)
  f <- switch(pdf,
    G=, gauss= function (mids, p) p$c * dnorm(mids, p$m, abs(p$s)),
    P=, poisson= function (mids, p) p$c * dpois(mids, p$m),
    C=, cauchy= function (mids, p) p$c / (((mids - p$m)/p$s)**2 + 1),
    emg= function (mids, p) p$c * demg(mids, p$m, abs(p$s), abs(p$l)),
    nemg= function (mids, p) p$c * dnemg(mids, p$m, abs(p$s), abs(p$l)),
    LE=, logexp= function (mids, p) p$c * do.call(dlogexp, c(list(x=mids), p[names(p) != "c"])),
    exp= function (mids, p) p$c * exp(mids/-abs(p$t)),
    U=, user= f,
    stop ("unknown core pdf"))
  f
}

norm.test.background.pdf <- function (background, modifier) {
  fs <- c()
  z <- function(mids, p) (mids - p$m) / abs(p$s)
  for (i in background) {
    g <- switch(i,
      N=, none= function (mids, p) rep(0, length(mids)),
      F=, flat= function (mids, p) rep(p$f.const, length(mids)),
      L=, linear= function (mids, p) p$l.const + p$l.slope * mids,
      E=, exponential= function (mids, p) p$e.const * exp(mids / p$e.slope),
      S=, step= function (mids, p) p$s.const * pnorm(mids, p$m, abs(p$s), lower.tail=F),
      G=, gauss= function (mids, p) p$g.c * dnorm(mids, p$g.m, abs(p$g.s)),
      stop ("unknown background"))
    if (nchar(modifier)) { 
      f <- switch(toupper(modifier), 
       "+"= function (mids, p) pmax(0, g(mids, p)),
       "!"= function (mids, p) { v <- g(mids, p); ifelse(v < 0, -1e40, v) },
       stop ("unknown background modifier"))
    } else f <- g

    fs <- c(fs, f)
  }

  function (mids, p) {
    r <- sapply(fs, function (x) x(mids, p))
    if (length(fs) == 1) return(r[,1])
    rowSums(r)
  }
}

norm.test.def.parameters <- function (b, core, background, start, fixed) {
  ms <- as.numeric(summary (b))
  par <- list(c=ms[4] * mean(diff(b$bins)), m=ms[2], s=ms[3])
  slope <- (b$counts[1] - tail(b$counts, 1))/(b$bins[1] - tail(b$bins, 1))

  par <- switch(core,
    G=, gauss= par,
    C=, cauchy= par,
    emg=, nemg= c(par, l=par[['s']]),
    P=, poisson= par[1:2],
    LE=, logexp= par[1:2],
    exp= c(par[1],t=unname(par[2])),
    U=, user= list(),
    stop ("unknown core pdf"))

  for (i in background) {
    p <- switch(i,
      N=, none= list(),
      F=, flat= list(f.const = min(b$counts)),
      L=, linear= list(l.const = b$counts[1], l.slope = slope),
      E=, exponential= list(e.const = b$counts[1]/exp(b$bins[1]/slope), e.slope = 1/slope),
      S=, step= list(s.const = min(b$counts)),
      G=, gauss= list(g.c = 1, g.s=1),
      stop ("unknown background"))

    par <- c(par, p)
  }

  par[names(start)] <- start
  par[names(fixed)] <- NULL

  par
}

fit.barrier <- function (..., x=c(...)) { if (all(x)) return(1); NA }
