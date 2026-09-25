###############################################################################
# EBC4008 Case 3 -- function library.
###############################################################################


## ---- 1. Infrastructure ----------------------------------------------------

# Paths + switches in one list. Creates the output dirs.
case3_config <- function(root = ".",
                         cache_on  = TRUE,
                         cache_tag = "v1",
                         save_figs = TRUE,
                         seed      = 1) {
  cfg <- list(
    root      = root,
    datadir   = file.path(root, "data"),
    resdir    = file.path(root, "results"),
    figdir    = file.path(root, "results", "figures"),
    tabdir    = file.path(root, "results", "tables"),
    cachedir  = file.path(root, "results", "cache"),
    logdir    = file.path(root, "results", "log"),
    cache_on  = cache_on,
    cache_tag = cache_tag,
    save_figs = save_figs,
    seed      = seed
  )
  for (d in c(cfg$datadir, cfg$resdir, cfg$figdir, cfg$tabdir, cfg$cachedir, cfg$logdir))
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  cfg
}

# cache("key", { slow })
# editing the body does NOT invalidate the cache -- bump cache_tag or delete
# results/cache when you change a model.
# The block runs in its own environment: scratch variables in it don't leak
# out. They used to, on a miss only -- so a cold run could overwrite globals
# (section 10's `s2` clobbered section 4's ADF object) and a warm run couldn't.
cache <- function(key, expr, force = FALSE, cfg = CFG) {
  path <- file.path(cfg$cachedir, paste0(key, "__", cfg$cache_tag, ".rds"))
  if (isTRUE(cfg$cache_on) && !force && file.exists(path)) {
    message(sprintf("  [cache HIT ] %s", key))
    return(readRDS(path))      # returns before touching expr
  }
  message(sprintf("  [cache MISS] %s -- computing", key))
  t0  <- Sys.time()
  val <- eval(substitute(expr), new.env(parent = parent.frame()))
  el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (isTRUE(cfg$cache_on)) saveRDS(val, path)
  message(sprintf("  [cache SAVE] %s (%.1f s)", key, el))
  val
}

# fig("name", { plot(...) }) -> results/figures/name.png
# Rscript has no screen device; without this, plots vanish into Rplots.pdf.
# save_figs = FALSE draws to the active device instead.
# par() gets reset here, so plot code can just call par(mfrow = ...). Don't put
# on.exit() in the { } -- it runs at top level, where Rscript ignores it and
# source() fires it after dev.off().
fig <- function(name, expr, width = 1500, height = 950, res = 130, cfg = CFG) {
  if (isTRUE(cfg$save_figs)) {
    path <- file.path(cfg$figdir, paste0(name, ".png"))
    grDevices::png(path, width = width, height = height, res = res)
    on.exit({
      grDevices::dev.off()
      message(sprintf("  [fig] %s", basename(path)))
    }, add = TRUE)
  } else {
    op <- graphics::par(no.readonly = TRUE)     # png: par dies with the device
    on.exit(graphics::par(op), add = TRUE)
  }
  force(expr)
  invisible(NULL)
}

# Print a table and write the csv. the csv keeps full precision.
tab <- function(name, x, digits = 4, caption = NULL, cfg = CFG) {
  if (!is.null(caption)) cat("\n", caption, "\n", sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.numeric, logical(1))
  pr  <- df
  pr[num] <- lapply(pr[num], function(z) round(z, digits))
  print(pr)
  utils::write.csv(df, file.path(cfg$tabdir, paste0(name, ".csv")), row.names = TRUE)
  message(sprintf("  [tab] %s.csv", name))
  invisible(x)
}

# Log banners.
say <- function(...) {
  msg <- paste0(...)
  bar <- strrep("=", max(78, nchar(msg) + 4))
  cat("\n", bar, "\n== ", msg, "\n", bar, "\n", sep = "")
}

note <- function(...) cat("\n-- ", paste0(...), "\n", sep = "")


## ---- 2. Data ---------------------------------------------------------------

# getSymbols with retries. Yahoo rate-limits and quantmod has no backoff.
# auto.assign = FALSE: the default assigns to globalenv and returns only the
# symbol name, which breaks lapply().
safe_get <- function(sym, from, to, tries = 3L, base_sleep = 2) {
  for (i in seq_len(tries)) {
    x <- tryCatch(
      quantmod::getSymbols(sym, src = "yahoo", from = from, to = to,
                           auto.assign = FALSE),
      error = function(e) {
        message(sprintf("  %s attempt %d/%d failed: %s", sym, i, tries,
                        conditionMessage(e)))
        NULL
      })
    # Yahoo can return a valid but EMPTY object when throttled, so !is.null
    # alone would let a 0-row series into the panel.
    if (!is.null(x) && NROW(x) > 0L) return(x)
    if (!is.null(x)) message(sprintf("  %s attempt %d/%d returned 0 rows", sym, i, tries))
    if (i < tries) Sys.sleep(base_sleep^i)      # no sleep after the last try
  }
  stop("failed to download ", sym, " after ", tries, " attempts")
}

# One row of integrity counts + range diagnostics per ticker.
# The Parkinson range is quadratic in log(H/L), so one bad tick
# wrecks everything downstream.
screen_quality <- function(x) {
  H <- as.numeric(quantmod::Hi(x)); L <- as.numeric(quantmod::Lo(x))
  C <- as.numeric(quantmod::Cl(x)); O <- as.numeric(quantmod::Op(x))
  V <- as.numeric(quantmod::Vo(x)); A <- as.numeric(quantmod::Ad(x))
  r <- diff(log(A))                     # adjusted c2c -- the return we model
  n <- length(C)
  # Decomposition below is an exact identity, so it needs RAW closes:
  #   log(C_t/C_{t-1}) = log(O_t/C_{t-1}) + log(C_t/O_t)
  rcc  <- log(C[-1] / C[-n])            # raw close-to-close, rows 2..n
  onr  <- log(O[-1] / C[-n])            # overnight
  intr <- log(C[-1] / O[-1])            # open-to-close
  vt   <- stats::var(rcc, na.rm = TRUE)
  data.frame(
    n            = length(H),
    na_any       = sum(!stats::complete.cases(cbind(O, H, L, C, V, A))),
    na_HL        = sum(is.na(H) | is.na(L)),
    zero_range   = sum(H == L, na.rm = TRUE),          # -> RR = 0 -> log(0)
    H_lt_L       = sum(H < L, na.rm = TRUE),           # Must be impossible
    C_outside_HL = sum(C > H | C < L, na.rm = TRUE),   # after-hours prints
    O_outside_HL = sum(O > H | O < L, na.rm = TRUE),
    zero_volume  = sum(V == 0, na.rm = TRUE),
    nonpos_price = sum(c(O, H, L, C, A) <= 0, na.rm = TRUE),
    max_logHL    = round(max(log(H / L), na.rm = TRUE), 4),
    max_absret   = round(max(abs(r), na.rm = TRUE), 4),
    # Under GBM, E[log(H/L)] / E|return| = 2 exactly, PROVIDED both cover the
    # same clock time. _cc divides by the close-to-close return, which also
    # spans the overnight gap -> biased to ~1.7. _oc uses open-to-close, the
    # same window as H/L -> lands on 2. Use _oc; _cc is kept for contrast.
    range_ratio_cc = round(mean(log(H / L), na.rm = TRUE) /
                           mean(abs(r), na.rm = TRUE), 3),
    range_ratio_oc = round(mean(log(H / L), na.rm = TRUE) /
                           mean(abs(log(C / O)), na.rm = TRUE), 3),
    # var(total) = var(overnight) + var(intraday) + 2cov, so these sum to 1.
    # Don't use 1 - var(intraday)/var(total) as the overnight share -- that
    # silently absorbs the covariance term.
    intraday_share  = round(stats::var(intr, na.rm = TRUE) / vt, 3),
    overnight_share = round(stats::var(onr,  na.rm = TRUE) / vt, 3),
    cov_share       = round(2 * stats::cov(onr, intr, use = "complete.obs") / vt, 3),
    # RR measures intraday variance, GARCH measures total -> GARCH should read
    # this much higher. Sections 8 and 10 test it against an ML fit.
    pred_garch_over_RR = round(1 / (stats::var(intr, na.rm = TRUE) / vt), 3),
    # Parkinson vs squared open-to-close return: same target, expected ~1.
    RR_over_intraday_var = round(
      mean((log(H / L))^2 / (4 * log(2)), na.rm = TRUE) /
      mean(log(C / O)^2, na.rm = TRUE), 3)
  )
}

# Stop the run if the data can't support a range-based analysis.
# zero_range: RR = 0 -> log(RR) = -Inf -> poisons HAR and the VAR.
assert_quality <- function(q) {
  stopifnot(all(q$zero_range == 0))
  stopifnot(all(q$H_lt_L == 0))
  stopifnot(all(q$na_HL == 0))
  stopifnot(all(q$nonpos_price == 0))
  s <- q$intraday_share + q$overnight_share + q$cov_share
  stopifnot(all(abs(s - 1) < 0.01))        # tolerance: the columns are rounded
  invisible(TRUE)
}


## ---- 3. Plotting -----------------------------------------------------------

# Stack k series in k panels.
# plot.xts ignores par(mfrow) -- it has its own layout engine,
# plot.zoo respects mfrow.
plot_panels <- function(series, titles = names(series), ylabs = NULL,
                        main = NULL, col = "grey20", ...) {
  k <- length(series)
  op <- graphics::par(mfrow = c(k, 1),
                      mar = c(2.4, 4.0, 2.0, 1.0),
                      oma = c(0, 0, if (is.null(main)) 0 else 2.2, 0))
  on.exit(graphics::par(op), add = TRUE)      # or the margins leak into later figures
  for (i in seq_len(k)) {
    graphics::plot(zoo::as.zoo(series[[i]]), xlab = "", col = col,
                   ylab = if (is.null(ylabs)) "" else ylabs[i],
                   main = titles[i], ...)
    graphics::grid(col = "grey88")
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, cex = 1.05, font = 2)
  invisible(NULL)
}

# Margins for a grid of acf()/pacf() panels. plot.acf puts its title 2 lines
# above the plot, so a top margin under ~3.5 pushes the title into the panel
# above (onto its "Lag" label).
ACF_MAR <- c(4, 4, 3.8, 1)

# acf() / pacf() plus the ARCH-robust 95% band in red (the dashed blue one is
# iid). Lag 0 (always 1) is dropped so the small lags and the band are readable.
# The PACF uses the ACF's robust se (the two coincide for a nearly white series).
# Uses robust_acf_se() from section 4.
acf_robust <- function(x, lag.max = 40, main = NULL, ...) {
  ra <- robust_acf_se(x, K = lag.max)
  a  <- stats::acf(as.numeric(x), lag.max = lag.max, plot = FALSE)
  a$acf <- a$acf[-1, , , drop = FALSE]; a$lag <- a$lag[-1, , , drop = FALSE]
  plot(a, main = main, ylim = range(a$acf, 2 * ra$se_robust, -2 * ra$se_robust), ...)
  graphics::lines(ra$lag,  1.96 * ra$se_robust, col = "red", lwd = 1.6)
  graphics::lines(ra$lag, -1.96 * ra$se_robust, col = "red", lwd = 1.6)
  invisible(ra)
}
pacf_robust <- function(x, lag.max = 40, main = NULL, ...) {
  ra <- robust_acf_se(x, K = lag.max)
  p  <- stats::pacf(as.numeric(x), lag.max = lag.max, plot = FALSE)
  plot(p, main = main, ylim = range(p$acf, 2 * ra$se_robust, -2 * ra$se_robust), ...)
  graphics::lines(ra$lag,  1.96 * ra$se_robust, col = "red", lwd = 1.6)
  graphics::lines(ra$lag, -1.96 * ra$se_robust, col = "red", lwd = 1.6)
  invisible(ra)
}

# Grey bands behind a plot. Call AFTER plotting, then redraw the line with
# lines() and restore the frame with box() -- rect() paints over both.
# Either pass xfun = as.POSIXct or convert the series with as_date_xts() first.
shade_episodes <- function(episodes, col = "grey60", alpha = 0.28, border = NA,
                           xfun = as.Date) {
  usr <- graphics::par("usr")
  for (nm in names(episodes)) {
    w <- episodes[[nm]]
    graphics::rect(xfun(w[1]), usr[3], xfun(w[2]), usr[4],
                   col = grDevices::adjustcolor(col, alpha), border = border)
  }
  invisible(NULL)
}

# Force a Date index. Convert once at the boundary rather than tracking units
# at every call site.
as_date_xts <- function(x) {
  idx <- zoo::index(x)
  if (inherits(idx, "Date")) return(x)
  xts::xts(zoo::coredata(x), order.by = as.Date(idx))
}

# Top n values with their dates
top_by_value <- function(x, n = 10, decreasing = TRUE, value_name = "value") {
  stopifnot(NCOL(x) == 1)
  v <- as.numeric(zoo::coredata(x))
  o <- order(v, decreasing = decreasing)[seq_len(min(n, length(v)))]
  out <- data.frame(rank = seq_along(o),
                    date = format(zoo::index(x)[o]),
                    value = v[o], row.names = NULL)
  names(out)[3] <- value_name
  out
}


## ---- 4. ACF / test utilities ----------------------------------------------

# Produce correct standard errors for the sample autocorrelations
#
# Are returns really autocorrelated, or does it only look that way because 
# volatility comes in clusters?
#
# ARCH-robust standard errors for the sample autocorrelations.
# The band acf() draws is 1.96/sqrt(n), which assumes iid. Returns are
# uncorrelated but heteroskedastic, so the real se is bigger and the iid band
# flags lags as significant that aren't. Diebold (1986), LNS (2001):
# var(rho_k) ~ (1/n) * mean(x_t^2 x_{t-k}^2) / mean(x_t^2)^2
# which collapses to 1/n if the squares are unrelated. 
# `inflation` is the ratio of the two, i.e. a direct read on volatility clustering.
robust_acf_se <- function(x, K = 40) {      
  x <- as.numeric(x) - mean(as.numeric(x))
  n <- length(x)
  d <- mean(x^2)^2
  se <- vapply(seq_len(K), function(k) {
    a <- x[(k + 1):n]; b <- x[1:(n - k)]     # a[i] and b[i] are k apart
    sqrt(mean(a^2 * b^2) / d / n)           # robust se of rho_k
  }, numeric(1))
  data.frame(lag = seq_len(K), acf = as.numeric(stats::acf(x, lag.max = K,
                                                           plot = FALSE)$acf)[-1],
             se_iid = 1 / sqrt(n), se_robust = se,
             inflation = se * sqrt(n))       
}

# Bundle the first m autocorrelations into one joint test of 
# "is the conditional mean predictable at all", using the corrected standard 
# errors rather than the iid ones.
# 
# Ljung-Box, textbook and ARCH-robust.
# The standard version uses the same iid se as the acf() band, so it
# over-rejects on returns: it detects heteroskedasticity and calls it mean predictability
lb_robust <- function(x, lags = 10) {
  x  <- as.numeric(x)
  n  <- length(x)
  ra <- robust_acf_se(x, K = lags)
  q_std <- n * (n + 2) * sum(ra$acf^2 / (n - ra$lag))
  q_rob <- sum((ra$acf / ra$se_robust)^2)
  data.frame(lags = lags,
             LB_standard = q_std,
             p_standard  = stats::pchisq(q_std, lags, lower.tail = FALSE),
             LB_robust   = q_rob,
             p_robust    = stats::pchisq(q_rob, lags, lower.tail = FALSE),
             mean_se_inflation = mean(ra$inflation))
}

# Compute R² on the variance scale by its actual definition, 
# rather than via var(), that quietly deletes part of the error.
# R^2 via sums of squares.
# 1 - var(resid)/var(y) is wrong: var() de-means the residuals, deleting error
# the model didn't explain, so the number comes out too high.
pseudo_r2 <- function(resid, y) {
  resid <- as.numeric(resid); y <- as.numeric(y)
  1 - sum(resid^2) / sum((y - mean(y))^2)
}


## ---- 5. Realized measures --------------------------------------------------

# Parkinson (1980) daily realized range: (log H - log L)^2 / (4 log 2).
# Uses RAW High/Low on purpose, any split/dividend factor multiplies both and
# cancels in the ratio, so the adjusted series is not needed here.
# 4*log(2) is the constant; with only H and L (m = 1) and
# open-to-close coverage, RR sits BELOW a close-to-close variance. That gap is
# the overnight return, quantified in screen_quality().
# H == L -> RR = 0 -> log(0) = -Inf, so those days become NA and are dropped
realized_range <- function(x, name = "RR") {
  rr <- (log(quantmod::Hi(x)) - log(quantmod::Lo(x)))^2 / (4 * log(2))
  bad <- !is.finite(rr) | rr <= 0
  rr[bad] <- NA
  colnames(rr) <- name
  rr
}

# Monthly RV = SUM of squared daily log returns (not mean, quadratic
# variation is additive over time). No demeaning: demeaning adds its own bias.
# Also returns the French-Schwert-Stambaugh correction. FSS fixes
# non-synchronous trading in INDEX returns, for one liquid stock the adjacent
# autocovariance is ~0 or negative (bid-ask bounce), so it can go negative.
# Use the plain version.
# min_obs filters BOTH rv and n_obs -- the sqrt(2/N) precision figure is
# computed from n_obs, so a mismatched mask makes it wrong.
# Stamped to month END (last date the measure is knowable).
monthly_rv <- function(r, min_obs = 15L, stamp = c("end", "start")) {
  stamp <- match.arg(stamp)
  rn  <- as.numeric(r)
  key <- format(zoo::index(r), "%Y-%m")

  rv    <- tapply(rn, key, function(z) sum(z^2))
  fss   <- tapply(rn, key, function(z) sum(z^2) + 2 * sum(z[-length(z)] * z[-1]))
  n_obs <- tapply(rn, key, length)
  dmin  <- tapply(zoo::index(r), key, min)
  dmax  <- tapply(zoo::index(r), key, max)

  keep  <- n_obs >= min_obs
  rv <- rv[keep]; fss <- fss[keep]; n_obs <- n_obs[keep]
  dmin <- dmin[keep]; dmax <- dmax[keep]

  when <- as.Date(if (stamp == "end") dmax else dmin, origin = "1970-01-01")
  xts::xts(cbind(RVm = as.numeric(rv),
                 RVm_fss = as.numeric(fss),
                 nobs = as.numeric(n_obs)),
           order.by = when)
}


## ---- 6. HAR ----------------------------------------------------------------

# HAR design (Corsi 2009): y = v[t+1], regressors = trailing means of v over
# each horizon in p. Averages not sums, so the coefficients stay comparable.
# sides = 1 is load-bearing: the default sides = 2 is CENTRED and would put
# future volatility on the RHS, giving a spurious R^2.
# Row i -> regressors at date[i + max(p) - 1], target at date[i + max(p)].
# na.omit drops max(p)-1 rows at the front and 1 at the back. Use
# har_fitted_dates() rather than hand-writing -(1:22); a one-day shift inflates
# the fit and is invisible in the output.
# check_gaps: filter() runs positionally and ignores dates, so if rows were
# dropped a "22-day" average can span more than 22 calendar days.
har_design <- function(v, p = c(1, 5, 22), check_gaps = TRUE) {
  stopifnot(is.numeric(p), all(p >= 1), !anyDuplicated(p))
  p  <- sort(as.integer(p))
  vv <- as.numeric(v)
  n  <- length(vv)
  stopifnot(n > max(p) + 1L)

  if (check_gaps && inherits(v, c("xts", "zoo"))) {
    idx <- zoo::index(v)
    if (inherits(idx, "Date")) {
      gaps <- as.numeric(diff(idx))
      nbig <- sum(gaps > 4)                 # >4 days = more than a long weekend
      if (nbig > 0)
        message(sprintf(paste0("  [har_design] %d gap(s) > 4 calendar days in the ",
                               "index; trailing averages span them positionally."), nbig))
    }
  }

  X <- vapply(p, function(k) as.numeric(stats::filter(vv, rep(1 / k, k), sides = 1)),
              numeric(n))
  colnames(X) <- paste0("h", p)
  d <- data.frame(y = c(vv[-1], NA_real_), X)
  stats::na.omit(d)
}

# Dates the HAR fitted values belong to. Derived so the caller can assert it.
har_fitted_dates <- function(v, p = c(1, 5, 22)) {
  idx <- if (inherits(v, c("xts", "zoo"))) zoo::index(v) else seq_along(v)
  m   <- max(as.integer(p))
  idx[(m + 1L):length(idx)]
}

# Fit HAR with Newey-West.
# HAC is needed because the residuals are autocorrelated
# NOT because the regressors overlap, which does not by itself induce it in a 
# one-step-ahead regression. 
# prewhite = FALSE: NeweyWest defaults to TRUE, which is Andrews-Monahan
# (1992), not Newey-West (1987).
# jensen: we model logs, so exp(fitted) is a MEDIAN. Set TRUE for anything
# reported as a forecast.
fit_har <- function(v, p = c(1, 5, 22), hac_lag = max(p), jensen = FALSE) {
  d   <- har_design(v, p = p)
  rhs <- paste(paste0("h", sort(as.integer(p))), collapse = " + ")
  f   <- stats::as.formula(paste("y ~", rhs))
  m   <- stats::lm(f, data = d)

  V   <- sandwich::NeweyWest(m, lag = hac_lag, prewhite = FALSE, adjust = TRUE)
  ct  <- lmtest::coeftest(m, vcov. = V)

  dts <- har_fitted_dates(v, p)
  stopifnot(length(stats::fitted(m)) == length(dts))     # alignment

  s2  <- sum(stats::residuals(m)^2) / stats::df.residual(m)
  lvl <- exp(stats::fitted(m) + if (jensen) s2 / 2 else 0)

  list(fit         = m,
       design      = d,
       hac         = ct,
       vcov_hac    = V,
       r2          = summary(m)$r.squared,
       persistence = sum(stats::coef(m)[-1]),
       sigma2      = s2,
       fitted_log  = xts::xts(stats::fitted(m), order.by = dts),
       fitted_lvl  = xts::xts(lvl, order.by = dts),
       p           = sort(as.integer(p)))
}


## ---- 7. Forecast evaluation ------------------------------------------------

# QLIKE loss, minimised at f = proxy. Robust to noise in the volatility proxy
# (Patton 2011) and much less dominated by crisis outliers than MSE.
# Needs strictly positive inputs; anything else comes back NA.
qlike <- function(proxy, f) {
  ok <- is.finite(proxy) & is.finite(f) & proxy > 0 & f > 0
  out <- rep(NA_real_, length(proxy))
  z   <- proxy[ok] / f[ok]
  out[ok] <- z - log(z) - 1
  out
}

mse_loss <- function(proxy, f) (proxy - f)^2

# Mincer-Zarnowitz: proxy = a + b*f, HAC, plus a joint Wald test of (a,b)=(0,1).
# b < 1 = the forecast swings too much (or is noisy); b > 1 = too timid.
# A multiplicative level error (f = c*truth) ALSO moves b (to 1/c), so rescale
# before reading b. a != 0 with b ~ 1 = additive bias.
mz_test <- function(proxy, f, hac_lag = 22) {
  ok <- is.finite(proxy) & is.finite(f)
  d  <- data.frame(y = proxy[ok], f = f[ok])
  m  <- stats::lm(y ~ f, data = d)
  V  <- sandwich::NeweyWest(m, lag = hac_lag, prewhite = FALSE, adjust = TRUE)
  ct <- lmtest::coeftest(m, vcov. = V)
  W  <- car_free_wald(m, V, c(0, 1))
  list(fit = m, hac = ct, r2 = summary(m)$r.squared,
       a = stats::coef(m)[1], b = stats::coef(m)[2],
       wald_stat = W$stat, wald_df = W$df, wald_p = W$p, n = nrow(d))
}

# Wald test of coef == value. Not worth a dependency
# on car just for linearHypothesis().
car_free_wald <- function(m, V, value) {
  b <- stats::coef(m)
  stopifnot(length(b) == length(value))
  d <- as.numeric(b - value)
  W <- as.numeric(t(d) %*% solve(V) %*% d)
  k <- length(b)
  list(stat = W, df = k, p = stats::pchisq(W, df = k, lower.tail = FALSE))
}

# Forecast encompassing: proxy = a + bG*f1 + bH*f2, HAC.
# Asks what a GARCH-X would ask (does the range add anything beyond GARCH?) but
# with OLS. Used because rugarch's variance-equation external.regressors left
# its coefficient pinned at exactly 0 with t = 0.
# vif is returned so the collinearity is visible.
encompass_test <- function(proxy, f1, f2, names2 = c("f1", "f2"), hac_lag = 22) {
  ok <- is.finite(proxy) & is.finite(f1) & is.finite(f2)
  d  <- data.frame(y = proxy[ok], a = f1[ok], b = f2[ok])
  m  <- stats::lm(y ~ a + b, data = d)
  V  <- sandwich::NeweyWest(m, lag = hac_lag, prewhite = FALSE, adjust = TRUE)
  ct <- lmtest::coeftest(m, vcov. = V)
  rho <- stats::cor(d$a, d$b)
  list(fit = m, hac = ct, r2 = summary(m)$r.squared,
       coef = stats::setNames(stats::coef(m), c("(Intercept)", names2)),
       p_drop_f2 = ct[3, 4],          # H0: coef on f2 is 0
       p_drop_f1 = ct[2, 4],          # H0: coef on f1 is 0
       corr_f1_f2 = rho,
       vif = 1 / (1 - rho^2),
       combined = stats::fitted(m),
       n = nrow(d))
}

# Diebold-Mariano on a loss differential, with the HLN small-sample correction.
# Hand-rolled because forecast::dm.test takes forecast ERRORS and applies its
# own power loss; QLIKE isn't a power of the error, so that would test the
# wrong hypothesis.
# Tells you which of two forecasts is better, not whether either is any good 
dm_test <- function(l1, l2, h = 1, hln = TRUE) {
  d  <- l1 - l2
  d  <- d[is.finite(d)]
  n  <- length(d)
  dbar <- mean(d)
  g0 <- stats::var(d) * (n - 1) / n          # Newey-West lrv, Bartlett, trunc h-1
  lr <- g0
  if (h > 1) for (k in seq_len(h - 1)) {
    gk <- sum((d[-seq_len(k)] - dbar) * (utils::head(d, -k) - dbar)) / n
    lr <- lr + 2 * (1 - k / h) * gk
  }
  stat <- dbar / sqrt(lr / n)
  if (hln && n > h) {
    corr <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
    stat <- stat * corr
  }
  list(stat = stat, n = n, mean_diff = dbar,
       p = 2 * stats::pt(-abs(stat), df = n - 1))
}


## ---- 8. VAR / spillovers ---------------------------------------------------

# Biggest p that leaves min_obs_per_param observations per coefficient.
# A VAR(p) on N series costs N*p + 1 params PER EQUATION out of T - p rows.
# Without a cap VARselect runs off a singular residual covariance 
# (criterion -> -Inf), picks the most overfitted order in the grid, and returns 
# a spillover index built on noise. Returns 0 when nothing is estimable
var_max_lag <- function(Tn, N, lag_max = 8L, min_obs_per_param = 10) {
  feasible <- vapply(seq_len(lag_max), function(k) {
    (N * k + 1) * min_obs_per_param <= (Tn - k)
  }, logical(1))
  if (!any(feasible)) return(0L)
  max(which(feasible))
}

# VARselect capped at what's estimable.
# unname(): VARselect()$selection returns a NAMED integer, and the name rides
# along into downstream data.frames as confusing row labels.
var_select_capped <- function(y, lag_max = 8L, ic = "SC(n)", type = "const",
                              min_obs_per_param = 10) {
  Tn <- nrow(y); N <- ncol(y)
  cap <- var_max_lag(Tn, N, lag_max = lag_max, min_obs_per_param = min_obs_per_param)
  if (cap == 0L) return(list(p = 0L, cap = 0L, ic_choice = NA_integer_, feasible = FALSE))
  sel <- vars::VARselect(y, lag.max = cap, type = type)$selection
  p   <- unname(as.integer(sel[ic]))
  list(p = p, cap = cap, ic_choice = p, feasible = TRUE, selection = sel)
}


## ---- 9. Formatting ---------------------------------------------------------

stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))

# urca S4 object -> one tidy row (stat next to its critical values).
# @teststat can hold several stats (tau, phi1, phi2, phi3); `which` picks one.
urca_row <- function(obj, label, which = 1L) {
  ts  <- obj@teststat
  cv  <- obj@cval
  nm  <- if (!is.null(rownames(cv))) rownames(cv)[min(which, nrow(cv))] else NA_character_
  data.frame(test = label, statistic = as.numeric(ts)[which], stat_name = nm,
             cv_1pct = cv[min(which, nrow(cv)), "1pct"],
             cv_5pct = cv[min(which, nrow(cv)), "5pct"],
             cv_10pct = cv[min(which, nrow(cv)), "10pct"],
             row.names = NULL)
}

# Crisis windows, each anchored to a published event date and fixed HERE --
# before the rolling spillover index is plotted. Choosing breakpoints off the
# figure you're about to interpret is circular.
# covid and banks23 are only ~94 days: var_max_lag() caps them at p = 1.
# Estimable, but treat those per-episode numbers as rough.
case3_episodes <- function() {
  list(
    dotcom  = c("2000-03-24", "2002-10-09"),  # S&P 500 peak to trough
    gfc     = c("2007-08-01", "2009-06-30"),  # quant quake / BNP freeze
    euro    = c("2011-07-01", "2012-07-26"),  # euro stress 
    covid   = c("2020-02-19", "2020-06-30"),  # pre-covid S&P peak
    banks23 = c("2023-02-15", "2023-06-30")   # SVB/Signature 
  )
}

# Seam-free crisis benchmark. The pooled crisis sample splices 5 non-contiguous
# blocks, which breaks the time ordering at each seam; GFC alone (483
# contiguous days) is the cross-check.
case3_clean_crisis <- function() list(gfc = c("2007-08-01", "2009-06-30"))
