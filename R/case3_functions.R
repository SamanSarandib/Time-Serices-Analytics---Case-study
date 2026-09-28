###############################################################################
# EBC4008 Case 3 -- function library (sourced by 00_setup.R, used by parts 01-13)
# Needs:  quantmod, xts, zoo, sandwich, lmtest, vars installed; CFG from case3_config()
# Makes:  helpers -- infra, data screen, plots, ACF/LB, RR/RV, HAR, forecast tests, VAR, episodes
# Writes: nothing when sourced; fig() -> results/figures/*.png, tab() -> results/tables/*.csv,
#         cache() -> results/cache/*.rds
###############################################################################


## ---- 1. Infrastructure ----------------------------------------------------

# paths + switches in one list; creates data/ and results/{figures,tables,cache,log}
# -> list(root, datadir, resdir, figdir, tabdir, cachedir, logdir,
#         cache_on, cache_tag, save_figs, seed)
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

# cache("key", { slow }) -> the block's value, saved to results/cache/<key>__<cache_tag>.rds
# editing the block does NOT refit: change the key or bump cache_tag. Block runs in its own env.
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

# fig("name", { plot(...) }) -> results/figures/name.png (save_figs = FALSE: current device)
# par() is reset afterwards, so plot code may set par(mfrow = ...); no on.exit() inside the { }
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

# print x rounded to `digits` + write results/tables/<name>.csv at full precision -> x (invisible)
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

# console banner: ==== line, "== msg", ==== line
say <- function(...) {
  msg <- paste0(...)
  bar <- strrep("=", max(78, nchar(msg) + 4))
  cat("\n", bar, "\n== ", msg, "\n", bar, "\n", sep = "")
}

# console sub-heading: "-- msg"
note <- function(...) cat("\n-- ", paste0(...), "\n", sep = "")


## ---- 2. Data ---------------------------------------------------------------

# quantmod::getSymbols (Yahoo) with retries, sleeping 2, 4, ... s; a 0-row reply counts as a failure
# -> xts of Open/High/Low/Close/Volume/Adjusted for sym; stops after `tries` failures
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
    if (!is.null(x) && NROW(x) > 0L) return(x)
    if (!is.null(x)) message(sprintf("  %s attempt %d/%d returned 0 rows", sym, i, tries))
    if (i < tries) Sys.sleep(base_sleep^i)      # no sleep after the last try
  }
  stop("failed to download ", sym, " after ", tries, " attempts")
}

# integrity + range checks, one OHLCV xts -> 1 row of 02_quality_screen: n | na_any | na_HL |
#  zero_range | H_lt_L | C_outside_HL | O_outside_HL | zero_volume | nonpos_price | max_logHL |
#  max_absret | range_ratio_cc/_oc | *_share (3) | pred_garch_over_RR | RR_over_intraday_var
screen_quality <- function(x) {
  H <- as.numeric(quantmod::Hi(x)); L <- as.numeric(quantmod::Lo(x))
  C <- as.numeric(quantmod::Cl(x)); O <- as.numeric(quantmod::Op(x))
  V <- as.numeric(quantmod::Vo(x)); A <- as.numeric(quantmod::Ad(x))
  r <- diff(log(A))                     # adjusted close-to-close log return
  n <- length(C)
  # RAW closes on purpose: overnight + open-to-close then add up exactly to close-to-close
  rcc  <- log(C[-1] / C[-n])            # raw close-to-close, rows 2..n
  onr  <- log(O[-1] / C[-n])            # overnight
  intr <- log(C[-1] / O[-1])            # open-to-close
  vt   <- stats::var(rcc, na.rm = TRUE)
  data.frame(
    n            = length(H),
    na_any       = sum(!stats::complete.cases(cbind(O, H, L, C, V, A))),
    na_HL        = sum(is.na(H) | is.na(L)),
    zero_range   = sum(H == L, na.rm = TRUE),          # -> RR = 0 -> log(0)
    H_lt_L       = sum(H < L, na.rm = TRUE),           # should be 0
    C_outside_HL = sum(C > H | C < L, na.rm = TRUE),   # after-hours prints
    O_outside_HL = sum(O > H | O < L, na.rm = TRUE),
    zero_volume  = sum(V == 0, na.rm = TRUE),
    nonpos_price = sum(c(O, H, L, C, A) <= 0, na.rm = TRUE),
    max_logHL    = round(max(log(H / L), na.rm = TRUE), 4),
    max_absret   = round(max(abs(r), na.rm = TRUE), 4),
    # mean log(H/L) / mean |return|; _cc = close-to-close, _oc = open-to-close (GBM benchmark: 2)
    range_ratio_cc = round(mean(log(H / L), na.rm = TRUE) /
                           mean(abs(r), na.rm = TRUE), 3),
    range_ratio_oc = round(mean(log(H / L), na.rm = TRUE) /
                           mean(abs(log(C / O)), na.rm = TRUE), 3),
    # shares of var(raw close-to-close): intraday + overnight + 2cov = 1
    intraday_share  = round(stats::var(intr, na.rm = TRUE) / vt, 3),
    overnight_share = round(stats::var(onr,  na.rm = TRUE) / vt, 3),
    cov_share       = round(2 * stats::cov(onr, intr, use = "complete.obs") / vt, 3),
    # 1 / intraday share = expected GARCH / RR variance ratio
    pred_garch_over_RR = round(1 / (stats::var(intr, na.rm = TRUE) / vt), 3),
    # mean Parkinson RR / mean squared open-to-close return
    RR_over_intraday_var = round(
      mean((log(H / L))^2 / (4 * log(2)), na.rm = TRUE) /
      mean(log(C / O)^2, na.rm = TRUE), 3)
  )
}

# stopifnot() on screen_quality(): zero_range, H_lt_L, na_HL, nonpos_price all 0; shares sum to 1
# -> invisible TRUE, or stops the run
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

# list of k series stacked in k panels via plot.zoo (plot.xts ignores mfrow); `main` = outer title
# -> draws only, returns NULL
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

# par(mar) for grids of acf()/pacf() panels; top 3.8 keeps each title off the panel above
ACF_MAR <- c(4, 4, 3.8, 1)

# acf() plot without lag 0 + ARCH-robust 95% band in red (dashed blue = iid band)
# -> draws; returns robust_acf_se() data.frame(lag, acf, se_iid, se_robust, inflation) invisibly
acf_robust <- function(x, lag.max = 40, main = NULL, ...) {
  ra <- robust_acf_se(x, K = lag.max)
  a  <- stats::acf(as.numeric(x), lag.max = lag.max, plot = FALSE)
  a$acf <- a$acf[-1, , , drop = FALSE]; a$lag <- a$lag[-1, , , drop = FALSE]
  plot(a, main = main, ylim = range(a$acf, 2 * ra$se_robust, -2 * ra$se_robust), ...)
  graphics::lines(ra$lag,  1.96 * ra$se_robust, col = "red", lwd = 1.6)
  graphics::lines(ra$lag, -1.96 * ra$se_robust, col = "red", lwd = 1.6)
  invisible(ra)
}
# same for pacf(); the red band reuses the ACF's robust se -> same invisible data.frame
pacf_robust <- function(x, lag.max = 40, main = NULL, ...) {
  ra <- robust_acf_se(x, K = lag.max)
  p  <- stats::pacf(as.numeric(x), lag.max = lag.max, plot = FALSE)
  plot(p, main = main, ylim = range(p$acf, 2 * ra$se_robust, -2 * ra$se_robust), ...)
  graphics::lines(ra$lag,  1.96 * ra$se_robust, col = "red", lwd = 1.6)
  graphics::lines(ra$lag, -1.96 * ra$se_robust, col = "red", lwd = 1.6)
  invisible(ra)
}

# grey full-height rect() per episode (named list of c(start, end)) on the open plot -> draws only
# call AFTER plotting, then redraw lines() + box() (rect covers them); POSIXct x: xfun = as.POSIXct
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

# force a Date index (POSIXct -> as.Date) -> same xts, Date-indexed
as_date_xts <- function(x) {
  idx <- zoo::index(x)
  if (inherits(idx, "Date")) return(x)
  xts::xts(zoo::coredata(x), order.by = as.Date(idx))
}

# n largest (or smallest) values of a 1-column xts -> data.frame(rank, date, <value_name>)
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

# sample ACF + ARCH-robust se of rho_k, k = 1..K (Diebold 1986)
# -> data.frame(lag, acf, se_iid, se_robust, inflation); inflation = se_robust / se_iid
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

# Ljung-Box Q(lags): textbook n(n+2) form + ARCH-robust sum((acf / se_robust)^2), both chi2(lags)
# -> 1-row data.frame(lags, LB_standard, p_standard, LB_robust, p_robust, mean_se_inflation)
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

# R^2 = 1 - SSR / SST from sums of squares (not 1 - var(resid)/var(y)) -> one number
pseudo_r2 <- function(resid, y) {
  resid <- as.numeric(resid); y <- as.numeric(y)
  1 - sum(resid^2) / sum((y - mean(y))^2)
}


## ---- 5. Realized measures --------------------------------------------------

# Parkinson (1980) daily range variance (log H - log L)^2 / (4 log 2), RAW H/L (adjustment cancels)
# -> 1-column xts named `name`; H == L days (RR = 0) set to NA
realized_range <- function(x, name = "RR") {
  rr <- (log(quantmod::Hi(x)) - log(quantmod::Lo(x)))^2 / (4 * log(2))
  bad <- !is.finite(rr) | rr <= 0
  rr[bad] <- NA
  colnames(rr) <- name
  rr
}

# monthly RV = sum of squared daily returns (no demeaning), plus French-Schwert-Stambaugh RVm_fss
# months with < min_obs days dropped; dated at the month's last (stamp = "end") or first trading day
# -> monthly xts with columns RVm | RVm_fss | nobs
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

# HAR design (Corsi 2009): y = v[t+1], regressors = trailing means of v over each horizon in p
# filter(sides = 1) = past only; messages if the index has gaps > 4 days
# -> data.frame(y, h1, h5, h22) (one h<k> per p), NA rows dropped: first max(p)-1 and the last
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

# dates of the HAR fitted values (target days: index from position max(p) + 1) -> index vector
har_fitted_dates <- function(v, p = c(1, 5, 22)) {
  idx <- if (inherits(v, c("xts", "zoo"))) zoo::index(v) else seq_along(v)
  m   <- max(as.integer(p))
  idx[(m + 1L):length(idx)]
}

# HAR by OLS (lm on har_design) + Newey-West HAC se (prewhite = FALSE = plain NW 1987)
# jensen = TRUE: fitted_lvl = exp(fit + sigma2/2) (mean) instead of exp(fit) (median)
# -> list(fit, design, hac, vcov_hac, r2, persistence, sigma2, fitted_log/_lvl (xts), p)
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

# QLIKE loss z - log(z) - 1, z = proxy / f (Patton 2011) -> vector; NA unless both finite and > 0
qlike <- function(proxy, f) {
  ok <- is.finite(proxy) & is.finite(f) & proxy > 0 & f > 0
  out <- rep(NA_real_, length(proxy))
  z   <- proxy[ok] / f[ok]
  out[ok] <- z - log(z) - 1
  out
}

# squared-error loss -> vector
mse_loss <- function(proxy, f) (proxy - f)^2

# Mincer-Zarnowitz: lm(proxy ~ f) with Newey-West se + Wald test of (a, b) = (0, 1)
# -> list(fit, hac, r2, a, b, wald_stat, wald_df, wald_p, n)
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

# Wald test of coef(m) == value with vcov V, chi2 with df = number of coefs -> list(stat, df, p)
car_free_wald <- function(m, V, value) {
  b <- stats::coef(m)
  stopifnot(length(b) == length(value))
  d <- as.numeric(b - value)
  W <- as.numeric(t(d) %*% solve(V) %*% d)
  k <- length(b)
  list(stat = W, df = k, p = stats::pchisq(W, df = k, lower.tail = FALSE))
}

# forecast encompassing: lm(proxy ~ f1 + f2) with Newey-West se; names2 labels the slopes in `coef`
# -> list(fit, hac, r2, coef, p_drop_f2, p_drop_f1, corr_f1_f2, vif, combined, n); combined = fitted
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

# Diebold-Mariano on the loss differential l1 - l2, HLN small-sample correction, t(n-1) p-value
# hand-rolled: forecast::dm.test takes errors, not QLIKE losses -> list(stat, n, mean_diff, p)
dm_test <- function(l1, l2, h = 1, hln = TRUE) {
  d  <- l1 - l2
  d  <- d[is.finite(d)]
  n  <- length(d)
  dbar <- mean(d)
  g0 <- stats::var(d) * (n - 1) / n          # long-run var: Bartlett, h-1 lags (h = 1: plain var)
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

# largest VAR lag p <= lag_max with Tn - p rows >= min_obs_per_param * (N*p + 1 coefs per equation)
# -> integer p, 0 if no lag is estimable
var_max_lag <- function(Tn, N, lag_max = 8L, min_obs_per_param = 10) {
  feasible <- vapply(seq_len(lag_max), function(k) {
    (N * k + 1) * min_obs_per_param <= (Tn - k)
  }, logical(1))
  if (!any(feasible)) return(0L)
  max(which(feasible))
}

# vars::VARselect up to the var_max_lag() cap, lag picked by `ic` ("SC(n)" = BIC)
# -> list(p, cap, ic_choice, feasible, selection); p = 0, feasible = FALSE when nothing is estimable
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

# significance stars: *** p < .01, ** < .05, * < .10 -> character vector
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .01, "***",
                     ifelse(p < .05, "**", ifelse(p < .10, "*", ""))))

# urca test object -> 1-row data.frame(test, statistic, stat_name, cv_1pct, cv_5pct, cv_10pct)
# `which` picks the stat from @teststat (e.g. tau, phi1, phi2, phi3)
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

# crisis windows from published event dates
# -> named list dotcom | gfc | euro | covid | banks23, each c(start, end) as "YYYY-MM-DD"
case3_episodes <- function() {
  list(
    dotcom  = c("2000-03-24", "2002-10-09"),  # S&P 500 peak to trough
    gfc     = c("2007-08-01", "2009-06-30"),  # quant quake / BNP freeze
    euro    = c("2011-07-01", "2012-07-26"),  # euro stress
    covid   = c("2020-02-19", "2020-06-30"),  # pre-covid S&P peak
    banks23 = c("2023-02-15", "2023-06-30")   # SVB/Signature
  )
}

# GFC alone as the contiguous (seam-free) crisis benchmark -> list(gfc = c(start, end))
case3_clean_crisis <- function() list(gfc = c("2007-08-01", "2009-06-30"))
