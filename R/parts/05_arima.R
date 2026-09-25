###############################################################################
# EBC4008 Case 3 -- part 05: ARIMA for the price, the log price and the returns
#
# Inputs come from need(),
# outputs go to results/parts/05_arima.rds via part_save().
#
# Needs : ly, rn, y <- 03_exploration
# Makes : arima_cmp, arima_coefs, arima_fits, equiv, resid_lb
# Writes: tables 05_correlogram_values, 05_arima_coefs, 05_model_comparison,
#         05_arima_equivalence; figures 05_correlograms, 05_residuals_<series>_<route>
#
# What the case asks, for y (price), log(y) and r (log return):
#   (i)  specify an ARIMA by hand from the correlogram
#   (ii) let auto.arima search the class (Hyndman-Khandakar 2008)
#   present the estimated models; then compare the fit on the returns with the
#   fit on the prices, and explain the difference.
# Returns are built in 03; their unit-root tests are in 04 (r is I(0), y and
# log(y) are I(1): so d = 1 for the levels, d = 0 for r).
# Package: forecast. Arima() is stats::arima with a drift option and the
# information criteria attached; auto.arima() is the Hyndman-Khandakar search;
# checkresiduals() plots the residuals and runs Ljung-Box with the right df.
###############################################################################

# bootstrap: find R/parts/
PARTS_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  for (e in sys.frames()) if (is.character(e$ofile)) f <- c(e$ofile, f)
  # Ctrl+Enter in RStudio sources no file: ask RStudio which one is open in the editor
  if (!length(f) && exists(".rs.api.getSourceEditorContext"))
    f <- as.character(.rs.api.getSourceEditorContext()$path)
  prev <- if (exists("PARTS_DIR")) get("PARTS_DIR") else character(0)
  cand <- c(dirname(normalizePath(f, mustWork = FALSE)), prev, "R/parts", "parts", ".")
  hit  <- cand[file.exists(file.path(cand, "00_setup.R"))]
  if (!length(hit)) stop("can't find R/parts/00_setup.R -- setwd() to the case3 folder")
  normalizePath(hit[1], winslash = "/")
})
source(file.path(PARTS_DIR, "00_setup.R"))

need("03_exploration", c("ly", "rn", "y"))

part_begin("05_arima")

## ---- 5. ARIMA for y, log(y) and r -----------------------------------------------

say("5. ARIMA FOR THE PRICE, THE LOG PRICE AND THE RETURNS")

## --- Route (i): by hand, from the correlogram
# Levels: the ACF of y and of log(y) starts near 1 and decays glacially -- the
# unit root 04 found -- so difference once (d = 1).
# Differences: read against the red ARCH-robust band (the blue iid band flags
# lags that are only volatility clustering, see 03):
#   diff(y): lag 1 = -0.04, inside the robust band (+-0.07); only lags 20 and
#            26 of 40 poke out (the chance rate) -> nothing to model beyond
#            the drift
#   r:       lag 1 = -0.08, outside the robust band (+-0.06), in both the ACF
#            and the PACF, nothing systematic after it -> an MA(1) or an AR(1);
#            ARMA(1,1) nests the two
# (exact values in 05_correlogram_values)
fig("05_correlograms", {
  par(mfrow = c(3, 2), mar = ACF_MAR)   # ACF_MAR: room for acf()'s titles
  acf(y,  lag.max = 40, main = "ACF of y -- decays glacially (unit root)")
  acf(ly, lag.max = 40, main = "ACF of log(y) -- same")
  acf_robust(diff(y),  lag.max = 40, main = "ACF of diff(y) (red = ARCH-robust band)")
  pacf_robust(diff(y), lag.max = 40, main = "PACF of diff(y) (red = ARCH-robust band)")
  acf_robust(rn,  lag.max = 40, main = "ACF of r = diff(log(y)) (red = ARCH-robust band)")
  pacf_robust(rn, lag.max = 40, main = "PACF of r (red = ARCH-robust band)")
}, height = 1500)

# the same correlograms as numbers, first 10 lags, with the robust 95% band
cg <- function(x, lab) {
  ra <- robust_acf_se(x, K = 10)
  setNames(data.frame(ra$acf, as.numeric(pacf(x, lag.max = 10, plot = FALSE)$acf),
                      1.96 * ra$se_robust),
           paste0(c("acf_", "pacf_", "band95_"), lab))
}
cg_tab <- data.frame(lag = 1:10, cg(diff(y), "dy"), cg(rn, "r"))
tab("05_correlogram_values", cg_tab, digits = 4,
    caption = "ACF and PACF of diff(y) and r, first 10 lags (band95 = ARCH-robust 1.96 se)")

# The hand models, (p, d, q) read off the correlograms:
#   y      ARIMA(0,1,0) + drift: a random walk with drift. (ARIMA(1,1,1) on y
#          puts AR and MA near 0 with no standard errors -- nothing to identify.)
#   log(y) ARIMA(1,1,1) + drift. Its AR (0.51) and MA (-0.59) nearly cancel
#          (correlation of the estimates -0.997; net lag-1 effect -0.07):
#          AICc prefers it to the MA(1) by 2.5, BIC prefers the MA(1) by 11.
#   r      ARMA(1,1) + mean: log(y)'s model written on r (differencing is part
#          of the model, see the equivalence table below)
# Drift / mean = the average daily change, i.e. the long-run upward trend.
HAND <- list(y = c(0, 1, 0), "log(y)" = c(1, 1, 1), r = c(1, 0, 1))
# include.drift for d = 1, include.mean for d = 0 (with d = 1 include.mean is
# silently ignored). method = "ML": exact likelihood -- comparable within a
# series and between log(y) and r (same data), NOT with y (dollars, no Jacobian).
fit_hand_y  <- Arima(y,  order = HAND$y,          include.drift = TRUE, method = "ML")
fit_hand_ly <- Arima(ly, order = HAND[["log(y)"]], include.drift = TRUE, method = "ML")
fit_hand_r  <- Arima(rn, order = HAND$r,          include.mean  = TRUE, method = "ML")

## --- Route (ii): auto.arima (Hyndman-Khandakar)
# How it works: d from repeated KPSS tests, then p and q by an information
# criterion. d is fixed at 04's Pantula answer: left free, KPSS picks d = 1 for
# log(y) and 0 for r, but d = 2 for y (KPSS rejects on diff(y), see 04 --
# over-differencing). stepwise = FALSE -> fit the whole grid (p, q <= 5)
# instead of the greedy hill-climb; approximation = FALSE -> exact ML, not the
# CSS shortcut. It drops any model with a root within 1% of the unit circle,
# whatever its likelihood (log(y)'s pick has AR roots of modulus 1.0100 --
# right on that cut-off; better-AICc models with drift were discarded by it).
# max.order = 10 so it can reach every (p, q) in the grid; the default 5 would
# rule out e.g. ARIMA(3,1,5).
# Criterion: AICc (auto.arima's default) for y and log(y); BIC for r, which
# penalises extra terms harder -- at n ~ 6,900 returns are nearly unpredictable
# and AICc would over-fit noise. log(y) with d = 1 and r with d = 0 are the
# same data, so their two auto picks differ only by the criterion: AICc on r
# also picks ARMA(3,5) (loglik 16239.57), BIC on log(y) also picks ARIMA(0,1,1)
# (ma1 -0.082).
fit_auto_y <- cache("arima_auto_y", {
  auto.arima(y, d = 1, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
             ic = "aicc", stepwise = FALSE, approximation = FALSE, trace = FALSE)
})
fit_auto_ly <- cache("arima_auto_logy", {
  auto.arima(ly, d = 1, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
             ic = "aicc", stepwise = FALSE, approximation = FALSE, trace = FALSE)
})
fit_auto_r <- cache("arima_auto_r", {
  auto.arima(rn, d = 0, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
             ic = "bic", stepwise = FALSE, approximation = FALSE, trace = FALSE)
})

arima_fits <- list(y_hand    = fit_hand_y,  y_auto    = fit_auto_y,
                   logy_hand = fit_hand_ly, logy_auto = fit_auto_ly,
                   r_hand    = fit_hand_r,  r_auto    = fit_auto_r)
SERIES <- c(y_hand = "y", y_auto = "y", logy_hand = "log(y)", logy_auto = "log(y)",
            r_hand = "r", r_auto = "r")
ROUTE  <- sub("^.*_", "", names(arima_fits))
mlab <- function(m) paste0("ARIMA(", paste(arimaorder(m), collapse = ","), ")",
                           if ("drift" %in% names(coef(m))) " + drift" else
                           if ("intercept" %in% names(coef(m))) " + mean" else "")

## --- Present the estimated models
# summary.Arima prints no p-values -> compute them.
# These are ML standard errors, which assume iid errors. Under ARCH they are
# too small (03: the robust se of the lag-1 autocorrelation is 2.6x the iid
# one), so read z and p as optimistic: e.g. the MA(1)'s z of -6.7 is about
# -2.6 with robust se.
arima_coefs <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; b <- coef(m); se <- sqrt(diag(m$var.coef))
  data.frame(series = SERIES[[k]], route = ROUTE[match(k, names(arima_fits))],
             model = mlab(m), term = names(b), estimate = unname(b), se = unname(se),
             z = unname(b / se), p = unname(2 * pnorm(-abs(b / se))), row.names = NULL)
}))
arima_coefs$sig <- stars(arima_coefs$p)
tab("05_arima_coefs", arima_coefs, digits = 5,
    caption = "Estimated ARIMA models: coefficients, ML (iid) standard errors, z and p")

# Residual checks. checkresiduals(): residuals over time, their ACF, their
# histogram, and a Ljung-Box test on 10 lags (p + q + 3 if more), df = lags - p - q.
resid_lb <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]
  fig(paste0("05_residuals_", k), checkresiduals(m, test = FALSE))
  invisible(capture.output(lbt <- checkresiduals(m, plot = FALSE, test = "LB")))
  data.frame(key = k, LB = unname(lbt$statistic), df = unname(lbt$parameter),
             p = unname(lbt$p.value), row.names = NULL)
}))

## --- How does the fit on returns compare with the fit on prices?
# R^2 by sums of squares (pseudo_r2), measured two ways:
#   R2_level  = share of the variance of the SERIES MODELLED that is explained
#   R2_change = share of the variance of the DAILY CHANGE that is explained:
#               diff(y) for the y models, r for the log(y) and r models
#               (for r the series modelled is already the daily change)
# With d = 1, residuals()[1] is a start-up value, not a forecast error: drop it.
# LB_p is the textbook Ljung-Box (checkresiduals'), which over-rejects under
# ARCH; LB_robust_p is the ARCH-robust version from 03, same lags and df.
arima_cmp <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; s <- SERIES[[k]]; d <- arimaorder(m)[2]
  x <- switch(s, y = y, "log(y)" = ly, r = rn)
  e <- as.numeric(residuals(m)); if (d == 1) e <- e[-1]
  lev <- if (d == 1) x[-1] else x
  chg <- if (d == 1) diff(x) else x
  lbdf <- resid_lb$df[resid_lb$key == k]
  lags <- lbdf + sum(arimaorder(m)[c(1, 3)])                 # the lags checkresiduals used
  data.frame(series = s, route = ROUTE[match(k, names(arima_fits))], model = mlab(m),
             n = m$nobs, k = length(coef(m)), loglik = m$loglik, aicc = m$aicc, bic = m$bic,
             sigma2 = m$sigma2,
             LB_p = resid_lb$p[resid_lb$key == k], LB_df = lbdf,
             LB_robust_p = pchisq(lb_robust(e, lags)$LB_robust, lbdf, lower.tail = FALSE),
             R2_level = pseudo_r2(e, lev), R2_change = pseudo_r2(e, chg), row.names = NULL)
}))
tab("05_model_comparison", arima_cmp, digits = 6,
    caption = "All six ARIMA models: information criteria, residual Ljung-Box (textbook and ARCH-robust), R^2 on the level and on the daily change")
# Residuals: the textbook Ljung-Box rejects for all six models, the ARCH-robust
# one for none -- what is left is volatility clustering (GARCH's job, 07), not mean
# dependence. The y models' loglik/AICc/BIC are in dollars: compare them only
# with each other.
# How to read it: on the LEVEL the price models explain ~99.9% -- but only
# because today's price is yesterday's price plus a small move (the unit
# root); the random walk (hand model for y) scores the same. On the daily
# CHANGE, which is what a forecaster needs, no model explains much more than
# 1%: returns are close to unpredictable in their mean. The return model's R^2
# is small because it is measured on the hard part, not because it is a worse
# model.
# auto.arima's picks carry AR and MA root pairs that nearly cancel (y: moduli
# 1.033 vs 1.028; log(y): 1.010 vs 1.017). For log(y) the terms are weakly
# identified (se 0.25-0.36); for y they are sharply estimated but add only 0.7%
# R^2 on diff(y). BIC, which charges more per parameter, prefers the hand model
# for log(y); for y it still prefers auto.arima's pick.

# Why the log(y) and r numbers match: the hand model on log(y) IS the hand
# model on r = diff(log(y)). Same data, same likelihood; coefficients equal up
# to the optimiser (they differ by ~0.001 against se ~0.13), and the drift is
# the mean of r. Only the scale R^2 is measured on differs.
equiv <- data.frame(
  quantity = c("log-likelihood", "sigma^2", "n used", names(coef(fit_hand_ly))),
  on_log_y_d1 = c(fit_hand_ly$loglik, fit_hand_ly$sigma2, fit_hand_ly$nobs,
                  unname(coef(fit_hand_ly))),
  on_r_d0     = c(fit_hand_r$loglik, fit_hand_r$sigma2, fit_hand_r$nobs,
                  unname(coef(fit_hand_r))))
equiv$difference <- equiv$on_log_y_d1 - equiv$on_r_d0
tab("05_arima_equivalence", equiv, digits = 8,
    caption = paste(mlab(fit_hand_ly), "on log(y) IS", mlab(fit_hand_r),
                    "on r (drift row = mean row)"))

part_save(c("arima_cmp", "arima_coefs", "arima_fits", "equiv", "resid_lb"))
