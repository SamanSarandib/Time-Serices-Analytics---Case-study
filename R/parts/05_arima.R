###############################################################################
# EBC4008 Case 3 -- part 05: ARIMA for the price, the log price and the returns
# Needs : ly, rn, y <- 03_exploration
# Makes : arima_cmp, arima_coefs, arima_fits, equiv, resid_lb
# Writes: tables 05_correlogram_values, 05_arima_coefs, 05_model_comparison,
#         05_arima_equivalence; figures 05_correlograms, 05_residuals_<series>_<route>
###############################################################################

# find R/parts/ and load the shared setup (Rscript, source() or RStudio Ctrl+Enter)
PARTS_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  for (e in sys.frames()) if (is.character(e$ofile)) f <- c(e$ofile, f)
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

## --- Route (i): by hand, from the correlogram
# acf() on the levels; acf_robust() / pacf_robust() (blue iid + red ARCH-robust band) on diff(y), r
# -> 05_correlograms: 3x2 panels, rows = ACF y, log(y) | ACF, PACF diff(y) | ACF, PACF r
fig("05_correlograms", {
  par(mfrow = c(3, 2), mar = ACF_MAR)   # ACF_MAR: room for acf()'s titles
  acf(y,  lag.max = 40, main = "ACF of y -- decays glacially (unit root)")
  acf(ly, lag.max = 40, main = "ACF of log(y) -- same")
  acf_robust(diff(y),  lag.max = 40, main = "ACF of diff(y) (red = ARCH-robust band)")
  pacf_robust(diff(y), lag.max = 40, main = "PACF of diff(y) (red = ARCH-robust band)")
  acf_robust(rn,  lag.max = 40, main = "ACF of r = diff(log(y)) (red = ARCH-robust band)")
  pacf_robust(rn, lag.max = 40, main = "PACF of r (red = ARCH-robust band)")
}, height = 1500)

# same correlograms as numbers: robust_acf_se() + pacf(), lags 1-10, band95 = 1.96 * robust se
# -> 05_correlogram_values: lag | acf_dy | pacf_dy | band95_dy | acf_r | pacf_r | band95_r
cg <- function(x, lab) {
  ra <- robust_acf_se(x, K = 10)
  setNames(data.frame(ra$acf, as.numeric(pacf(x, lag.max = 10, plot = FALSE)$acf),
                      1.96 * ra$se_robust),
           paste0(c("acf_", "pacf_", "band95_"), lab))
}
cg_tab <- data.frame(lag = 1:10, cg(diff(y), "dy"), cg(rn, "r"))
tab("05_correlogram_values", cg_tab, digits = 4,
    caption = "ACF and PACF of diff(y) and r, first 10 lags (band95 = ARCH-robust 1.96 se)")

# hand (p, d, q), read off the correlograms
HAND <- list(y = c(0, 1, 0), "log(y)" = c(1, 1, 1), r = c(1, 0, 1))
# forecast::Arima(), exact ML; drift for d = 1, mean for d = 0 (include.mean is ignored when d = 1)
# -> 3 Arima objects
fit_hand_y  <- Arima(y,  order = HAND$y,          include.drift = TRUE, method = "ML")
fit_hand_ly <- Arima(ly, order = HAND[["log(y)"]], include.drift = TRUE, method = "ML")
fit_hand_r  <- Arima(rn, order = HAND$r,          include.mean  = TRUE, method = "ML")

## --- Route (ii): auto.arima (Hyndman-Khandakar)
# forecast::auto.arima(), d fixed from 04, full grid p, q <= 5 (stepwise = FALSE), exact ML;
# max.order = 10 so every (p, q) is reachable; AICc for y, log(y), BIC for r
# -> 3 Arima objects
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

# -> arima_fits: named list of 6 Arima objects (y_hand, y_auto, logy_hand, logy_auto, r_hand, r_auto)
arima_fits <- list(y_hand    = fit_hand_y,  y_auto    = fit_auto_y,
                   logy_hand = fit_hand_ly, logy_auto = fit_auto_ly,
                   r_hand    = fit_hand_r,  r_auto    = fit_auto_r)
# series and route per fit; mlab() -> label like "ARIMA(1,1,1) + drift" / "ARIMA(1,0,1) + mean"
SERIES <- c(y_hand = "y", y_auto = "y", logy_hand = "log(y)", logy_auto = "log(y)",
            r_hand = "r", r_auto = "r")
ROUTE  <- sub("^.*_", "", names(arima_fits))
mlab <- function(m) paste0("ARIMA(", paste(arimaorder(m), collapse = ","), ")",
                           if ("drift" %in% names(coef(m))) " + drift" else
                           if ("intercept" %in% names(coef(m))) " + mean" else "")

## --- Present the estimated models
# coef() with ML (iid) se from var.coef, z and normal p (summary.Arima prints no p)
# -> 05_arima_coefs: series | route | model | term | estimate | se | z | p | sig
arima_coefs <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; b <- coef(m); se <- sqrt(diag(m$var.coef))
  data.frame(series = SERIES[[k]], route = ROUTE[match(k, names(arima_fits))],
             model = mlab(m), term = names(b), estimate = unname(b), se = unname(se),
             z = unname(b / se), p = unname(2 * pnorm(-abs(b / se))), row.names = NULL)
}))
arima_coefs$sig <- stars(arima_coefs$p)
tab("05_arima_coefs", arima_coefs, digits = 5,
    caption = "Estimated ARIMA models: coefficients, ML (iid) standard errors, z and p")

# forecast::checkresiduals(): plot + Ljung-Box on 10 lags (p + q + 3 if more), df = lags - p - q
# -> figure 05_residuals_<key> (residuals over time, ACF, histogram); resid_lb: key | LB | df | p
resid_lb <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]
  fig(paste0("05_residuals_", k), checkresiduals(m, test = FALSE))
  invisible(capture.output(lbt <- checkresiduals(m, plot = FALSE, test = "LB")))
  data.frame(key = k, LB = unname(lbt$statistic), df = unname(lbt$parameter),
             p = unname(lbt$p.value), row.names = NULL)
}))

## --- How does the fit on returns compare with the fit on prices?
# pseudo_r2() on the level and on the daily change (diff(y) for y, r for log(y) and r);
# LB_p = checkresiduals' Ljung-Box, LB_robust_p = lb_robust() (ARCH-robust), same lags and df
# -> 05_model_comparison: series | route | model | n | k | loglik | aicc | bic | sigma2 |
#    LB_p | LB_df | LB_robust_p | R2_level | R2_change
arima_cmp <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; s <- SERIES[[k]]; d <- arimaorder(m)[2]
  x <- switch(s, y = y, "log(y)" = ly, r = rn)
  e <- as.numeric(residuals(m)); if (d == 1) e <- e[-1]      # d = 1: residual 1 is a start-up value
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

# hand log(y) model (d = 1) vs hand r model (d = 0): loglik, sigma2, n and coefs side by side
# -> 05_arima_equivalence: quantity | on_log_y_d1 | on_r_d0 | difference
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
