###############################################################################
# EBC4008 Case 3 -- part 05: ARIMA on y, log(y), r
# Reads : results/rds/03_exploration.rds (y, ly, rn)
# Writes: tables 05_correlogram_values, 05_arima_coefs, 05_model_comparison,
#         05_arima_equivalence; figures 05_correlograms, 05_residuals_<key>;
#         results/rds/05_arima.rds
###############################################################################

library(forecast)   # Arima(), auto.arima(), arimaorder(), checkresiduals()

d03 <- readRDS("results/rds/03_exploration.rds")
y <- d03$y; ly <- d03$ly; rn <- d03$rn          # price, log price, log return (numeric)

# ACF lags 1..K + ARCH-robust se (Diebold 1986) -> data.frame; inflation = se_robust / se_iid
robust_acf_se <- function(x, K) {
  x <- x - mean(x)
  n <- length(x)
  d <- mean(x^2)^2
  se <- sapply(1:K, function(k) sqrt(mean(x[(k + 1):n]^2 * x[1:(n - k)]^2) / d / n))
  data.frame(lag = 1:K, acf = as.numeric(acf(x, lag.max = K, plot = FALSE)$acf)[-1],
             se_iid = 1 / sqrt(n), se_robust = se, inflation = se * sqrt(n))
}

# acf() without lag 0 (or pacf()) of x, lags 1-40, + ARCH-robust 95% band (red; blue = iid)
plot_robust <- function(x, partial, main) {
  se <- robust_acf_se(x, 40)$se_robust
  if (partial) {
    a <- pacf(x, lag.max = 40, plot = FALSE)
  } else {
    a <- acf(x, lag.max = 40, plot = FALSE)
    a$acf <- a$acf[-1, , , drop = FALSE]; a$lag <- a$lag[-1, , , drop = FALSE]
  }
  plot(a, main = main, ylim = range(a$acf, 2 * se, -2 * se))
  lines(1:40,  1.96 * se, col = "red", lwd = 1.6)
  lines(1:40, -1.96 * se, col = "red", lwd = 1.6)
}

# R^2 = 1 - SSR/SST (not 1 - var(resid)/var(x))
r2 <- function(e, x) 1 - sum(e^2) / sum((x - mean(x))^2)

## --- Route (i): by hand
# acf() levels; ARCH-robust acf/pacf on diff(y), r -> 05_correlograms
png("results/figures/05_correlograms.png", width = 1500, height = 1500, res = 130)
par(mfrow = c(3, 2), mar = c(4, 4, 3.8, 1))   # room for acf() titles
acf(y,  lag.max = 40, main = "ACF of y -- decays glacially (unit root)")
acf(ly, lag.max = 40, main = "ACF of log(y) -- same")
plot_robust(diff(y), FALSE, "ACF of diff(y) (red = ARCH-robust band)")
plot_robust(diff(y), TRUE,  "PACF of diff(y) (red = ARCH-robust band)")
plot_robust(rn, FALSE, "ACF of r = diff(log(y)) (red = ARCH-robust band)")
plot_robust(rn, TRUE,  "PACF of r (red = ARCH-robust band)")
dev.off()

# ACF and PACF of diff(y) and r, lags 1-10; band95 = ARCH-robust 1.96 se -> 05_correlogram_values
cg <- function(x, lab) {
  ra <- robust_acf_se(x, 10)
  setNames(data.frame(ra$acf, as.numeric(pacf(x, lag.max = 10, plot = FALSE)$acf),
                      1.96 * ra$se_robust),
           paste0(c("acf_", "pacf_", "band95_"), lab))
}
cg_tab <- data.frame(lag = 1:10, cg(diff(y), "dy"), cg(rn, "r"))
write.csv(cg_tab, "results/tables/05_correlogram_values.csv")
print(cg_tab, digits = 4)

# Arima(), exact ML; drift if d = 1 (include.mean ignored), mean if d = 0 -> fit_hand_*
HAND <- list(y = c(0, 1, 0), "log(y)" = c(1, 1, 1), r = c(1, 0, 1))   # hand-picked (p,d,q)
fit_hand_y  <- Arima(y,  order = HAND$y,          include.drift = TRUE, method = "ML")
fit_hand_ly <- Arima(ly, order = HAND[["log(y)"]], include.drift = TRUE, method = "ML")
fit_hand_r  <- Arima(rn, order = HAND$r,          include.mean  = TRUE, method = "ML")

## --- Route (ii): auto.arima (Hyndman-Khandakar)
# auto.arima(), d fixed, full grid p, q <= 5, exact ML; AICc y, log(y); BIC r -> fit_auto_*
# max.order = 10: default 5 would cap p + q
fit_auto_y  <- auto.arima(y, d = 1, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
                          ic = "aicc", stepwise = FALSE, approximation = FALSE, trace = FALSE)
fit_auto_ly <- auto.arima(ly, d = 1, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
                          ic = "aicc", stepwise = FALSE, approximation = FALSE, trace = FALSE)
fit_auto_r  <- auto.arima(rn, d = 0, max.p = 5, max.q = 5, max.order = 10, seasonal = FALSE,
                          ic = "bic", stepwise = FALSE, approximation = FALSE, trace = FALSE)

arima_fits <- list(y_hand    = fit_hand_y,  y_auto    = fit_auto_y,
                   logy_hand = fit_hand_ly, logy_auto = fit_auto_ly,
                   r_hand    = fit_hand_r,  r_auto    = fit_auto_r)
# series per fit; route = key after "_"; mlab() -> label "ARIMA(p,d,q) + drift / mean"
SERIES <- c(y_hand = "y", y_auto = "y", logy_hand = "log(y)", logy_auto = "log(y)",
            r_hand = "r", r_auto = "r")
mlab <- function(m) paste0("ARIMA(", paste(arimaorder(m), collapse = ","), ")",
                           if ("drift" %in% names(coef(m))) " + drift" else
                           if ("intercept" %in% names(coef(m))) " + mean" else "")

## --- Estimated models
# coef(), ML (iid) se from var.coef, z, normal p, stars -> arima_coefs, 05_arima_coefs
arima_coefs <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; b <- coef(m); se <- sqrt(diag(m$var.coef))
  data.frame(series = SERIES[[k]], route = sub("^.*_", "", k),
             model = mlab(m), term = names(b), estimate = unname(b), se = unname(se),
             z = unname(b / se), p = unname(2 * pnorm(-abs(b / se))), row.names = NULL)
}))
arima_coefs$sig <- with(arima_coefs, ifelse(is.na(p), "", ifelse(p < .01, "***",
                        ifelse(p < .05, "**", ifelse(p < .10, "*", "")))))
write.csv(arima_coefs, "results/tables/05_arima_coefs.csv")
print(arima_coefs, digits = 5)

# checkresiduals(): plot -> 05_residuals_<key>; Ljung-Box, lags = max(10, p + q + 3),
# df = lags - p - q -> resid_lb
resid_lb <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]
  png(paste0("results/figures/05_residuals_", k, ".png"), width = 1500, height = 950, res = 130)
  checkresiduals(m, test = FALSE)
  dev.off()
  capture.output(lbt <- checkresiduals(m, plot = FALSE, test = "LB"))   # silence its print
  data.frame(key = k, LB = unname(lbt$statistic), df = unname(lbt$parameter),
             p = unname(lbt$p.value), row.names = NULL)
}))

## --- Compare six fits
# info criteria; Ljung-Box textbook + ARCH-robust (same lags/df); R^2 on level and change
# -> arima_cmp, 05_model_comparison
arima_cmp <- do.call(rbind, lapply(names(arima_fits), function(k) {
  m <- arima_fits[[k]]; s <- SERIES[[k]]; d <- arimaorder(m)[2]
  x <- switch(s, y = y, "log(y)" = ly, r = rn)
  e <- as.numeric(residuals(m)); if (d == 1) e <- e[-1]      # d = 1: e[1] is a start-up value
  lev <- if (d == 1) x[-1] else x
  chg <- if (d == 1) diff(x) else x
  lb  <- resid_lb[resid_lb$key == k, ]
  ra  <- robust_acf_se(e, lb$df + sum(arimaorder(m)[c(1, 3)]))   # lags checkresiduals used
  data.frame(series = s, route = sub("^.*_", "", k), model = mlab(m),
             n = m$nobs, k = length(coef(m)), loglik = m$loglik, aicc = m$aicc, bic = m$bic,
             sigma2 = m$sigma2, LB_p = lb$p, LB_df = lb$df,
             LB_robust_p = pchisq(sum((ra$acf / ra$se_robust)^2), lb$df, lower.tail = FALSE),
             R2_level = r2(e, lev), R2_change = r2(e, chg), row.names = NULL)
}))
write.csv(arima_cmp, "results/tables/05_model_comparison.csv")
print(arima_cmp, digits = 6)

# hand log(y) (d = 1) IS hand r (d = 0): loglik, sigma2, n, coefs (drift row = mean row)
# -> equiv, 05_arima_equivalence
equiv <- data.frame(
  quantity = c("log-likelihood", "sigma^2", "n used", names(coef(fit_hand_ly))),
  on_log_y_d1 = c(fit_hand_ly$loglik, fit_hand_ly$sigma2, fit_hand_ly$nobs,
                  unname(coef(fit_hand_ly))),
  on_r_d0     = c(fit_hand_r$loglik, fit_hand_r$sigma2, fit_hand_r$nobs,
                  unname(coef(fit_hand_r))))
equiv$difference <- equiv$on_log_y_d1 - equiv$on_r_d0
write.csv(equiv, "results/tables/05_arima_equivalence.csv")
print(equiv, digits = 8)

saveRDS(mget(c("arima_cmp", "arima_coefs", "arima_fits", "equiv", "resid_lb")),
        "results/rds/05_arima.rds")
