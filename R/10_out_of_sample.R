###############################################################################
# Case 3 -- part 10: GARCH vs HAR one-step-ahead forecast evaluation
#
# Copy this file to R/parts/10_out_of_sample.R in the actual case3 project.
# It uses the shared 00_setup.R / part toolkit and saves outputs in results/.
#
# Both models forecast the same target: squared daily log returns. The GARCH
# forecast is close-to-close conditional variance; the HAR forecast uses the
# previous 1, 5, and 22 trading days of log squared returns. The final 20% of
# dates are held out chronologically. Parameters are estimated only on the
# training sample; forecasts then update one day at a time using observed data.
###############################################################################

PARTS_DIR <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  for (e in sys.frames()) if (is.character(e$ofile)) f <- c(e$ofile, f)
  prev <- if (exists("PARTS_DIR")) get("PARTS_DIR") else character(0)
  cand <- c(dirname(normalizePath(f, mustWork = FALSE)), prev, "R/parts", "parts", ".")
  hit  <- cand[file.exists(file.path(cand, "00_setup.R"))]
  if (!length(hit)) stop("can't find R/parts/00_setup.R -- put this file in R/parts/")
  normalizePath(hit[1], winslash = "/")
})
source(file.path(PARTS_DIR, "00_setup.R"))

need("03_exploration", c("r"))
part_begin("10_out_of_sample")

say("10. OUT-OF-SAMPLE GARCH VS HAR FORECAST EVALUATION")

TRAIN_SHARE <- 0.80
HAR_LAGS <- c(1L, 5L, 22L)
GARCH_ORDER_GRID <- list(c(1L, 1L), c(1L, 2L), c(2L, 1L), c(2L, 2L))

ret <- as.numeric(r)
ret_dates <- zoo::index(r)
stopifnot(all(is.finite(ret)), length(ret) > 100L)

n_all <- length(ret)
n_test <- max(1L, floor(n_all * (1 - TRAIN_SHARE)))
train_n <- n_all - n_test
test_idx <- seq.int(train_n + 1L, n_all)
train_ret <- ret[seq_len(train_n)]
proxy <- ret^2

if (train_n <= max(HAR_LAGS) + 30L)
  stop("The training sample is too short for the 22-day HAR model.")

## ---- Choose the GARCH order and family using training data only --------------

fit_garch_order <- function(order, x) {
  spec <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = order),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std")

  tryCatch({
    fit <- rugarch::ugarchfit(spec = spec, data = x, solver = "hybrid")
    if (fit@fit$convergence != 0L) return(NULL)
    bic <- as.numeric(rugarch::infocriteria(fit)["Bayes", 1L])
    list(order = order, spec = spec, fit = fit, BIC = bic)
  }, error = function(e) {
    message("  GARCH(", order[1], ",", order[2], ") failed: ", conditionMessage(e))
    NULL
  })
}

garch_candidates <- Filter(Negate(is.null), lapply(
  GARCH_ORDER_GRID, fit_garch_order, x = train_ret))
if (!length(garch_candidates)) stop("No candidate GARCH model converged on the training sample.")

garch_order_selection <- do.call(rbind, lapply(garch_candidates, function(z) {
  data.frame(p = z$order[1], q = z$order[2], BIC = z$BIC)
}))
garch_selected_order <- garch_candidates[[which.min(vapply(
  garch_candidates, `[[`, numeric(1), "BIC"))]]

tab("10_oos_garch_order_selection", garch_order_selection, digits = 5,
    caption = "GARCH order selection using training observations only (Student-t errors)")

# Match the specifications compared in part 07, but repeat selection using
# training data so no holdout information influences the chosen GARCH model.
selected_order <- garch_selected_order$order
garch_model_specs <- list(
  sGARCH_norm = rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = selected_order),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "norm"),
  sGARCH_std = rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = selected_order),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"),
  eGARCH_std = rugarch::ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = selected_order),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"),
  gjrGARCH_std = rugarch::ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = selected_order),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
    distribution.model = "std"))

fit_garch_model <- function(label) {
  tryCatch({
    spec <- garch_model_specs[[label]]
    fit <- rugarch::ugarchfit(spec = spec, data = train_ret, solver = "hybrid")
    if (fit@fit$convergence != 0L) return(NULL)
    bic <- as.numeric(rugarch::infocriteria(fit)["Bayes", 1L])
    list(model = label, spec = spec, fit = fit, BIC = bic)
  }, error = function(e) {
    message("  ", label, " failed: ", conditionMessage(e))
    NULL
  })
}

garch_model_candidates <- Filter(Negate(is.null), lapply(
  names(garch_model_specs), fit_garch_model))
if (!length(garch_model_candidates))
  stop("No candidate GARCH specification converged on the training sample.")
garch_model_selection <- do.call(rbind, lapply(garch_model_candidates, function(z) {
  data.frame(model = z$model, BIC = z$BIC)
}))
garch_selected_model <- garch_model_candidates[[which.min(vapply(
  garch_model_candidates, `[[`, numeric(1), "BIC"))]]
tab("10_oos_garch_model_selection", garch_model_selection, digits = 5,
    caption = "GARCH family selection using training observations only")
say("Training-sample BIC selected ", garch_selected_model$model,
    " at GARCH(", selected_order[1], ",", selected_order[2], ")")

## ---- Rolling one-day GARCH forecasts -----------------------------------------

# Reserve the test tail in ugarchfit. ugarchforecast then rolls the variance
# state forward using each newly observed return while keeping fitted parameters
# fixed, which produces genuinely one-step-ahead forecasts without refitting
# thousands of times.
garch_fit_oos <- rugarch::ugarchfit(
  spec = garch_selected_model$spec, data = ret, out.sample = n_test,
  solver = "hybrid")
if (garch_fit_oos@fit$convergence != 0L)
  stop("The selected GARCH model did not converge on the training sample.")

garch_fc <- rugarch::ugarchforecast(
  garch_fit_oos, n.ahead = 1L, n.roll = n_test - 1L)
garch_sigma_fc <- as.numeric(sigma(garch_fc)[1L, ])
if (length(garch_sigma_fc) != n_test)
  stop("Unexpected number of rolling GARCH forecasts.")
garch_var_fc <- garch_sigma_fc^2

## ---- HAR forecasts of the same variance proxy --------------------------------

log_proxy <- rep(NA_real_, n_all)
positive_proxy <- is.finite(proxy) & proxy > 0
log_proxy[positive_proxy] <- log(proxy[positive_proxy])

har_training_rows <- lapply(seq.int(max(HAR_LAGS) + 1L, train_n), function(target) {
  origin <- target - 1L
  i1 <- origin
  i5 <- seq.int(origin - 4L, origin)
  i22 <- seq.int(origin - 21L, origin)
  needed <- c(target, i1, i5, i22)
  if (!all(is.finite(log_proxy[needed]))) return(NULL)
  data.frame(y = log_proxy[target],
             h1 = log_proxy[i1],
             h5 = mean(log_proxy[i5]),
             h22 = mean(log_proxy[i22]))
})
har_training <- do.call(rbind, Filter(Negate(is.null), har_training_rows))
if (is.null(har_training) || nrow(har_training) < 30L)
  stop("Too few valid training rows for the HAR model.")

har_fit <- stats::lm(y ~ h1 + h5 + h22, data = har_training)
har_sigma2 <- sum(stats::residuals(har_fit)^2) / stats::df.residual(har_fit)

har_var_fc <- rep(NA_real_, n_test)
for (j in seq_len(n_test)) {
  target <- train_n + j
  origin <- target - 1L
  i1 <- origin
  i5 <- seq.int(origin - 4L, origin)
  i22 <- seq.int(origin - 21L, origin)
  if (!all(is.finite(log_proxy[c(i1, i5, i22)]))) next

  newdata <- data.frame(h1 = log_proxy[i1],
                        h5 = mean(log_proxy[i5]),
                        h22 = mean(log_proxy[i22]))
  pred_log <- as.numeric(stats::predict(har_fit, newdata = newdata))[1L]
  # Jensen correction turns the log-scale forecast into a conditional mean,
  # which is the target that QLIKE is designed to score.
  har_var_fc[j] <- exp(pred_log + har_sigma2 / 2)
}

oos_forecasts <- data.frame(
  date = as.Date(ret_dates[test_idx]),
  actual_squared_return = proxy[test_idx],
  garch_variance = garch_var_fc,
  har_variance = har_var_fc,
  row.names = NULL)

## ---- Score the common holdout -----------------------------------------------

score_forecast <- function(y, f, model_name) {
  q <- qlike(y, f)
  ok <- is.finite(q) & is.finite(y) & is.finite(f)
  data.frame(model = model_name, n = sum(ok),
             RMSE = sqrt(mean((y[ok] - f[ok])^2)),
             MAE = mean(abs(y[ok] - f[ok])),
             QLIKE = mean(q[ok]), row.names = NULL)
}

oos_metrics <- rbind(
  score_forecast(oos_forecasts$actual_squared_return,
                 oos_forecasts$garch_variance, "GARCH"),
  score_forecast(oos_forecasts$actual_squared_return,
                 oos_forecasts$har_variance, "HAR (1/5/22)")
)

qlike_garch <- qlike(oos_forecasts$actual_squared_return,
                     oos_forecasts$garch_variance)
qlike_har <- qlike(oos_forecasts$actual_squared_return,
                   oos_forecasts$har_variance)
dm <- dm_test(qlike_garch, qlike_har, h = 1L)
oos_dm <- data.frame(
  comparison = "GARCH minus HAR QLIKE loss",
  mean_loss_difference = dm$mean_diff,
  DM_statistic = dm$stat,
  p_value = dm$p,
  n = dm$n,
  interpretation = if (dm$mean_diff > 0) "HAR has lower average QLIKE" else
    "GARCH has lower or equal average QLIKE",
  row.names = NULL)

tab("10_oos_metrics", oos_metrics, digits = 6,
    caption = "One-step-ahead forecast accuracy on the same chronological holdout")
tab("10_oos_dm_test", oos_dm, digits = 6,
    caption = "Diebold-Mariano comparison of QLIKE losses (positive difference favors HAR)")

fig("10_oos_forecast_comparison", {
  annualised_actual <- sqrt(252 * oos_forecasts$actual_squared_return) * 100
  annualised_garch <- sqrt(252 * oos_forecasts$garch_variance) * 100
  annualised_har <- sqrt(252 * oos_forecasts$har_variance) * 100
  plot(oos_forecasts$date, annualised_actual, type = "l", col = "grey55",
       xlab = "", ylab = "annualised volatility (%)",
       main = "One-step forecasts vs absolute-return volatility proxy")
  lines(oos_forecasts$date, annualised_garch, col = "steelblue", lwd = 1.2)
  lines(oos_forecasts$date, annualised_har, col = "firebrick", lwd = 1.2)
  legend("topright", bty = "n", lwd = c(1, 1.2, 1.2),
         col = c("grey55", "steelblue", "firebrick"),
         legend = c("observed |return| proxy", "GARCH", "HAR (1/5/22)"))
  grid(col = "grey88")
})

part_save(c("garch_model_selection", "garch_order_selection",
            "garch_selected_model", "garch_selected_order", "garch_fit_oos",
            "har_fit", "oos_dm", "oos_forecasts", "oos_metrics"))
