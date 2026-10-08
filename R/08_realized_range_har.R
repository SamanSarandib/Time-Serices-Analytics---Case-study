###############################################################################
# EBC4008 Case 3 -- part 08: realized range, HAR-RR
# Reads : results/rds/01_data.rds (jpm)
#         results/rds/02_quality_screen.rds (qual)
#         results/rds/07_garch.rds (sigma_garch)
# Writes: tables 08_rr_distribution, 08_rr_unit_root, 08_har_coefs, 08_hac_decomposition,
#         08_hac_bandwidth_profile, 08_har_horizon_robustness, 08_scale_decomposition;
#         figures 08_rr_series, 08_rr_histograms, 08_har_vs_garch;
#         results/rds/08_realized_range_har.rds
###############################################################################

library(xts)        # xts(), index(), merge(); also attaches zoo
library(sandwich)   # NeweyWest()
library(bootUR)     # boot_adf()

d01 <- readRDS("results/rds/01_data.rds")
jpm <- d01$jpm
d02 <- readRDS("results/rds/02_quality_screen.rds")
qual <- d02$qual
d07 <- readRDS("results/rds/07_garch.rds")
sigma_garch <- d07$sigma_garch

# HAR (Corsi 2009): v[t+1] on trailing means of v over p periods, OLS, NW lag = max(p) -> list
fit_har <- function(v, p) {
  vv <- as.numeric(v)
  X  <- sapply(p, function(k) as.numeric(stats::filter(vv, rep(1 / k, k), sides = 1)))
  colnames(X) <- paste0("h", p)
  d <- na.omit(data.frame(y = c(vv[-1], NA), X))      # drops first max(p) - 1 rows and the last
  m <- lm(y ~ ., data = d)                             # means are positional (span holidays)
  V <- NeweyWest(m, lag = max(p), prewhite = FALSE, adjust = TRUE)
  dts <- index(v)[-(1:max(p))]                         # fitted value dated at its target y
  stopifnot(length(fitted(m)) == length(dts))          # alignment
  list(fit = m, design = d, hac = lmtest::coeftest(m, vcov. = V), vcov_hac = V,
       r2 = summary(m)$r.squared, persistence = sum(coef(m)[-1]),
       sigma2 = sum(residuals(m)^2) / df.residual(m),
       fitted_log = xts(fitted(m), order.by = dts),
       fitted_lvl = xts(exp(fitted(m)), order.by = dts),     # exp of a log fit = a median
       p = as.integer(p))
}

## --- Realized range
# Parkinson (1980) (log H - log L)^2 / (4 log 2) -> RR_raw (H == L days NA), RR (NAs dropped)
RR_raw <- (log(quantmod::Hi(jpm)) - log(quantmod::Lo(jpm)))^2 / (4 * log(2))
RR_raw[!is.finite(RR_raw) | RR_raw <= 0] <- NA
colnames(RR_raw) <- "RR"
RR <- na.omit(RR_raw)
rr <- as.numeric(RR)

# skewness, excess kurtosis of RR vs log(RR) -> 08_rr_distribution
rr_dist <- as.data.frame(t(c(skew_level = moments::skewness(rr),
                             skew_log   = moments::skewness(log(rr)),
                             kurt_level = moments::kurtosis(rr) - 3,
                             kurt_log   = moments::kurtosis(log(rr)) - 3)))
write.csv(rr_dist, "results/tables/08_rr_distribution.csv")
print(rr_dist, digits = 3)

# RR, log(RR) -> 08_rr_series
png("results/figures/08_rr_series.png", width = 1500, height = 900, res = 130)
par(mfrow = c(2, 1), mar = c(2.4, 4.0, 2.0, 1.0))
plot(zoo::as.zoo(RR), xlab = "", ylab = "", col = "grey20",
     main = "Daily realized range (Parkinson), level -- spikes dominate")
grid(col = "grey88")
plot(zoo::as.zoo(log(RR)), xlab = "", ylab = "", col = "grey20",
     main = "log realized range -- much closer to Gaussian, and this is what we model")
grid(col = "grey88")
dev.off()

# histograms RR, log(RR) + fitted normal (red) -> 08_rr_histograms
png("results/figures/08_rr_histograms.png", width = 1500, height = 950, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
hist(rr, breaks = 200, col = "grey65", border = NA, main = "RR (level)", xlab = "RR")
hist(log(rr), breaks = 100, col = "grey65", border = NA,
     main = "log RR", xlab = "log RR", freq = FALSE)
curve(dnorm(x, mean(log(rr)), sd(log(rr))), add = TRUE, col = "red", lwd = 2)
dev.off()

## --- Is log(RR) I(0)?
# bootstrap ADF: intercept, MAIC, AWB p (p = 0: below 1/B) -> ur_rr, 08_rr_unit_root
maxlag <- trunc(12 * (nrow(RR) / 100)^(1 / 4))
set.seed(1)
ur_rr <- boot_adf(log(rr), data_name = "log(RR)", bootstrap = "AWB", B = 1999,
                  deterministics = "intercept", max_lag = maxlag, criterion = "MAIC",
                  criterion_scale = FALSE, show_progress = FALSE)
rr_ur <- data.frame(series = "log(RR)", deterministics = "intercept", criterion = "MAIC",
                    ceiling = maxlag, lags_used = as.integer(ur_rr$details$`selected lags`),
                    tau = unname(ur_rr$statistic), p_boot = unname(ur_rr$p.value))
rr_ur$reject_5pct <- rr_ur$p_boot < 0.05
write.csv(rr_ur, "results/tables/08_rr_unit_root.csv")
print(rr_ur, digits = 3)

## --- HAR-RR
# log RR on 1/5/22-day means, Newey-West 22 lags -> har_rr
har_rr <- fit_har(log(RR), p = c(1, 5, 22))

# coefficients with Newey-West se + stars (*** < .01, ** < .05, * < .10) -> 08_har_coefs
har_coefs <- data.frame(term = rownames(har_rr$hac), estimate = har_rr$hac[, 1],
                        se_hac = har_rr$hac[, 2], t_hac = har_rr$hac[, 3],
                        p_hac = har_rr$hac[, 4], row.names = NULL)
pv <- har_coefs$p_hac
har_coefs$sig <- ifelse(pv < .01, "***", ifelse(pv < .05, "**", ifelse(pv < .10, "*", "")))
write.csv(har_coefs, "results/tables/08_har_coefs.csv")
print(har_coefs, digits = 4)

# bgtest() orders 1, 5, 22 -> har_bg; dwtest() contrast only (lagged y) -> dw
hr <- residuals(har_rr$fit)
har_bg <- do.call(rbind, lapply(c(1, 5, 22), function(k) {
  bg <- lmtest::bgtest(har_rr$fit, order = k)
  data.frame(order = k, LM = unname(bg$statistic), p = bg$p.value)
}))
dw <- lmtest::dwtest(har_rr$fit)

# Ljung-Box(10) on HAC scores x_t e_t -> sc: p-values (Intercept) | h1 | h5 | h22
X  <- model.matrix(har_rr$fit)
sc <- sapply(seq_len(ncol(X)), function(j) Box.test(X[, j] * hr, lag = 10,
                                                    type = "Ljung-Box")$p.value)

# White (HC0), Newey-West se / OLS se -> 08_hac_decomposition
se_ols <- sqrt(diag(vcov(har_rr$fit)))
se_hc  <- sqrt(diag(sandwich::vcovHC(har_rr$fit, type = "HC0")))
se_hac <- sqrt(diag(har_rr$vcov_hac))
hac_dec <- data.frame(term = names(se_ols), se_OLS = se_ols,
                      ratio_HC_over_OLS = se_hc / se_ols,
                      ratio_HAC_over_OLS = se_hac / se_ols, row.names = NULL)
write.csv(hac_dec, "results/tables/08_hac_decomposition.csv")
print(hac_dec, digits = 4)

# Newey-West se / OLS se per lag -> 08_hac_bandwidth_profile: hac_lag | coefs
bw_prof <- do.call(rbind, lapply(c(1, 5, 10, 22, 40), function(L) {
  V <- NeweyWest(har_rr$fit, lag = L, prewhite = FALSE, adjust = TRUE)
  data.frame(hac_lag = L, t(round(sqrt(diag(V)) / se_ols, 4)))
}))
write.csv(bw_prof, "results/tables/08_hac_bandwidth_profile.csv")
print(bw_prof, digits = 4)

# Newey-West se per lag -> har_rr_se: hac_lag | se per coef
har_rr_se <- do.call(rbind, lapply(c(10, 22, 44), function(L) {
  V <- NeweyWest(har_rr$fit, lag = L, prewhite = FALSE, adjust = TRUE)
  data.frame(hac_lag = L, t(sqrt(diag(V))), check.names = FALSE)
}))

# refit other horizon sets, HAC lag = longest horizon -> 08_har_horizon_robustness
har_alt <- list(
  "1,5,22 (Corsi)" = har_rr,
  "1,5"            = fit_har(log(RR), p = c(1, 5)),
  "1,5,22,66"      = fit_har(log(RR), p = c(1, 5, 22, 66)),
  "1,22"           = fit_har(log(RR), p = c(1, 22)))
har_rob <- data.frame(spec = names(har_alt),
                      r2   = sapply(har_alt, `[[`, "r2"),
                      persistence = sapply(har_alt, `[[`, "persistence"),
                      n = sapply(har_alt, function(z) nrow(z$design)),
                      row.names = NULL)
write.csv(har_rob, "results/tables/08_har_horizon_robustness.csv")
print(har_rob, digits = 4)

# alignment checks (1-day shift inflates R^2 silently): 1st fit = day 23 -> first_fit_date
first_fit_date <- index(RR)[23]
stopifnot(length(fitted(har_rr$fit)) == nrow(RR) - 22)
stopifnot(index(har_rr$fitted_lvl)[1] == first_fit_date)
stopifnot(identical(as.numeric(index(har_rr$fitted_lvl)),
                    as.numeric(index(RR)[-(1:22)])))

## --- HAR vs GARCH
# Jensen: exp(fit + sigma^2/2) = mean (exp(fit) = median); merge() aligns dates -> cmp
har_jensen <- exp(har_rr$fitted_log + har_rr$sigma2 / 2)
cmp <- na.omit(merge(har_rr$fitted_lvl, har_jensen, sigma_garch))
colnames(cmp) <- c("HAR_naive_exp", "HAR_jensen", "GARCH_variance")

# GARCH/HAR level gap: Jensen factor + overnight share -> 08_scale_decomposition
jensen_factor <- exp(har_rr$sigma2 / 2)
overnight     <- qual["JPM", "overnight_share"]
pred_ratio    <- qual["JPM", "pred_garch_over_RR"]
ratio_tab <- data.frame(
  quantity = c("mean GARCH variance",
               "mean HAR fitted, naive exp() (a MEDIAN)",
               "mean HAR fitted, Jensen-corrected (a MEAN)",
               "ratio GARCH / HAR_naive",
               "ratio GARCH / HAR_jensen",
               "Jensen factor exp(sigma^2/2), predicted",
               "overnight variance share (02)",
               "predicted GARCH/HAR = 1 / intraday variance share"),
  value = c(mean(cmp$GARCH_variance), mean(cmp$HAR_naive_exp), mean(cmp$HAR_jensen),
            mean(cmp$GARCH_variance) / mean(cmp$HAR_naive_exp),
            mean(cmp$GARCH_variance) / mean(cmp$HAR_jensen),
            jensen_factor, overnight, pred_ratio))
write.csv(ratio_tab, "results/tables/08_scale_decomposition.csv")
print(ratio_tab, digits = 5)

# top: annualised vol GARCH vs HAR_jensen; bottom: log gap -> 08_har_vs_garch
png("results/figures/08_har_vs_garch.png", width = 1500, height = 950, res = 130)
par(mfrow = c(2, 1), mar = c(2.6, 4.2, 2.2, 1))
av <- sqrt(252 * cmp) * 100
plot(zoo::as.zoo(av$GARCH_variance), col = "steelblue", lwd = 0.9, xlab = "",
     ylab = "annualised sigma (%)",
     main = "GARCH conditional volatility (blue) vs Jensen-corrected HAR-RR (red)")
lines(zoo::as.zoo(av$HAR_jensen), col = "firebrick", lwd = 0.9)
grid(col = "grey88")
lg <- log(cmp$GARCH_variance) - log(cmp$HAR_jensen)
plot(zoo::as.zoo(lg), col = "grey25", xlab = "", ylab = "log GARCH - log HAR",
     main = "Log ratio: a roughly CONSTANT level gap, not a timing disagreement")
abline(h = mean(lg), col = "red", lty = 2)
abline(h = log(pred_ratio), col = "blue", lty = 3, lwd = 2)
grid(col = "grey88")
legend("topright", c("mean log gap", "predicted from overnight share"),
       lty = c(2, 3), col = c("red", "blue"), lwd = c(1, 2), bty = "n", cex = 0.8)
dev.off()

saveRDS(mget(c("cmp", "dw", "first_fit_date", "har_bg", "har_rr", "har_rr_se", "hr",
               "jensen_factor", "pred_ratio", "RR", "RR_raw", "sc", "ur_rr", "X")),
        "results/rds/08_realized_range_har.rds")
