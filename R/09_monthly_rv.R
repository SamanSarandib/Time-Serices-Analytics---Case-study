###############################################################################
# EBC4008 Case 3 -- part 09: monthly RV, HAR-RV
# Reads : results/rds/03_exploration.rds (r)
#         results/rds/08_realized_range_har.rds (har_rr)
# Writes: tables 09_rv_fss_comparison, 09_rv_unit_root, 09_har_coefs, 09_monthly_har_robustness,
#         09_har_daily_vs_monthly; figures 09_monthly_rv, 09_monthly_har_fit;
#         results/rds/09_monthly_rv.rds
###############################################################################

library(xts)        # xts(), index(), merge(); also attaches zoo
library(sandwich)   # NeweyWest()
library(bootUR)     # boot_adf()

d03 <- readRDS("results/rds/03_exploration.rds")
r <- d03$r
d08 <- readRDS("results/rds/08_realized_range_har.rds")
har_rr <- d08$har_rr

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

## --- Monthly RV
# RV_m = sum_{d in m} r_d^2 (sum, not mean; not demeaned)
# per month: RV, FSS-corrected RV, #days; months < 15 days dropped, dated at last day -> RVm_all
rn  <- as.numeric(r)
key <- format(index(r), "%Y-%m")
rv    <- tapply(rn, key, function(z) sum(z^2))
fss   <- tapply(rn, key, function(z) sum(z^2) + 2 * sum(z[-length(z)] * z[-1]))
n_obs <- tapply(rn, key, length)
dmax  <- tapply(index(r), key, max)                   # tapply drops the Date class
keep  <- n_obs >= 15
RVm_all <- xts(cbind(RVm = as.numeric(rv[keep]), RVm_fss = as.numeric(fss[keep]),
                     nobs = as.numeric(n_obs[keep])),
               order.by = as.Date(dmax[keep], origin = "1970-01-01"))
RV_m <- RVm_all[, "RVm"]
dropped  <- names(n_obs)[!keep]
nobs_rng <- range(RVm_all$nobs)

# relative se of RV_m = sqrt(2 / N), iid normal, N = mean days/month -> rv_rel_se
rv_rel_se <- sqrt(2 / mean(RVm_all$nobs))

# French-Schwert-Stambaugh adds 2 sum r_t r_{t+1}; plain RV used -> 09_rv_fss_comparison
fsscmp <- data.frame(
  measure = c("RV plain", "RV with FSS correction"),
  mean = c(mean(RVm_all$RVm), mean(RVm_all$RVm_fss)),
  sd   = c(sd(RVm_all$RVm),   sd(RVm_all$RVm_fss)),
  min  = c(min(RVm_all$RVm),  min(RVm_all$RVm_fss)),
  n_nonpositive = c(sum(RVm_all$RVm <= 0), sum(RVm_all$RVm_fss <= 0)))
write.csv(fsscmp, "results/tables/09_rv_fss_comparison.csv")
print(fsscmp, digits = 6)

# RV_m, log(RV_m) -> 09_monthly_rv
png("results/figures/09_monthly_rv.png", width = 1500, height = 900, res = 130)
par(mfrow = c(2, 1), mar = c(2.4, 4, 2, 1))
plot(zoo::as.zoo(RV_m), xlab = "", ylab = "", col = "grey20",
     main = "Monthly realized variance (sum of squared daily log returns)")
grid(col = "grey88")
plot(zoo::as.zoo(log(RV_m)), xlab = "", ylab = "", col = "grey20",
     main = "log monthly realized variance -- the modelling scale")
grid(col = "grey88")
dev.off()

## --- Is log(RV_m) I(0)?
# bootstrap ADF: intercept, MAIC, AWB p -> ur_rvm, 09_rv_unit_root
maxlag <- trunc(12 * (nrow(RV_m) / 100)^(1 / 4))
set.seed(1)
ur_rvm <- boot_adf(as.numeric(log(RV_m)), data_name = "log(RV_m)", bootstrap = "AWB", B = 1999,
                   deterministics = "intercept", max_lag = maxlag, criterion = "MAIC",
                   criterion_scale = FALSE, show_progress = FALSE)
rvm_ur <- data.frame(series = "log(RV_m)", deterministics = "intercept", criterion = "MAIC",
                     ceiling = maxlag, lags_used = as.integer(ur_rvm$details$`selected lags`),
                     tau = unname(ur_rvm$statistic), p_boot = unname(ur_rvm$p.value))
rvm_ur$reject_5pct <- rvm_ur$p_boot < 0.05
write.csv(rvm_ur, "results/tables/09_rv_unit_root.csv")
print(rvm_ur, digits = 3)

## --- Monthly HAR
# log RV_m on 1/3/12-month means, Newey-West 12 lags -> har_m; se at other lags -> har_m_se
har_m <- fit_har(log(RV_m), p = c(1, 3, 12))
har_m_se <- do.call(rbind, lapply(c(3, 6, 12), function(L) {
  V <- NeweyWest(har_m$fit, lag = L, prewhite = FALSE, adjust = TRUE)
  data.frame(hac_lag = L, t(sqrt(diag(V))), check.names = FALSE)
}))

# coefficients, Newey-West se, stars -> 09_har_coefs
harm_coefs <- data.frame(term = rownames(har_m$hac), estimate = har_m$hac[, 1],
                         se_hac = har_m$hac[, 2], t_hac = har_m$hac[, 3],
                         p_hac = har_m$hac[, 4], row.names = NULL)
pv <- harm_coefs$p_hac
harm_coefs$sig <- ifelse(pv < .01, "***", ifelse(pv < .05, "**", ifelse(pv < .10, "*", "")))
write.csv(harm_coefs, "results/tables/09_har_coefs.csv")
print(harm_coefs, digits = 4)

# refit other horizon sets, HAC lag = longest horizon -> mh_tab, 09_monthly_har_robustness
har_m_alt <- list(
  "1,3,12 (analogy)" = har_m,
  "1,3"              = fit_har(log(RV_m), p = c(1, 3)),
  "1,6"              = fit_har(log(RV_m), p = c(1, 6)),
  "1,3,6"            = fit_har(log(RV_m), p = c(1, 3, 6)),
  "1,12"             = fit_har(log(RV_m), p = c(1, 12)))
mh_tab <- data.frame(spec = names(har_m_alt),
                     n    = sapply(har_m_alt, function(z) nrow(z$design)),
                     r2   = sapply(har_m_alt, `[[`, "r2"),
                     adj_r2 = sapply(har_m_alt, function(z) summary(z$fit)$adj.r.squared),
                     persistence = sapply(har_m_alt, `[[`, "persistence"),
                     row.names = NULL)
write.csv(mh_tab, "results/tables/09_monthly_har_robustness.csv")
print(mh_tab, digits = 4)

# daily HAR-RR vs monthly HAR-RV -> 09_har_daily_vs_monthly
har_cmp <- data.frame(
  model = c("daily HAR-RR (1/5/22)", "monthly HAR-RV (1/3/12)"),
  n = c(nrow(har_rr$design), nrow(har_m$design)),
  r2 = c(har_rr$r2, har_m$r2),
  persistence = c(har_rr$persistence, har_m$persistence),
  proxy_rel_se_pct = c(NA, rv_rel_se * 100))
write.csv(har_cmp, "results/tables/09_har_daily_vs_monthly.csv")
print(har_cmp, digits = 4)

# actual vs one-step fitted log RV_m -> 09_monthly_har_fit
png("results/figures/09_monthly_har_fit.png", width = 1500, height = 950, res = 130)
par(mar = c(3, 4.2, 2.5, 1))
ac <- na.omit(merge(log(RV_m), har_m$fitted_log))
colnames(ac) <- c("actual", "fitted")
plot(zoo::as.zoo(ac[, 1]), col = "grey45", lwd = 1, xlab = "", ylab = "log RV",
     main = "Monthly HAR: actual (grey) vs one-step-ahead fitted (red)")
lines(zoo::as.zoo(ac[, 2]), col = "firebrick", lwd = 1.6)
grid(col = "grey88")
dev.off()

saveRDS(mget(c("dropped", "har_m", "har_m_se", "mh_tab", "nobs_rng", "rv_rel_se",
               "RVm_all", "ur_rvm")), "results/rds/09_monthly_rv.rds")
