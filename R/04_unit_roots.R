###############################################################################
# EBC4008 Case 3 -- part 04: unit roots (Pantula): y, log(y), r
# Reads : results/rds/03_exploration.rds (y, ly, rn)
# Writes: tables 04_lag_selection, 04_pantula, 04_kpss_bandwidth, 04_unitroot_summary,
#         04_integration_order; results/rds/04_unit_roots.rds
###############################################################################

library(bootUR)    # boot_adf()
library(tseries)   # kpss.test()

d03 <- readRDS("results/rds/03_exploration.rds")
y <- d03$y; ly <- d03$ly; rn <- d03$rn

SER       <- list(y = y, "log(y)" = ly, r = rn)
DET       <- c(y = "trend", "log(y)" = "trend", r = "intercept")   # boot_adf deterministics
KPSS_NULL <- c(y = "Trend", "log(y)" = "Trend", r = "Level")       # KPSS H0 on level

maxlag <- trunc(12 * (length(y) / 100)^(1 / 4))   # Schwert (1989) ceiling

# bootstrap ADF: AWB, 1999 reps, lags by `crit`, OLS detrending (default) -> htest
badf <- function(x, det, crit, label) {
  set.seed(1)
  boot_adf(x, data_name = label, bootstrap = "AWB", B = 1999, deterministics = det,
           max_lag = maxlag, criterion = crit, criterion_scale = FALSE, show_progress = FALSE)
}

# lags boot_adf picked -> integer
lags_used <- function(a) as.integer(a$details$`selected lags`)

# boot_adf on the level per lag criterion (MAIC/MBIC: Ng-Perron) -> lagfit, 04_lag_selection
lagfit <- list()
lagsel <- NULL
for (s in names(SER)) {
  lagfit[[s]] <- list()
  for (cr in c("AIC", "BIC", "MAIC", "MBIC")) {
    a <- badf(SER[[s]], DET[[s]], cr, label = s)
    lagfit[[s]][[cr]] <- a
    lagsel <- rbind(lagsel, data.frame(
      series = s, deterministics = DET[[s]], criterion = cr, ceiling = a$specifications$max_lag,
      lags_used = lags_used(a), tau = unname(a$statistic), p_boot = unname(a$p.value)))
  }
}
write.csv(lagsel, "results/tables/04_lag_selection.csv")
print(lagsel, digits = 3)

# Pantula: s1 = ADF on diff (H0 I(2)), s2 = ADF on level (H0 I(1)), MAIC lags -> ur[[s]]
# reject = bootstrap p < 5% (p = 0: below 1/1999) -> 04_pantula
ur <- list()
pantula <- NULL
for (s in names(SER)) {
  a1 <- badf(diff(SER[[s]]), "intercept", "MAIC", label = paste0("diff(", s, ")"))
  a2 <- lagfit[[s]]$MAIC
  ur[[s]] <- list(s1 = a1, s2 = a2)
  pantula <- rbind(pantula, data.frame(
    series     = s,
    step       = c("1: I(2) vs I(1)", "2: I(1) vs I(0)"),
    regression = c(paste0("ADF on diff(", s, "), intercept"), paste0("ADF on ", s, ", ", DET[[s]])),
    lags_used  = c(lags_used(a1), lags_used(a2)),
    tau        = unname(c(a1$statistic, a2$statistic)),
    p_boot     = unname(c(a1$p.value, a2$p.value))))
}
pantula$reject_5pct <- pantula$p_boot < 0.05
write.csv(pantula, "results/tables/04_pantula.csv")
print(pantula, digits = 3)

# KPSS on level and diff, lags trunc(4|12 (T/100)^(1/4)); short fits -> ur[[s]]$kp_lvl, $kp_dif
# p capped to [0.01, 0.10] (its warning muted) -> 04_kpss_bandwidth
kpss_tab <- NULL
for (s in names(SER)) for (L in c("short", "long")) {
  a <- suppressWarnings(kpss.test(SER[[s]],       null = KPSS_NULL[[s]], lshort = L == "short"))
  b <- suppressWarnings(kpss.test(diff(SER[[s]]), null = "Level",        lshort = L == "short"))
  a$data.name <- s; b$data.name <- paste0("diff(", s, ")")   # else prints "data: SER[[s]]"
  kpss_tab <- rbind(kpss_tab, data.frame(
    series = s, bandwidth = L, lags = unname(a$parameter), level_null = KPSS_NULL[[s]],
    kpss_level = unname(a$statistic), p_level = a$p.value,
    kpss_diff  = unname(b$statistic), p_diff  = b$p.value))
  if (L == "short") { ur[[s]]$kp_lvl <- a; ur[[s]]$kp_dif <- b }
}
write.csv(kpss_tab, "results/tables/04_kpss_bandwidth.csv")
print(kpss_tab, digits = 4)

# all tests, all series; reject = p < 5% -> ur_summary, 04_unitroot_summary
ur_summary <- NULL
for (s in names(SER)) {
  u <- ur[[s]]
  ur_summary <- rbind(ur_summary, data.frame(
    series    = s,
    test      = c(paste0("ADF diff(", s, "), intercept"), paste0("ADF ", s, ", ", DET[[s]]),
                  paste0("KPSS ", s, " (", KPSS_NULL[[s]], ")"),
                  paste0("KPSS diff(", s, ") (Level)")),
    statistic = unname(c(u$s1$statistic, u$s2$statistic, u$kp_lvl$statistic, u$kp_dif$statistic)),
    p_value   = unname(c(u$s1$p.value, u$s2$p.value, u$kp_lvl$p.value, u$kp_dif$p.value)),
    null      = c("2 unit roots", "1 unit root",
                  if (KPSS_NULL[[s]] == "Trend") "trend stationary" else "level stationary",
                  "difference is level stationary")))
}
ur_summary$reject_5pct <- ur_summary$p_value < 0.05
write.csv(ur_summary, "results/tables/04_unitroot_summary.csv")
print(ur_summary, digits = 4)

# d from Pantula; does short-bandwidth KPSS agree (NA for I(2)) -> integ, 04_integration_order
integ <- NULL
for (s in names(SER)) {
  p  <- pantula[pantula$series == s, ]
  k  <- kpss_tab[kpss_tab$series == s & kpss_tab$bandwidth == "short", ]
  d  <- if (!p$reject_5pct[1]) 2L else if (!p$reject_5pct[2]) 1L else 0L
  kl <- k$p_level < 0.05; kd <- k$p_diff < 0.05
  integ <- rbind(integ, data.frame(
    series = s, pantula_step1_reject = p$reject_5pct[1], pantula_step2_reject = p$reject_5pct[2],
    order_adf = paste0("I(", d, ")"), kpss_level_rejects = kl, kpss_diff_rejects = kd,
    kpss_agrees = if (d == 1L) kl && !kd else if (d == 0L) !kl else NA))
}
write.csv(integ, "results/tables/04_integration_order.csv")
print(integ, digits = 4)

saveRDS(mget(c("integ", "maxlag", "ur", "ur_summary")), "results/rds/04_unit_roots.rds")
