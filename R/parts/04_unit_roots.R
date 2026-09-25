###############################################################################
# EBC4008 Case 3 -- part 04: Unit roots: the Pantula principle
#
# Inputs come from need(),
# outputs go to results/parts/04_unit_roots.rds via part_save().
#
# Needs : ly, P, r, rn, y <- 03_exploration
# Makes : integ, maxlag, ur, ur_summary
# Writes: tables 04_lag_selection, 04_pantula, 04_deterministic_ladder,
#         04_kpss_bandwidth, 04_dfgls_lag_sensitivity, 04_zivot_andrews,
#         04_za_lag_sensitivity, 04_unitroot_summary, 04_integration_order; figures
#         04_zivot_andrews
#
# Every test runs on all three series: y (price), log(y) (log price) and
# r (log return).

# y is JPMorgan's (JPM) daily adjusted closing price, log(y) its natural log,
# r = diff(log(y)). Adjusted means past prices are corrected for stock splits
# and dividends.

# Two stories for a rising price
# 1. random walk with drift (a "unit root") -> Every shock stays in the price forever. There's no "normal" level to return to, and the range of possible prices keeps widening.
# 2. trend-stationary -> The price wobbles around a fixed upward line and is always pulled back to it, so shocks fade.

# under 2, a crashed price is expected to climb back to the trend line (predictable), under 1 the crash is permanent

# I(0) = the series has a home, shocks die out, forecast drifts back to the mean (uncertainty levels off)
# I(1) = memory never fades, spread keeps growing with the horizon (shocks are permanent) -> ARIMA gets d=1, model returns
# I(2) = even the returns have a random walk, so shocks even permanently change the slope
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

need("03_exploration", c("ly", "P", "r", "rn", "y"))

part_begin("04_unit_roots")

## ---- 4. Unit roots: the Pantula principle ---------------------------------------

say("4. UNIT ROOTS AND THE PANTULA PRINCIPLE")

# Pantula: start at the highest plausible order and work down. DF critical
# values assume exactly one unit root under H0. Each difference drops one
# deterministic term (trend -> constant -> nothing).
# test I(2) vs I(1), on reject I(1) vs I(0)

# The three series, and the deterministic terms their LEVEL regression gets:
# y and log(y) rise over 27 years -> start general, with a trend; r -> drift.
# The ladder further down shows the unit-root verdict doesn't depend on that:
# no tau rejects at any rung. For log(y) phi3/phi2/phi1 don't reject either
# (drift t = 1.37), so the trend is the conservative choice. For y the phi
# rejections come from a POSITIVE coefficient on the lagged level (tau3 +2.02,
# tau2 +4.21: exponential growth), not from a trend (trend t = 0.71).
SER       <- list(y = y, "log(y)" = ly, r = rn)
DET       <- c(y = "trend", "log(y)" = "trend", r = "drift")     # ur.df type for the level
KPSS_TYPE <- c(y = "tau",   "log(y)" = "tau",   r = "mu")        # KPSS H0: trend- / level-stationary
ERS_MODEL <- c(y = "trend", "log(y)" = "trend", r = "constant")  # DF-GLS deterministic part
DATES     <- list(y = index(P), "log(y)" = index(P), r = index(r))
KEY       <- c(y = "y", "log(y)" = "logy", r = "r")              # file-safe names (cache keys)

maxlag <- trunc(12 * (length(y) / 100)^(1 / 4))   # Schwert (1989) upper bound
stopifnot(trunc(12 * (length(rn) / 100)^(1 / 4)) == maxlag)   # r is one obs shorter: same ceiling

# lags an ur.df fit actually used: coefficient rows minus the deterministic
# terms and z.lag.1 (@lags holds the ceiling, not the choice)
lags_used <- function(a) nrow(a@testreg$coefficients) -
  switch(a@model, trend = 3L, drift = 2L, none = 1L)

# What does selectlags actually pick? At n ~ 6,900 AIC's penalty is nothing.
# Augmented Dickey-Fuller (ADF) unit root test, once with AIC
# ADF measures pullback with a regression (does yesterdays price help predict today's change)
#     - If the price is unusually high yesterday and tends to fall today (and vice versa), that's pull-back.
#     - If yesterday's level tells you nothing about today's change, there's no pull-back.
# lags need to be chosen carefully: too few  => pattern in the noise
#                                   too many => weaker at detecting pullback
#
# Where we choose the lag by AIC vs BIC
# - fit a regression 34 times 1lag, .... , 34lags
# - Give each fit a score, where lower is better (score = (how badly the model fits) + (penalty) × (number of coefficients))
# - Keeps the version with the lowest score.
lagsel <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c("AIC", "BIC"), function(crit) {
  # The result a is an object holding the test statistics, critical values and the winning regression.
  # ur.df: it fits the ADF regression 34 times (1 lag, 2 lags, … 34 lags), scores each fit with BIC (or AIC), keeps the best one, and calculates the test statistics from it.
  a <- ur.df(SER[[s]], type = DET[[s]], lags = maxlag, selectlags = crit)
  data.frame(series = s, type = DET[[s]], criterion = crit, ceiling = maxlag,
             lags_used = lags_used(a),
             tau = as.numeric(a@teststat)[1],   # the unit-root test (is the coefficient on the lagged level 0?)
             cv_5pct = a@cval[1, "5pct"])
}))))
tab("04_lag_selection", lagsel, digits = 3,
    caption = "Lags actually selected by ur.df under AIC vs BIC, per series")

# The fits every later table reads from, per series (all lags by AIC):
#   s1        step 1: ADF on the first difference, drift only
#   trend / drift / none: ADF on the level with each deterministic set
#   s2        step 2 = the level fit with DET's terms
# (diff(r) keeps a constant too: harmless, it only costs power -- tau -28.1
# against -2.86)
ur <- lapply(setNames(names(SER), names(SER)), function(s) {
  x <- SER[[s]]
  f <- list(s1    = ur.df(diff(x), type = "drift", lags = maxlag, selectlags = "AIC"),
            trend = ur.df(x, type = "trend", lags = maxlag, selectlags = "AIC"),
            drift = ur.df(x, type = "drift", lags = maxlag, selectlags = "AIC"),
            none  = ur.df(x, type = "none",  lags = maxlag, selectlags = "AIC"))
  f$s2 <- f[[DET[[s]]]]
  f
})

# Step 1: H0 I(2) vs I(1) -> ADF on the first difference, drift only
# Can we rule out I(2)? tau2 must be significantly NEGATIVE (below -2.86): ADF is left-tailed
# Step 2: H0 I(1) vs I(0) -> ADF on the level. For y and log(y) WITH trend (a
# random walk with drift has a trend in its mean, so the regression has to nest
# one); for r with drift only.
# I(1) or stationary around a trend? (again check tau against critical-value matrix)
pantula <- do.call(rbind, lapply(names(SER), function(s) {
  a1 <- ur[[s]]$s1; a2 <- ur[[s]]$s2
  data.frame(series = s,
             step = c("1: I(2) vs I(1)", "2: I(1) vs I(0)"),
             regression = c(paste0("ADF on diff(", s, "), drift"),
                            paste0("ADF on ", s, ", ", DET[[s]])),
             lags_used = c(lags_used(a1), lags_used(a2)),
             tau = c(as.numeric(a1@teststat)[1], as.numeric(a2@teststat)[1]),
             cv_5pct = c(a1@cval[1, "5pct"], a2@cval[1, "5pct"]),
             row.names = NULL)
}))
pantula$reject_5pct <- pantula$tau < pantula$cv_5pct
tab("04_pantula", pantula, digits = 3,
    caption = "Pantula: step 1 (two unit roots?) then step 2 (one unit root?), per series")

# Deterministic ladder (Dolado-Jenkinson-Sosvilla-Rivero): trend -> drift ->
# none, each tau read against its OWN critical values.
ladder <- do.call(rbind, lapply(names(SER), function(s) {
  ft <- ur[[s]]$trend; fd <- ur[[s]]$drift; fn <- ur[[s]]$none
  data.frame(series = s,
    rung = c("trend (tau3)", "trend (phi3)", "trend (phi2)",
             "drift (tau2)", "drift (phi1)", "none (tau1)"),
    statistic = c(as.numeric(ft@teststat)[1], as.numeric(ft@teststat)[3],
                  as.numeric(ft@teststat)[2], as.numeric(fd@teststat)[1],
                  as.numeric(fd@teststat)[2], as.numeric(fn@teststat)[1]),
    cv_5pct = c(ft@cval[1, "5pct"], ft@cval[3, "5pct"], ft@cval[2, "5pct"],
                fd@cval[1, "5pct"], fd@cval[2, "5pct"], fn@cval[1, "5pct"]),
    row.names = NULL)
}))
ladder$reject_5pct <- ifelse(grepl("phi", ladder$rung),
                             ladder$statistic > ladder$cv_5pct,   # F-type: right tail
                             ladder$statistic < ladder$cv_5pct)   # tau: left tail
tab("04_deterministic_ladder", ladder, digits = 3,
    caption = "The full deterministic-terms ladder (Dolado-Jenkinson-Sosvilla-Rivero), per series")

# KPSS flips the null (H0 = stationary), so agreement with ADF actually means
# something. ur.kpss lags default "short" (~11 here) vs Schwert 34 in the ADFs,
# and the bandwidth scales the statistic -> show both.
# KPSS fits a straight line through the series, H0 says it is stationary around that line
# add the gaps up over time. This running total is S(t) = gap(1) + gap(2) + … + gap(t).
#     If the series is stationary, gaps above and below the line alternate quickly and cancel out, so the running total stays small.
#     If it has a unit root, the price stays above (or below) the line for years at a stretch. The gaps all have the same sign, so the running total balloons.
# KPSS = (sum of S(t)² / n²)  ÷  long-run variance
# Level: trend-stationary (tau) for y and log(y), level-stationary (mu) for r.
# Difference: level-stationary (mu) for all three.
kpss_tab <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c("short", "long"), function(L) {
  a <- ur.kpss(SER[[s]],       type = KPSS_TYPE[[s]], lags = L)
  b <- ur.kpss(diff(SER[[s]]), type = "mu",           lags = L)
  data.frame(series = s, bandwidth = L, level_type = KPSS_TYPE[[s]],
             kpss_level = as.numeric(a@teststat)[1], cv5_level = a@cval[1, "5pct"],
             kpss_diff  = as.numeric(b@teststat)[1], cv5_diff  = b@cval[1, "5pct"])
}))))
tab("04_kpss_bandwidth", kpss_tab, digits = 4,
    caption = "KPSS under both urca bandwidth rules, level and first difference, per series")
for (s in names(SER)) {
  ur[[s]]$kp_lvl <- ur.kpss(SER[[s]],       type = KPSS_TYPE[[s]])   # H0: trend/level-stationary
  ur[[s]]$kp_dif <- ur.kpss(diff(SER[[s]]), type = "mu")             # H0: level-stationary
}

# DF-GLS (Elliott-Rothenberg-Stock): more power than ADF near the unit root.
# ur.ers lag.max is IMPOSED, not searched: passing 34 fits 34 lags and throws
# the power advantage away. Show a range instead.
# another unit-root test
# ADF has low power near the unit root, DF-GLS splits the job:
# 1. est the trendline (or the mean, for r) in an alternative way
# 2. remove that line and run DF test
# H0 = unit root vs H1 = stationary around a rising line (around a mean, for r).
# r: GLS demeaning puts almost all weight on the first obs, so the GLS "mean"
# is ~r[1] (-0.025, 1.1 sd below the mean) and the statistic hinges on that one
# day: with r[1] set to the mean, 4 lags gives -39.9 instead of -12.6, and 34
# lags -14.7 instead of -1.8. The 34-lag non-rejection is that artefact; ADF and
# KPSS settle r.
ers_tab <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c(1L, 2L, 4L, 8L, 12L, maxlag), function(L) {
  e <- ur.ers(SER[[s]], type = "DF-GLS", model = ERS_MODEL[[s]], lag.max = L)
  data.frame(series = s, model = ERS_MODEL[[s]], lags_imposed = L,
             statistic = as.numeric(e@teststat)[1],
             cv_1pct = e@cval[1, "1pct"], cv_5pct = e@cval[1, "5pct"],
             reject_5pct = as.numeric(e@teststat)[1] < e@cval[1, "5pct"])
}))))
tab("04_dfgls_lag_sensitivity", ers_tab, digits = 4,
    caption = "DF-GLS: lag.max is IMPOSED, not selected -- sensitivity to it, per series")
for (s in names(SER))
  ur[[s]]$ers <- ur.ers(SER[[s]], type = "DF-GLS", model = ERS_MODEL[[s]], lag.max = 4)

# Zivot-Andrews: unit-root test that allows one break. Five crises in the
# sample, and a break can make a stationary series look I(1) to ADF (Perron 1989).
# for every possible break date, fit an ADF-style equation that includes the break
# repeat for every date. urca tries every observation from 1 to 6,886 as the break date, which is 6,886 regressions.
# keep the most negative statistic
# compare with ZA's critical values
# H0: a unit root, with no break needed. H1: stationary around a trend line that breaks once.
# Model "both" (break in intercept and trend) for y and log(y); "intercept" for
# r, which has no trend. r's ZA profile is flat anyway, so its break date isn't
# identified -- only the rejection matters there.
# lag is imposed (ur.za doesn't select it): 12, with a lag sensitivity for
# log(y) below. ~45 s per fit, so cached; the key spells out the spec and n so
# a changed spec or data can't pick up an old fit.
ZA_MODEL <- c(y = "both", "log(y)" = "both", r = "intercept")
ZA_LAG   <- 12L
za_fit <- function(s, lag) cache(sprintf("za_%s_%s_L%d_n%d", KEY[[s]], ZA_MODEL[[s]], lag, length(SER[[s]])),
                                 ur.za(SER[[s]], model = ZA_MODEL[[s]], lag = lag))
for (s in names(SER)) ur[[s]]$za <- za_fit(s, ZA_LAG)
stopifnot(all(sapply(names(SER), function(s) length(ur[[s]]$za@y) == length(SER[[s]]))))
za_tab <- do.call(rbind, lapply(names(SER), function(s) {
  z <- ur[[s]]$za
  data.frame(series = s, model = z@model, lag = z@lag, statistic = as.numeric(z@teststat)[1],
             cv_1pct = z@cval[1], cv_5pct = z@cval[2], cv_10pct = z@cval[3],  # unnamed (1%, 5%, 10%)
             reject_5pct = as.numeric(z@teststat)[1] < z@cval[2],
             break_obs = z@bpoint, break_date = format(DATES[[s]][z@bpoint]),
             row.names = NULL)
}))
tab("04_zivot_andrews", za_tab, digits = 3,
    caption = sprintf("Zivot-Andrews (%d lags; model %s), per series", ZA_LAG,
                      paste(ZA_MODEL, collapse = " / ")))

# The imposed lag matters for log(y): it rejects at 5% only with 1 lag (at 10%
# with 1 and 4), not with 8 or more; the break date stays 2008-10-02. So ZA is
# borderline: a break-stationary alternative can't be firmly ruled out.
za_lag <- do.call(rbind, lapply(c(1L, 4L, 8L, ZA_LAG, maxlag), function(L) {
  z <- if (L == ZA_LAG) ur[["log(y)"]]$za else za_fit("log(y)", L)
  data.frame(lag = L, statistic = as.numeric(z@teststat)[1], cv_5pct = z@cval[2],
             cv_10pct = z@cval[3], reject_5pct = as.numeric(z@teststat)[1] < z@cval[2],
             reject_10pct = as.numeric(z@teststat)[1] < z@cval[3],
             break_date = format(index(P)[z@bpoint]))
}))
tab("04_za_lag_sensitivity", za_lag, digits = 3,
    caption = "Zivot-Andrews on log(y), model both: sensitivity to the imposed lag")

# Drawn by hand: urca's plot(ur.za) labels its lines "1% / 2.5% / 5%" but draws
# the 1%, 5% and 10% values, and has no dates on the x-axis.
# The line is the ZA statistic for each candidate break date; the test takes
# its minimum, and rejects if that minimum falls below a critical line.
ZA_COL <- c("red", "darkorange", "forestgreen")     # 1%, 5%, 10% critical values
fig("04_zivot_andrews", {
  par(mfrow = c(3, 1), mar = c(3, 4.5, 3, 1), oma = c(0, 0, 2, 0))
  for (s in names(SER)) {
    z <- ur[[s]]$za
    d <- DATES[[s]][seq_along(z@tstats)]
    plot(d, z@tstats, type = "l", xlab = "", ylab = "ZA t-statistic",
         ylim = range(c(z@tstats, z@cval), na.rm = TRUE),
         main = sprintf("%s: minimum %.2f on %s (5%% critical value %.2f)",
                        s, as.numeric(z@teststat)[1], format(DATES[[s]][z@bpoint]), z@cval[2]))
    for (i in 1:3) abline(h = z@cval[i], col = ZA_COL[i], lty = c(2, 5, 3)[i], lwd = 2)
    abline(v = DATES[[s]][z@bpoint], col = "grey50", lty = 3)
  }
  # one legend for all three panels, in the outer top margin
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("top", c("1% critical value", "5% critical value", "10% critical value",
                  "chosen break date"),
         col = c(ZA_COL, "grey50"), lty = c(2, 5, 3, 3), lwd = c(2, 2, 2, 1), horiz = TRUE, bty = "n")
}, height = 1300)

# every test on every series in one place
ur_summary <- do.call(rbind, lapply(names(SER), function(s) {
  u <- ur[[s]]
  data.frame(series = s,
    test = c(paste0("ADF diff(", s, "), drift (tau2)"), paste0("ADF ", s, ", trend (tau3)"),
             paste0("ADF ", s, ", drift (tau2)"), paste0("KPSS ", s, " (", KPSS_TYPE[[s]], ")"),
             paste0("KPSS diff(", s, ") (mu)"), paste0("DF-GLS ", s, ", ", ERS_MODEL[[s]], ", 4 lags"),
             paste0("Zivot-Andrews ", s, ", ", ZA_MODEL[[s]])),
    statistic = c(as.numeric(u$s1@teststat)[1], as.numeric(u$trend@teststat)[1],
                  as.numeric(u$drift@teststat)[1], as.numeric(u$kp_lvl@teststat)[1],
                  as.numeric(u$kp_dif@teststat)[1], as.numeric(u$ers@teststat)[1],
                  as.numeric(u$za@teststat)[1]),
    cv_5pct = c(u$s1@cval[1, "5pct"], u$trend@cval[1, "5pct"], u$drift@cval[1, "5pct"],
                u$kp_lvl@cval[1, "5pct"], u$kp_dif@cval[1, "5pct"], u$ers@cval[1, "5pct"],
                u$za@cval[2]),
    null = c("2 unit roots", "1 unit root", "1 unit root",
             if (KPSS_TYPE[[s]] == "tau") "trend stationary" else "level stationary",
             "difference is level stationary", "1 unit root",
             "1 unit root, no break (H1: stationary + 1 break)"),
    row.names = NULL)
}))
ur_summary$reject_5pct <- ifelse(
  grepl("KPSS", ur_summary$test),
  ur_summary$statistic > ur_summary$cv_5pct,      # KPSS rejects on the RIGHT tail
  ur_summary$statistic < ur_summary$cv_5pct)      # ADF-type reject on the LEFT tail
tab("04_unitroot_summary", ur_summary, caption = "All unit-root evidence in one place, per series")

# Order of integration: Pantula's answer, then whether KPSS agrees with it.
# KPSS agrees with I(1) if it rejects stationarity of the level but not of the
# difference; with I(0) if it doesn't reject the level.
# y is the odd one out: ADF can't reject a unit root, but KPSS also rejects
# for diff(y). Not a second unit root: a dollar move is a % move times the
# price, so the MEAN of diff(y) grows with the level ($0.002/day in the first
# half, $0.09 in the second), and KPSS(mu) tests for a constant mean. The
# growing spread alone doesn't do it: KPSS on the zero-mean part
# P[t-1]*(u - mean(u)) is 0.146 < 0.463. log(y) removes both (its difference,
# r, passes ADF and KPSS) -- the reason the models use log(y) and r.
# And y's taus are POSITIVE (+2.02 / +4.21 / +5.03; DF-GLS about +1): a high
# price predicts a RISE, exponential growth in dollars. ADF only tests the left
# tail, so "I(1)" for y means "not stationary"; in dollars the price is
# explosive. Another reason to model log(y).
integ <- do.call(rbind, lapply(names(SER), function(s) {
  p  <- pantula[pantula$series == s, ]
  k  <- kpss_tab[kpss_tab$series == s & kpss_tab$bandwidth == "short", ]
  d  <- if (!p$reject_5pct[1]) 2L else if (!p$reject_5pct[2]) 1L else 0L
  kl <- k$kpss_level > k$cv5_level; kd <- k$kpss_diff > k$cv5_diff
  data.frame(series = s, pantula_step1_reject = p$reject_5pct[1],
             pantula_step2_reject = p$reject_5pct[2],
             order_adf = paste0("I(", d, ")"),
             kpss_level_rejects = kl, kpss_diff_rejects = kd,
             kpss_agrees = if (d == 1L) kl && !kd else if (d == 0L) !kl else NA,
             row.names = NULL)
}))
tab("04_integration_order", integ,
    caption = "Order of integration per series: Pantula (ADF), checked against KPSS")

# urca fits keep the function frame they ran in (lag matrix, helper fits): drop it before saving
ur <- lapply(ur, function(u) lapply(u, slim))

part_save(c("integ", "maxlag", "ur", "ur_summary"))
