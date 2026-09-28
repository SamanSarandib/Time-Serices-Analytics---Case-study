###############################################################################
# EBC4008 Case 3 -- part 04: unit roots (Pantula principle) on y, log(y), r
# Needs : ly, rn, y <- 03_exploration
# Makes : ur (per-series list of urca fits), ur_summary, integ (data.frames),
#         maxlag (integer ADF lag ceiling)
# Writes: tables 04_lag_selection, 04_pantula, 04_deterministic_ladder,
#         04_kpss_bandwidth, 04_unitroot_summary, 04_integration_order
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

part_begin("04_unit_roots")

# the three series + per-series test specs (y, log(y): trend; r: constant only)
SER       <- list(y = y, "log(y)" = ly, r = rn)
DET       <- c(y = "trend", "log(y)" = "trend", r = "drift")     # ur.df type for the level
KPSS_TYPE <- c(y = "tau",   "log(y)" = "tau",   r = "mu")        # KPSS H0: trend- / level-stationary

maxlag <- trunc(12 * (length(y) / 100)^(1 / 4))   # Schwert (1989) upper bound

# lags an ur.df fit actually picked (@lags holds the ceiling, not the choice) -> integer
lags_used <- function(a) nrow(a@testreg$coefficients) -
  switch(a@model, trend = 3L, drift = 2L, none = 1L)

# ADF (ur.df) with lags chosen by AIC vs BIC up to maxlag, per series
# -> 04_lag_selection: series | type | criterion | ceiling | lags_used | tau | cv_5pct
lagsel <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c("AIC", "BIC"), function(crit) {
  a <- ur.df(SER[[s]], type = DET[[s]], lags = maxlag, selectlags = crit)
  data.frame(series = s, type = DET[[s]], criterion = crit, ceiling = maxlag,
             lags_used = lags_used(a),
             tau = as.numeric(a@teststat)[1],   # [1] = tau, the rest are phi F-stats
             cv_5pct = a@cval[1, "5pct"])
}))))
tab("04_lag_selection", lagsel, digits = 3,
    caption = "Lags actually selected by ur.df under AIC vs BIC, per series")

# ur.df per series, lags by AIC; s2 = whichever level fit DET picks
# -> ur[[s]]: list(s1 = ADF on diff (drift), trend, drift, none = ADF on the level, s2)
ur <- lapply(setNames(names(SER), names(SER)), function(s) {
  x <- SER[[s]]
  f <- list(s1    = ur.df(diff(x), type = "drift", lags = maxlag, selectlags = "AIC"),
            trend = ur.df(x, type = "trend", lags = maxlag, selectlags = "AIC"),
            drift = ur.df(x, type = "drift", lags = maxlag, selectlags = "AIC"),
            none  = ur.df(x, type = "none",  lags = maxlag, selectlags = "AIC"))
  f$s2 <- f[[DET[[s]]]]
  f
})

# Pantula: step 1 ADF on the diff (H0 I(2)), step 2 ADF on the level (H0 I(1)), left-tailed
# -> 04_pantula: series | step | regression | lags_used | tau | cv_5pct | reject_5pct
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

# Dolado-Jenkinson-Sosvilla-Rivero ladder: trend -> drift -> none, each stat vs its own cv
# -> 04_deterministic_ladder: series | rung | statistic | cv_5pct | reject_5pct
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

# KPSS (ur.kpss, H0 stationary) on the level and the diff, under both urca bandwidths
# -> 04_kpss_bandwidth: series | bandwidth | level_type | kpss_level | cv5_level | kpss_diff | cv5_diff
kpss_tab <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c("short", "long"), function(L) {
  a <- ur.kpss(SER[[s]],       type = KPSS_TYPE[[s]], lags = L)
  b <- ur.kpss(diff(SER[[s]]), type = "mu",           lags = L)
  data.frame(series = s, bandwidth = L, level_type = KPSS_TYPE[[s]],
             kpss_level = as.numeric(a@teststat)[1], cv5_level = a@cval[1, "5pct"],
             kpss_diff  = as.numeric(b@teststat)[1], cv5_diff  = b@cval[1, "5pct"])
}))))
tab("04_kpss_bandwidth", kpss_tab, digits = 4,
    caption = "KPSS under both urca bandwidth rules, level and first difference, per series")
# default ("short") KPSS fits -> ur[[s]]$kp_lvl, ur[[s]]$kp_dif
for (s in names(SER)) {
  ur[[s]]$kp_lvl <- ur.kpss(SER[[s]],       type = KPSS_TYPE[[s]])   # H0: trend/level-stationary
  ur[[s]]$kp_dif <- ur.kpss(diff(SER[[s]]), type = "mu")             # H0: level-stationary
}

# every test on every series in one table (ADF, KPSS)
# -> 04_unitroot_summary: series | test | statistic | cv_5pct | null | reject_5pct
ur_summary <- do.call(rbind, lapply(names(SER), function(s) {
  u <- ur[[s]]
  data.frame(series = s,
    test = c(paste0("ADF diff(", s, "), drift (tau2)"), paste0("ADF ", s, ", trend (tau3)"),
             paste0("ADF ", s, ", drift (tau2)"), paste0("KPSS ", s, " (", KPSS_TYPE[[s]], ")"),
             paste0("KPSS diff(", s, ") (mu)")),
    statistic = c(as.numeric(u$s1@teststat)[1], as.numeric(u$trend@teststat)[1],
                  as.numeric(u$drift@teststat)[1], as.numeric(u$kp_lvl@teststat)[1],
                  as.numeric(u$kp_dif@teststat)[1]),
    cv_5pct = c(u$s1@cval[1, "5pct"], u$trend@cval[1, "5pct"], u$drift@cval[1, "5pct"],
                u$kp_lvl@cval[1, "5pct"], u$kp_dif@cval[1, "5pct"]),
    null = c("2 unit roots", "1 unit root", "1 unit root",
             if (KPSS_TYPE[[s]] == "tau") "trend stationary" else "level stationary",
             "difference is level stationary"),
    row.names = NULL)
}))
ur_summary$reject_5pct <- ifelse(
  grepl("KPSS", ur_summary$test),
  ur_summary$statistic > ur_summary$cv_5pct,      # KPSS rejects on the RIGHT tail
  ur_summary$statistic < ur_summary$cv_5pct)      # ADF rejects on the LEFT tail
tab("04_unitroot_summary", ur_summary, caption = "All unit-root evidence in one place, per series")

# order of integration: d from Pantula, then whether short-bandwidth KPSS agrees (NA for I(2))
# -> 04_integration_order: series | pantula_step1_reject | pantula_step2_reject | order_adf |
#    kpss_level_rejects | kpss_diff_rejects | kpss_agrees
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

# -> results/parts/04_unit_roots.rds
part_save(c("integ", "maxlag", "ur", "ur_summary"))
