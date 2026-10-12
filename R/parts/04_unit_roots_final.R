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

SER    <- list(y = y, "log(y)" = ly, r = rn)
maxlag <- trunc(12 * (length(y) / 100)^(1 / 4))        # Schwert (1989) = bootUR's default max_lag
LEVEL  <- 0.05

# Start at the highest order and work down, one unit root at a time:
#   test diff(x, 2)  H0: I(3) vs I(2)    intercept
#   test diff(x, 1)  H0: I(2) vs I(1)    intercept
#   test x           H0: I(1) vs I(0)    trend for y, log(y) (they rise); intercept for r
# The order is the first step from the top where the unit root is NOT rejected
# (all three rejected -> I(0)).
# bootUR::adf = standard ADF, lags by modified AIC (MAIC), MacKinnon p-value.
LEVEL_DET <- c(y = "trend", "log(y)" = "trend", r = "intercept")

ur <- lapply(setNames(names(SER), names(SER)), function(s) {
  x <- SER[[s]]
  list(d2 = adf(diff(x, differences = 2), deterministics = "intercept"),
       d1 = adf(diff(x, differences = 1), deterministics = "intercept"),
       d0 = adf(x,                        deterministics = LEVEL_DET[[s]]))
})

# -> 04_pantula: series | step | H0 | tested_on | deterministics | statistic | p_value | reject_5pct
pantula <- do.call(rbind, lapply(names(SER), function(s) {
  data.frame(series = s, step = 1:3,
             H0 = c("I(3) vs I(2)", "I(2) vs I(1)", "I(1) vs I(0)"),
             tested_on = c(paste0("diff(", s, ", 2)"), paste0("diff(", s, ")"), s),
             deterministics = c("intercept", "intercept", LEVEL_DET[[s]]),
             statistic = sapply(ur[[s]], function(o) unname(o$statistic)),
             p_value   = sapply(ur[[s]], function(o) unname(o$p.value)),
             row.names = NULL)
}))
pantula$reject_5pct <- pantula$p_value < LEVEL
tab("04_pantula", pantula, digits = 4,
    caption = "Pantula principle with bootUR::adf: from the second difference down to the level")

# order = first non-rejection from the top
# -> integ: series | order | d   (KPSS columns added in section 2)
integ <- do.call(rbind, lapply(names(SER), function(s) {
  rej <- pantula$reject_5pct[pantula$series == s]          # steps 1, 2, 3
  d   <- if (!rej[1]) 3L else if (!rej[2]) 2L else if (!rej[3]) 1L else 0L
  data.frame(series = s, order = paste0("I(", d, ")"), d = d, row.names = NULL)
}))

## ---- 2. KPSS: the opposite null as a check ---------------------------------
# KPSS flips the hypotheses (H0 = stationary), so agreement with Pantula means
# something. Level: trend-stationary (tau) for y, log(y); level-stationary (mu)
# for r. First difference: level-stationary (mu). Both urca bandwidths, since
# the bandwidth scales the statistic. KPSS rejects on the RIGHT tail.
# -> 04_kpss_bandwidth: series | bandwidth | level_type | kpss_level | cv5_level | kpss_diff | cv5_diff
KPSS_TYPE <- c(y = "tau", "log(y)" = "tau", r = "mu")
kpss_tab <- do.call(rbind, lapply(names(SER), function(s) do.call(rbind, lapply(c("short", "long"), function(L) {
  a <- ur.kpss(SER[[s]],       type = KPSS_TYPE[[s]], lags = L)
  b <- ur.kpss(diff(SER[[s]]), type = "mu",           lags = L)
  data.frame(series = s, bandwidth = L, level_type = KPSS_TYPE[[s]],
             kpss_level = as.numeric(a@teststat)[1], cv5_level = a@cval[1, "5pct"],
             kpss_diff  = as.numeric(b@teststat)[1], cv5_diff  = b@cval[1, "5pct"])
}))))
tab("04_kpss_bandwidth", kpss_tab, digits = 4,
    caption = "KPSS under both urca bandwidth rules, level and first difference, per series")

# does KPSS (short bandwidth) agree with Pantula's order?
# I(1): rejects stationarity of the level but not of the difference; I(0): doesn't reject the level
k_short <- kpss_tab[kpss_tab$bandwidth == "short", ]
integ$kpss_level_rejects <- k_short$kpss_level > k_short$cv5_level
integ$kpss_diff_rejects  <- k_short$kpss_diff  > k_short$cv5_diff
integ$kpss_agrees <- ifelse(integ$d == 1L, integ$kpss_level_rejects & !integ$kpss_diff_rejects,
                            ifelse(integ$d == 0L, !integ$kpss_level_rejects, NA))
tab("04_integration_order", integ,
    caption = "Order of integration per series: Pantula (ADF), checked against KPSS")
# y: KPSS also rejects for diff(y). Not a second unit root: a dollar move is a
# % move times the price, so the mean of diff(y) grows with the level and
# KPSS(mu) tests for a constant mean. log(y) removes this -- why we model log(y).

# Pantula used deterministics CHOSEN by reasoning. Now that d is known, we will check them.

# (a) Test them directly on the stationary transform (diff(x) for I(1), x for I(0)).
# No unit root left there, so ordinary t-tests are valid; HAC (Newey-West) se
# because volatility clustering makes the iid se too small.
#   const: I(1) -> the drift of x (mean of diff(x));  I(0) -> the mean of x
#   trend: does that drift / mean change over time? (time rescaled to (0, 1])
# log(y) and r give identical rows: diff(log(y)) IS r.
# Newey-West lag by the standard rule of thumb floor(4 (T/100)^(2/9)) = 10 here.
# (Checked: lags 0 and 22 give the same verdict, t = 1.37 and 1.65 for the drift.)
# -> 04_deterministics: series | order | tested_on | const | t_const | p_const |
#    trend | t_trend | p_trend | sig_const | sig_trend
HAC_LAG <- floor(4 * (length(rn) / 100)^(2 / 9))
d_of <- function(s) integ$d[integ$series == s]

det_tab <- do.call(rbind, lapply(names(SER), function(s) {
  d   <- d_of(s)
  z   <- if (d >= 1L) diff(SER[[s]], differences = d) else SER[[s]]
  tt  <- seq_along(z) / length(z)
  m0  <- lm(z ~ 1)                     # constant only: drift / mean
  m1  <- lm(z ~ tt)                    # constant + time: trend in it?
  se0 <- sqrt(diag(NeweyWest(m0, lag = HAC_LAG, prewhite = FALSE, adjust = TRUE)))
  se1 <- sqrt(diag(NeweyWest(m1, lag = HAC_LAG, prewhite = FALSE, adjust = TRUE)))
  t0  <- unname(coef(m0)[1] / se0[1])
  t1  <- unname(coef(m1)[2] / se1[2])
  data.frame(series = s, order = paste0("I(", d, ")"),
             tested_on = if (d >= 1L) paste0("diff(", s, ")") else s,
             const = unname(coef(m0)[1]), t_const = t0, p_const = 2 * pnorm(-abs(t0)),
             trend = unname(coef(m1)[2]), t_trend = t1, p_trend = 2 * pnorm(-abs(t1)),
             row.names = NULL)
}))
det_tab$sig_const <- stars(det_tab$p_const)
det_tab$sig_trend <- stars(det_tab$p_trend)
tab("04_deterministics", det_tab, digits = 5,
    caption = "Deterministic terms tested on the stationary transform, HAC (Newey-West) t-stats")

# Retest I(1) vs I(0) with the deterministics we now know.
# Pantula's level test used terms chosen by reasoning (trend for y, log(y)).
# (a) tells us which terms the data actually support:
#   I(1): drift significant -> trend, else intercept
#   I(0): trend significant -> trend, mean significant -> intercept, else none
# An unneeded term doesn't invalidate ADF but costs power, so redo the last
# Pantula step (the level, H0: I(1) vs H1: I(0)) with the supported terms and
# check that the conclusion is the same.
# -> 04_deterministics_adf: series | det_pantula | p_pantula | det_supported |
#    stat_supported | p_supported | reject_5pct | same_conclusion
det_adf <- do.call(rbind, lapply(names(SER), function(s) {
  k   <- det_tab[det_tab$series == s, ]
  cs  <- k$p_const < LEVEL; ts <- k$p_trend < LEVEL
  sup <- if (d_of(s) >= 1L) { if (cs) "trend" else "intercept" } else
  { if (ts) "trend" else if (cs) "intercept" else "none" }
  o   <- adf(SER[[s]], deterministics = sup)
  p0  <- pantula$p_value[pantula$series == s & pantula$step == 3]
  data.frame(series = s,
             det_pantula = LEVEL_DET[[s]], p_pantula = p0,
             det_supported = sup,
             stat_supported = unname(o$statistic), p_supported = unname(o$p.value),
             reject_5pct = unname(o$p.value) < LEVEL,
             same_conclusion = (p0 < LEVEL) == (unname(o$p.value) < LEVEL),
             row.names = NULL)
}))
tab("04_deterministics_adf", det_adf, digits = 4,
    caption = "Retest of I(1) vs I(0) on the level with the data-supported deterministics")

# -> results/parts/04_unit_roots.rds
part_save(c("det_adf", "det_tab", "integ", "kpss_tab", "maxlag", "pantula", "ur"))