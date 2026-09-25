###############################################################################
# EBC4008 Case 3 -- part 03: Data exploration
#
# Inputs come from need(),
# outputs go to results/parts/03_exploration.rds via part_save().
#
# Needs : jpm <- 01_data
# Makes : arch, lbr, ly, P, r, ra_r, rn, y
# Writes: tables 03_descriptive_stats, 03_stylized_facts, 03_robust_acf,
#         03_ljungbox_robust, 03_calendar_effects, 03_annual_volatility;
#         figures 03_price_return_panels, 03_logprice_episodes,
#         03_distribution_series, 03_boxplot_by_year, 03_rolling_stats,
#         03_acf_contrast, 03_distribution, 03_boxplot_calendar,
#         03_annual_volatility
#
# Notation from here on: y = price (adjusted close), log(y) = log price (ly in
# the code), r = log return = diff(log(y)).

# This script turns JPMs price into series usable by the models later on:
# log-price, log-returns, plain numbers (instead of xtc)

# Further it establishes 4 facts with evidence:
# 1. Prices trend; returns sit around zero (price and return panels)
# 2. Returns are close to uncorrelated once you use the right standard errors: robust Ljung-Box(10) p = 0.35; only lag 1 (-0.08, robust t = -2.6) survives -- hence the ARMA(1,1)/AR(1) mean in 05/07
# 3. Volatility clusters and persists (	ACF of r2 and |r|, ARCH-LM, yearly volatility)
# 4. Fat tails (kurtosis, Jarque–Bera, Q-Q plot, 3-sd count)

# We also work with robust variance, because:
# normal variance assumes the size of today's move has nothing to do with the size of the move k days earlier
# JPM's returns are nearly uncorrelated (a small lag 1 only) but not independent, because big days sit next to big days
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

need("01_data", c("jpm"))

part_begin("03_exploration")

## ---- 3. Data exploration --------------------------------------------------------

say("3. DATA EXPLORATION")

P  <- Ad(jpm); colnames(P) <- "adjclose"            # adjusted price, xts (keeps the dates)
r  <- na.omit(diff(log(P))); colnames(r) <- "ret"   # log return, xts
y  <- as.numeric(P)                                 # y = price, plain
ly <- log(y)                                        # log(y) = log price, plain
rn <- as.numeric(r)                                 # r = log return, plain

# Log returns: they add up over time (monthly RV in section 9 needs that), are
# symmetric, and a random walk in the log price keeps prices positive.

# The purpose is to show the two facts
# 1. Price wanders with no fixed level, while return stays around 0
# 2. The direction of returns is barely predictable (lag 1 only), but size is
# P, log(P), r, abs(r) 
fig("03_price_return_panels", plot_panels(
  list(P, log(P), r, abs(r)),
  titles = c(paste(PRIMARY, "adjusted close (levels)"),
             paste(PRIMARY, "log adjusted close"),
             "Daily log returns -- note the volatility clustering",
             "|log returns| -- persistence is in the SECOND moment"),
  main = paste0(PRIMARY, ", ", format(start(P)), " to ", format(end(P)))),
  height = 1250)

# plot -> shade -> redraw the line on top (rect() paints over it)
# Log(p) through the crisis 
fig("03_logprice_episodes", {
  par(mar = c(3, 4, 2.5, 1))
  plot(zoo::as.zoo(log(P)), type = "n", xlab = "", ylab = "log adjusted close",
       main = paste(PRIMARY, "log price with documented crisis episodes shaded"))
  shade_episodes(EPISODES)
  lines(zoo::as.zoo(log(P)), lwd = 1.6, col = "grey15")
  grid(col = "grey88")
  legend("topleft", bty = "n", cex = 0.85,
         legend = paste(names(EPISODES),
                        sapply(EPISODES, function(w) paste(substr(w, 1, 7), collapse = " to "))))
})

## --- y, log(y) and r side by side
# Same statistics for all three, so the contrast is in one table:
# y and log(y) have huge autocorrelation at lag 1 AND lag 10 (no fixed level,
# a unit root -- tested in 04); r has almost none. Skewness/kurtosis/Jarque-Bera
# of a trending level describe the path taken, not a distribution -- they are
# only meaningful for r.
sl <- list(y = y, "log(y)" = ly, r = rn)
desc <- do.call(cbind, lapply(sl, function(x) {
  jbx <- tseries::jarque.bera.test(x)
  lbx <- Box.test(x, lag = 10, type = "Ljung-Box")
  ac  <- as.numeric(acf(x, lag.max = 10, plot = FALSE)$acf)
  c(n = length(x), mean = mean(x), median = median(x), sd = sd(x),
    min = min(x), max = max(x),
    skewness = moments::skewness(x), excess_kurtosis = moments::kurtosis(x) - 3,
    jarque_bera = unname(jbx$statistic), jarque_bera_p = jbx$p.value,
    acf_lag1 = ac[2], acf_lag10 = ac[11],
    ljung_box10_iid = unname(lbx$statistic), ljung_box10_iid_p = lbx$p.value)
}))
desc <- data.frame(statistic = rownames(desc), desc, row.names = NULL, check.names = FALSE)
tab("03_descriptive_stats", desc, digits = 4,
    caption = paste("Descriptive statistics of y (price), log(y) and r (log return);",
                    "Ljung-Box = textbook iid version, the robust one is in 03_ljungbox_robust"))

# Histogram (top) over a horizontal boxplot (bottom) on the same x-axis, per
# series. y is right-skewed (the 2020s run-up sits in the upper tail); log(y)
# pulls that in; r is centred on 0, symmetric and peaked, with many points
# beyond the whiskers -- the fat tails.
fig("03_distribution_series", {
  layout(matrix(1:6, 2, 3), heights = c(3, 1.4))
  for (nm in names(sl)) {
    x <- sl[[nm]]; xl <- range(x)
    par(mar = c(1, 4, 3, 1))
    hist(x, breaks = 80, col = "grey70", border = NA, xlim = xl, xaxt = "n",
         main = paste("Distribution of", nm), xlab = "", ylab = "count")
    par(mar = c(4, 4, 0.5, 1))
    boxplot(x, horizontal = TRUE, ylim = xl, col = "grey85", outcex = 0.4,
            outcol = "grey40", xlab = nm)
  }
}, height = 900)

# One box per calendar year. y and log(y): the box moves, so the mean isn't
# constant (non-stationary). r: every box sits on 0, but its height changes
# with the regime -- constant mean, changing CONDITIONAL variance (volatility
# clustering; the unconditional variance can still be constant).
# 1999 (from May) and 2026 (to September) are partial years.
yr_y <- format(index(P), "%Y")
yr_r <- format(index(r), "%Y")
fig("03_boxplot_by_year", {
  par(mfrow = c(3, 1), mar = c(3.5, 4.5, 3, 1))
  boxplot(y ~ yr_y, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "y",
          main = "y (price) by year: the level drifts")
  boxplot(ly ~ yr_y, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "log(y)",
          main = "log(y) by year: the level still drifts")
  boxplot(rn ~ yr_r, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "r",
          main = "r by year: centred on 0, spread changes with the regime")
  abline(h = 0, col = "red", lty = 2)
}, height = 1400)

# Rolling 250-day (~1 year) mean and sd: a visual stationarity check.
# Stationary = the rolling stats keep coming back to one level; a unit root
# lets the mean wander off. y fails, and its sd grows with the price level --
# why we take logs. log(y) still has a wandering mean. r: no trend in the mean,
# it stays within ~2.5 se (sd(r)/sqrt(250) = 0.0015, red dotted = 2 se) of 0 --
# noise, not a moving level. Its sd swings ~8x with the crises: volatility
# clustering (conditional variance, what GARCH models), not a unit root.
ROLL <- 250
sx <- list(y = P, "log(y)" = log(P), r = r)
fig("03_rolling_stats", {
  par(mfrow = c(3, 2), mar = c(2.5, 4.5, 3, 1))
  for (nm in names(sx)) {
    x <- sx[[nm]]
    plot(zoo::as.zoo(zoo::rollapply(x, ROLL, mean, align = "right")), xlab = "",
         ylab = nm, main = paste0("Rolling ", ROLL, "-day mean of ", nm))
    if (nm == "r") {
      abline(h = 0, col = "red", lty = 2)
      abline(h = c(-2, 2) * sd(rn) / sqrt(ROLL), col = "red", lty = 3)
    }
    grid(col = "grey88")
    plot(zoo::as.zoo(zoo::rollapply(x, ROLL, sd, align = "right")), xlab = "",
         ylab = nm, main = paste0("Rolling ", ROLL, "-day sd of ", nm))
    grid(col = "grey88")
  }
}, height = 1300)

# Stylized facts as a table, one row per fact.
# stylized facts table provides statistics used for the rest of the case's 
# modelling choices (stylized-facts = patterns that daily stock-return series show)
# Values:
# 1. Sample and scale
# mean (daily, %)    = how small the avg return is next to the noise (t-stat shows its not significantly different from 0; Newey-West t = 1.5-1.6, still < 1.96)
# sd (daily, %)      = size of the daily moves
# sd (annualised, %) = same number the unit volatility is quoted on
# min (%) / max (%)  = the worst and best days 10 sd (should come every 10^20 years)

# 2. Shape of the distribution:
# skewness                     = checks whether big up days and big down days are equally common
# excess kurtosis              = measures how heavy the tails are
# Jarque–Bera stat             = the formal normality test, combining the two rows above 
# share of days |r| > 3 sd (%) = the kurtosis result in plain language
# normal benchmark (%)         = the share of days beyond 3 sd under a normal distribution

# 3. can past returns predict returns?:
# Ljung-Box(10) on r

# 4. Can past volatility predict volatility?
# Ljung–Box(10) on r2 and Ljung–Box(22) on |r|
# ARCH-LM(12) = Engle's formal test for ARCH effects, the test usually run before fitting a GARCH. (ARCH-LM is the regression form that GARCH itself is built on)
jb  <- tseries::jarque.bera.test(rn)
lb  <- function(x, k) Box.test(x, lag = k, type = "Ljung-Box")
arch <- FinTS::ArchTest(rn, lags = 12, demean = TRUE)

facts <- data.frame(
  statistic = c("n", "mean (daily, %)", "t-stat of the mean (iid)", "sd (daily, %)", "sd (annualised, %)",
                "min (%)", "max (%)", "skewness", "excess kurtosis",
                "Jarque-Bera stat", "Jarque-Bera p",
                "Ljung-Box(10) on r, stat", "Ljung-Box(10) on r, p",
                "Ljung-Box(10) on r^2, stat", "Ljung-Box(10) on r^2, p",
                "Ljung-Box(22) on |r|, stat", "Ljung-Box(22) on |r|, p",
                "ARCH-LM(12) stat", "ARCH-LM(12) p",
                "share of days |r| > 3 sd (%)", "normal benchmark (%)"),
  value = c(length(rn), mean(rn) * 100, mean(rn) / (sd(rn) / sqrt(length(rn))),
            sd(rn) * 100, sd(rn) * sqrt(252) * 100,
            min(rn) * 100, max(rn) * 100,
            moments::skewness(rn), moments::kurtosis(rn) - 3,
            unname(jb$statistic), jb$p.value,
            unname(lb(rn, 10)$statistic), lb(rn, 10)$p.value,
            unname(lb(rn^2, 10)$statistic), lb(rn^2, 10)$p.value,
            unname(lb(abs(rn), 22)$statistic), lb(abs(rn), 22)$p.value,
            unname(arch$statistic), arch$p.value,
            mean(abs(rn) > 3 * sd(rn)) * 100, (2 * pnorm(-3)) * 100))
tab("03_stylized_facts", facts, digits = 4,
    caption = "Stylized facts of JPM daily log returns")


# acf()'s dashed band is 1.96/sqrt(n), i.e. iid. Returns are nearly uncorrelated but
# heteroskedastic, so that band is too narrow and flags fake significant lags.
# Red = ARCH-robust band (Diebold 1986; Lobato-Nankervis-Savin 2001).
# se_robust = 1/sqrt(n)  ×  √( mean(a^2 · b^2) / mean(x^2)^2 )
# 2 weaknesses:
# 1. luck adds up over number of tests
# 2. small effects add up
ra_r <- robust_acf_se(rn, K = 40)
tab("03_robust_acf", head(ra_r, 12), digits = 5,
    caption = "Return ACF with iid vs heteroskedasticity-robust standard errors (first 12 lags)")

# the correction matters for the ACF of r; for r^2 the dependence IS the signal
# Outside the red band: lags 1 (-0.080), 33 and 34 -- 3 of 40, ~2 expected by
# chance, so only lag 1 counts. The PACF gets the same robust band.
fig("03_acf_contrast", {
  par(mfrow = c(2, 2), mar = ACF_MAR)   # ACF_MAR: room for acf()'s titles
  acf_robust(rn, lag.max = 40, main = "ACF of returns (red = ARCH-robust band)")
  acf(rn^2,    lag.max = 40, main = "ACF of SQUARED returns")
  acf(abs(rn), lag.max = 40, main = "ACF of |returns|")
  pacf_robust(rn, lag.max = 40, main = "PACF of returns (red = ARCH-robust band)")
}, height = 1050)

# same correction for the Ljung-Box test itself
# Ljung-Box: For each lag, measure how many standard errors the autocorrelation is from zero. That is its t-ratio. then squares each t-ratio
# se used in LB is still too narrow
# LBr:
# mean_se_inflation is the average correction factor over the lags, which tells you how much the correction mattered for this series.
lbr <- rbind(cbind(series = "r",   lb_robust(rn, 10)),
             cbind(series = "r",   lb_robust(rn, 22)),
             cbind(series = "r^2", lb_robust(rn^2, 10)),
             cbind(series = "|r|", lb_robust(abs(rn), 22)))
tab("03_ljungbox_robust", lbr, digits = 4,
    caption = "Ljung-Box: textbook versus heteroskedasticity-robust")

fig("03_distribution", {
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
  qqnorm(rn, pch = 16, cex = 0.35, col = "grey35",
         main = "Normal Q-Q: the tails bend away"); qqline(rn, col = "red", lwd = 2)
  h <- hist(rn, breaks = 200, plot = FALSE)
  plot(h$mids, h$density, type = "h", col = "grey60", xlab = "return", ylab = "density",
       main = "Return density vs fitted normal", xlim = c(-0.15, 0.15))
  curve(dnorm(x, mean(rn), sd(rn)), add = TRUE, col = "red", lwd = 2)
  legend("topright", c("empirical", "normal"), lty = 1, lwd = c(1, 2),
         col = c("grey60", "red"), bty = "n")
})

# Calendar effects: does the average return depend on the month or the
# weekday? Boxes = median and IQR, red dot = mean. Outliers aren't drawn, or the
# 2008-09 days squash every box flat.
# Test: regress r on month (weekday) dummies, H0 = all means equal. With the
# iid F, and with White (HC1) standard errors, which don't let the crisis
# days' volatility pass for a calendar effect.
mon <- factor(month.abb[as.integer(format(index(r), "%m"))], levels = month.abb)
wkd <- factor(c("Mon", "Tue", "Wed", "Thu", "Fri")[as.integer(format(index(r), "%u"))],
              levels = c("Mon", "Tue", "Wed", "Thu", "Fri"))
fig("03_boxplot_calendar", {
  par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
  for (g in list(list(f = mon, lab = "month"), list(f = wkd, lab = "weekday"))) {
    boxplot(rn * 100 ~ g$f, outline = FALSE, col = "grey85", xlab = "",
            ylab = "daily log return (%)", las = 2,
            main = paste("r by", g$lab, "(outliers not drawn)"))
    points(seq_along(levels(g$f)), tapply(rn * 100, g$f, mean), pch = 19, col = "red")
    abline(h = 0, lty = 2, col = "grey40")
  }
})
cal_test <- do.call(rbind, lapply(list(month = mon, weekday = wkd), function(f) {
  m  <- lm(rn ~ f)
  wi <- lmtest::waldtest(m)                                             # iid F
  wr <- lmtest::waldtest(m, vcov = sandwich::vcovHC(m, type = "HC1"))   # White
  data.frame(groups = nlevels(f), df = abs(wi$Df[2]),
             F_iid = wi$F[2], p_iid = wi$`Pr(>F)`[2],
             F_robust = wr$F[2], p_robust = wr$`Pr(>F)`[2])
}))
cal_test <- data.frame(effect = rownames(cal_test), cal_test, row.names = NULL)
tab("03_calendar_effects", cal_test, digits = 4,
    caption = "Calendar effects in r: H0 = equal mean return in every month / weekday")

# annual vol: clustering and the crises in one picture
yr    <- format(index(r), "%Y")
annvol <- tapply(rn, yr, function(z) sd(z) * sqrt(252) * 100)
fig("03_annual_volatility", {
  par(mar = c(4, 4, 2.5, 1))
  barplot(annvol, las = 2, col = "grey65", border = NA,
          ylab = "annualised sd of daily returns (%)",
          main = paste(PRIMARY, "realised volatility by calendar year"))
  abline(h = mean(annvol), lty = 2, col = "red")
})
tab("03_annual_volatility", data.frame(year = names(annvol),
                                       ann_vol_pct = round(as.numeric(annvol), 2)),
    caption = "Annualised realised volatility by year (%)")

part_save(c("arch", "lbr", "ly", "P", "r", "ra_r", "rn", "y"))
