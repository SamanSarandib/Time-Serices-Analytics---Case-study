###############################################################################
# EBC4008 Case 3 -- part 03: Data exploration
# Needs : jpm <- 01_data
# Makes : P, r (xts adjusted close, log return); y, ly, rn (plain numeric price,
#         log price, log return); ra_r, lbr (data.frames); arch (ArchTest htest)
# Writes: tables 03_descriptive_stats, 03_stylized_facts, 03_robust_acf,
#         03_ljungbox_robust, 03_calendar_effects, 03_annual_volatility;
#         figures 03_price_return_panels, 03_logprice_episodes,
#         03_distribution_series, 03_boxplot_by_year, 03_rolling_stats,
#         03_acf_contrast, 03_distribution, 03_boxplot_calendar,
#         03_annual_volatility
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

need("01_data", c("jpm"))

part_begin("03_exploration")

P  <- Ad(jpm); colnames(P) <- "adjclose"            # adjusted price, xts (keeps the dates)
r  <- na.omit(diff(log(P))); colnames(r) <- "ret"   # log return, xts
y  <- as.numeric(P)                                 # y = price, plain
ly <- log(y)                                        # log(y) = log price, plain
rn <- as.numeric(r)                                 # r = log return, plain

# -> 03_price_return_panels: 4 rows P | log(P) | r | |r|
fig("03_price_return_panels", plot_panels(
  list(P, log(P), r, abs(r)),
  titles = c(paste(PRIMARY, "adjusted close (levels)"),
             paste(PRIMARY, "log adjusted close"),
             "Daily log returns",
             "|log returns|"),
  main = paste0(PRIMARY, ", ", format(start(P)), " to ", format(end(P)))),
  height = 1250)

# -> 03_logprice_episodes: 3 rows P | log(P) | r, legend in the top panel
fig("03_logprice_episodes", plot_panels(
  list(P, log(P), r),
  titles = c(paste(PRIMARY, "adjusted close (levels)"),
             paste(PRIMARY, "log adjusted close"),
             "Daily log returns"),
  main = paste(PRIMARY, "with documented crisis episodes shaded"),
  panel = function(x, y, ...) {
    shade_episodes(EPISODES)                  # first, so the bands sit behind the series
    lines(x, y, ...)
    if (par("mfg")[1] == 1)                   # mfg[1] = row of the panel being drawn: legend in the top one
      legend("topleft", bty = "n", cex = 0.85,
             legend = paste(names(EPISODES),
                            sapply(EPISODES, function(w) paste(substr(w, 1, 7), collapse = " to "))))
  }),
  height = 1100)

# Data summary
# -> 03_descriptive_stats: statistic | y | log(y) | r (14 rows)
sl <- list(y = y, "log(y)" = ly, r = rn)
desc <- do.call(cbind, lapply(sl, function(x) {
  jbx <- tseries::jarque.bera.test(x)              # H0: normal (p = prob of normal)
  lbx <- Box.test(x, lag = 10, type = "Ljung-Box") # H0: iid (p = prob of iid)
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

# -> 03_distribution_series: 2x3 grid, rows histogram | boxplot, columns y | log(y) | r
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

# one box per calendar year (1999 and 2026 are partial years)
# -> 03_boxplot_by_year: 3 rows y | log(y) | r (r with a red dashed 0 line)
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

# visual stationarity check
# -> 03_rolling_stats: 3x2 grid, rows y | log(y) | r, columns mean | sd
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

# stylized facts of r: mean (+ iid t), sd, min/max, skew, kurtosis, Jarque-Bera, iid Ljung-Box
# on r(10), r^2(10), |r|(22), ARCH-LM(12) (FinTS::ArchTest), share of |r| > 3 sd vs the normal's
# -> 03_stylized_facts: statistic | value (21 rows)
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


# ACF of r with iid vs ARCH-robust se (robust_acf_se)
# -> ra_r / 03_robust_acf: lag | acf | se_iid | se_robust | inflation (= se_robust / se_iid)
ra_r <- robust_acf_se(rn, K = 40)
tab("03_robust_acf", head(ra_r, 12), digits = 5,
    caption = "Return ACF with iid vs heteroskedasticity-robust standard errors (first 12 lags)")

# acf_robust()/pacf_robust() add the ARCH-robust band in red; r^2 and |r| keep plain acf()
# -> 03_acf_contrast: 2x2 grid, ACF r | ACF r^2 | ACF |r| | PACF r (40 lags)
fig("03_acf_contrast", {
  par(mfrow = c(2, 2), mar = ACF_MAR)   # ACF_MAR: room for acf()'s titles
  acf_robust(rn, lag.max = 40, main = "ACF of returns (red = ARCH-robust band)")
  acf(rn^2,    lag.max = 40, main = "ACF of SQUARED returns")
  acf(abs(rn), lag.max = 40, main = "ACF of |returns|")
  pacf_robust(rn, lag.max = 40, main = "PACF of returns (red = ARCH-robust band)")
}, height = 1050)

# Ljung-Box, textbook vs ARCH-robust (lb_robust), on r (10, 22 lags), r^2 (10), |r| (22)
# -> lbr / 03_ljungbox_robust: series | lags | LB_standard | p_standard | LB_robust | p_robust | mean_se_inflation
lbr <- rbind(cbind(series = "r",   lb_robust(rn, 10)),
             cbind(series = "r",   lb_robust(rn, 22)),
             cbind(series = "r^2", lb_robust(rn^2, 10)),
             cbind(series = "|r|", lb_robust(abs(rn), 22)))
tab("03_ljungbox_robust", lbr, digits = 4,
    caption = "Ljung-Box: textbook versus heteroskedasticity-robust")

# qqnorm() + qqline(); 200-bin histogram density vs dnorm(mean(r), sd(r)), x cut to +-0.15
# -> 03_distribution: 2 panels, normal Q-Q of r | density vs fitted normal
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

# boxplots of r (%) by month and weekday, red dot = mean, outliers not drawn
# -> 03_boxplot_calendar: 2 panels month | weekday
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
# lm(rn ~ f) for f = month, weekday; joint F test of the dummies (lmtest::waldtest), iid and White HC1 vcov
# -> 03_calendar_effects: effect | groups | df | F_iid | p_iid | F_robust | p_robust
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

# annualised sd of r per calendar year (sd * sqrt(252) * 100)
# -> figure 03_annual_volatility: one bar per year, red dashed = mean; table: year | ann_vol_pct
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
