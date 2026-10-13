###############################################################################
# EBC4008 Case 3 -- part 03: data exploration
# Reads : results/rds/01_data.rds (jpm, EPISODES)
# Writes: tables 03_descriptive_stats, 03_stylized_facts, 03_robust_acf, 03_ljungbox_robust,
#         03_calendar_effects, 03_annual_volatility; figures 03_price_return_panels,
#         03_logprice_episodes, 03_distribution_series, 03_boxplot_by_year, 03_rolling_stats,
#         03_acf_contrast, 03_distribution, 03_boxplot_calendar, 03_annual_volatility;
#         results/rds/03_exploration.rds
###############################################################################

library(quantmod)   # Ad(); also attaches xts and zoo

d01 <- readRDS("results/rds/01_data.rds")
jpm <- d01$jpm; EPISODES <- d01$EPISODES

# k series in k stacked panels (plot.zoo), bold outer title
plot_panels <- function(series, titles, main, ...) {
  par(mfrow = c(length(series), 1), mar = c(2.4, 4, 2, 1), oma = c(0, 0, 2.2, 0))
  for (i in seq_along(series)) {
    plot(zoo::as.zoo(series[[i]]), xlab = "", ylab = "", col = "grey20", main = titles[i], ...)
    grid(col = "grey88")
  }
  mtext(main, outer = TRUE, cex = 1.05, font = 2)
}

# ACF lags 1..K + ARCH-robust se (Diebold 1986) -> data.frame; inflation = se_robust / se_iid
robust_acf_se <- function(x, K) {
  x <- x - mean(x)
  n <- length(x)
  d <- mean(x^2)^2
  se <- sapply(1:K, function(k) sqrt(mean(x[(k + 1):n]^2 * x[1:(n - k)]^2) / d / n))
  data.frame(lag = 1:K, acf = as.numeric(acf(x, lag.max = K, plot = FALSE)$acf)[-1],
             se_iid = 1 / sqrt(n), se_robust = se, inflation = se * sqrt(n))
}

# Ljung-Box Q(lags), textbook + ARCH-robust, both chi2(lags) -> 1-row data.frame
lb_robust <- function(x, lags) {
  n  <- length(x)
  ra <- robust_acf_se(x, lags)
  q_std <- n * (n + 2) * sum(ra$acf^2 / (n - ra$lag))
  q_rob <- sum((ra$acf / ra$se_robust)^2)
  data.frame(lags = lags,
             LB_standard = q_std, p_standard = pchisq(q_std, lags, lower.tail = FALSE),
             LB_robust   = q_rob, p_robust   = pchisq(q_rob, lags, lower.tail = FALSE),
             mean_se_inflation = mean(ra$inflation))
}

P  <- Ad(jpm); colnames(P) <- "adjclose"            # adjusted close, xts
r  <- na.omit(diff(log(P))); colnames(r) <- "ret"   # log return, xts
y  <- as.numeric(P)                                 # price, numeric
ly <- log(y)                                        # log price, numeric
rn <- as.numeric(r)                                 # log return, numeric

titles <- c("JPM adjusted close (levels)", "JPM log adjusted close",
            "Daily log returns", "|log returns|")

# P, log(P), r, |r| -> 03_price_return_panels
png("results/figures/03_price_return_panels.png", width = 1500, height = 1250, res = 130)
plot_panels(list(P, log(P), r, abs(r)), titles,
            main = paste0("JPM, ", format(start(P)), " to ", format(end(P))))
dev.off()

# P, log(P), r, EPISODES shaded -> 03_logprice_episodes
png("results/figures/03_logprice_episodes.png", width = 1500, height = 1100, res = 130)
plot_panels(list(P, log(P), r), titles[1:3],
            main = "JPM with documented crisis episodes shaded",
            panel = function(x, y, ...) {
              usr <- par("usr")                        # bands first: behind the series
              for (w in EPISODES)
                rect(as.Date(w[1]), usr[3], as.Date(w[2]), usr[4],
                     col = adjustcolor("grey60", 0.28), border = NA)
              lines(x, y, ...)
              if (par("mfg")[1] == 1)                  # mfg[1] = panel row; legend on top
                legend("topleft", bty = "n", cex = 0.85,
                       legend = paste(names(EPISODES), sapply(EPISODES, function(w)
                         paste(substr(w, 1, 7), collapse = " to "))))
            })
dev.off()

# summary stats of y, log(y), r -> 03_descriptive_stats
# (Ljung-Box = textbook iid version, the robust one is in 03_ljungbox_robust)
sl <- list(y = y, "log(y)" = ly, r = rn)
desc <- do.call(cbind, lapply(sl, function(x) {
  jbx <- tseries::jarque.bera.test(x)              # H0: normal
  lbx <- Box.test(x, lag = 10, type = "Ljung-Box") # H0: iid
  ac  <- as.numeric(acf(x, lag.max = 10, plot = FALSE)$acf)
  c(n = length(x), mean = mean(x), median = median(x), sd = sd(x),
    min = min(x), max = max(x),
    skewness = moments::skewness(x), excess_kurtosis = moments::kurtosis(x) - 3,
    jarque_bera = unname(jbx$statistic), jarque_bera_p = jbx$p.value,
    acf_lag1 = ac[2], acf_lag10 = ac[11],
    ljung_box10_iid = unname(lbx$statistic), ljung_box10_iid_p = lbx$p.value)
}))
desc <- data.frame(statistic = rownames(desc), desc, row.names = NULL, check.names = FALSE)
write.csv(desc, "results/tables/03_descriptive_stats.csv")
print(desc, digits = 4)

# histogram over boxplot per series -> 03_distribution_series
png("results/figures/03_distribution_series.png", width = 1500, height = 900, res = 130)
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
dev.off()

# boxplots by year (first/last partial) -> 03_boxplot_by_year
yr_y <- format(index(P), "%Y")
yr_r <- format(index(r), "%Y")
png("results/figures/03_boxplot_by_year.png", width = 1500, height = 1400, res = 130)
par(mfrow = c(3, 1), mar = c(3.5, 4.5, 3, 1))
boxplot(y ~ yr_y, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "y",
        main = "y (price) by year: the level drifts")
boxplot(ly ~ yr_y, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "log(y)",
        main = "log(y) by year: the level still drifts")
boxplot(rn ~ yr_r, col = "grey85", outcex = 0.4, las = 2, xlab = "", ylab = "r",
        main = "r by year: centred on 0, spread changes with the regime")
abline(h = 0, col = "red", lty = 2)
dev.off()

# rolling mean, sd -> 03_rolling_stats
ROLL <- 250
sx <- list(y = P, "log(y)" = log(P), r = r)
png("results/figures/03_rolling_stats.png", width = 1500, height = 1300, res = 130)
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
dev.off()

# stylized facts of JPM daily log returns -> 03_stylized_facts
jb    <- tseries::jarque.bera.test(rn)
lb_r  <- Box.test(rn,      lag = 10, type = "Ljung-Box")
lb_r2 <- Box.test(rn^2,    lag = 10, type = "Ljung-Box")
lb_a  <- Box.test(abs(rn), lag = 22, type = "Ljung-Box")
arch  <- FinTS::ArchTest(rn, lags = 12, demean = TRUE)

facts <- data.frame(
  statistic = c("n", "mean (daily, %)", "t-stat of the mean (iid)", "sd (daily, %)",
                "sd (annualised, %)", "min (%)", "max (%)", "skewness", "excess kurtosis",
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
            unname(lb_r$statistic),  lb_r$p.value,
            unname(lb_r2$statistic), lb_r2$p.value,
            unname(lb_a$statistic),  lb_a$p.value,
            unname(arch$statistic), arch$p.value,
            mean(abs(rn) > 3 * sd(rn)) * 100, (2 * pnorm(-3)) * 100))
write.csv(facts, "results/tables/03_stylized_facts.csv")
print(facts, digits = 4)

# return ACF, iid vs ARCH-robust se (first 12 lags) -> ra_r, 03_robust_acf
ra_r <- robust_acf_se(rn, K = 40)
write.csv(head(ra_r, 12), "results/tables/03_robust_acf.csv")
print(head(ra_r, 12), digits = 5)

# ACF/PACF; ARCH-robust band (red, from ra_r) for r only -> 03_acf_contrast
png("results/figures/03_acf_contrast.png", width = 1500, height = 1050, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3.8, 1))   # room for acf() titles
a <- acf(rn, lag.max = 40, plot = FALSE)
a$acf <- a$acf[-1, , , drop = FALSE]; a$lag <- a$lag[-1, , , drop = FALSE]   # drop lag 0
plot(a, main = "ACF of returns (red = ARCH-robust band)",
     ylim = range(a$acf, 2 * ra_r$se_robust, -2 * ra_r$se_robust))
lines(ra_r$lag,  1.96 * ra_r$se_robust, col = "red", lwd = 1.6)
lines(ra_r$lag, -1.96 * ra_r$se_robust, col = "red", lwd = 1.6)
acf(rn^2,    lag.max = 40, main = "ACF of SQUARED returns")
acf(abs(rn), lag.max = 40, main = "ACF of |returns|")
p <- pacf(rn, lag.max = 40, plot = FALSE)
plot(p, main = "PACF of returns (red = ARCH-robust band)",
     ylim = range(p$acf, 2 * ra_r$se_robust, -2 * ra_r$se_robust))
lines(ra_r$lag,  1.96 * ra_r$se_robust, col = "red", lwd = 1.6)
lines(ra_r$lag, -1.96 * ra_r$se_robust, col = "red", lwd = 1.6)
dev.off()

# Ljung-Box, textbook vs ARCH-robust -> lbr, 03_ljungbox_robust
lbr <- rbind(cbind(series = "r",   lb_robust(rn, 10)),
             cbind(series = "r",   lb_robust(rn, 22)),
             cbind(series = "r^2", lb_robust(rn^2, 10)),
             cbind(series = "|r|", lb_robust(abs(rn), 22)))
write.csv(lbr, "results/tables/03_ljungbox_robust.csv")
print(lbr, digits = 4)

# normal Q-Q of r; density vs fitted normal -> 03_distribution
png("results/figures/03_distribution.png", width = 1500, height = 950, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
qqnorm(rn, pch = 16, cex = 0.35, col = "grey35",
       main = "Normal Q-Q: the tails bend away"); qqline(rn, col = "red", lwd = 2)
h <- hist(rn, breaks = 200, plot = FALSE)
plot(h$mids, h$density, type = "h", col = "grey60", xlab = "return", ylab = "density",
     main = "Return density vs fitted normal", xlim = c(-0.15, 0.15))
curve(dnorm(x, mean(rn), sd(rn)), add = TRUE, col = "red", lwd = 2)
legend("topright", c("empirical", "normal"), lty = 1, lwd = c(1, 2),
       col = c("grey60", "red"), bty = "n")
dev.off()

# r (%) by month, weekday; red dot = mean -> 03_boxplot_calendar
mon <- factor(month.abb[as.integer(format(index(r), "%m"))], levels = month.abb)
wkd <- factor(c("Mon", "Tue", "Wed", "Thu", "Fri")[as.integer(format(index(r), "%u"))],
              levels = c("Mon", "Tue", "Wed", "Thu", "Fri"))
cal <- list(month = mon, weekday = wkd)
png("results/figures/03_boxplot_calendar.png", width = 1500, height = 950, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))
for (nm in names(cal)) {
  f <- cal[[nm]]
  boxplot(rn * 100 ~ f, outline = FALSE, col = "grey85", xlab = "",
          ylab = "daily log return (%)", las = 2,
          main = paste("r by", nm, "(outliers not drawn)"))
  points(seq_along(levels(f)), tapply(rn * 100, f, mean), pch = 19, col = "red")
  abline(h = 0, lty = 2, col = "grey40")
}
dev.off()

# joint F on month / weekday dummies, iid vs HC1; H0 = equal mean return -> 03_calendar_effects
cal_test <- do.call(rbind, lapply(cal, function(f) {
  m  <- lm(rn ~ f)
  wi <- lmtest::waldtest(m)                                             # iid F
  wr <- lmtest::waldtest(m, vcov = sandwich::vcovHC(m, type = "HC1"))   # White
  data.frame(groups = nlevels(f), df = abs(wi$Df[2]),
             F_iid = wi$F[2], p_iid = wi$`Pr(>F)`[2],
             F_robust = wr$F[2], p_robust = wr$`Pr(>F)`[2])
}))
cal_test <- data.frame(effect = rownames(cal_test), cal_test, row.names = NULL)
write.csv(cal_test, "results/tables/03_calendar_effects.csv")
print(cal_test, digits = 4)

# annualised sd of r (%) by year -> 03_annual_volatility (figure + table)
annvol <- tapply(rn, yr_r, function(z) sd(z) * sqrt(252) * 100)
png("results/figures/03_annual_volatility.png", width = 1500, height = 950, res = 130)
par(mar = c(4, 4, 2.5, 1))
barplot(annvol, las = 2, col = "grey65", border = NA,
        ylab = "annualised sd of daily returns (%)",
        main = "JPM realised volatility by calendar year")
abline(h = mean(annvol), lty = 2, col = "red")
dev.off()
annvol_tab <- data.frame(year = names(annvol), ann_vol_pct = round(as.numeric(annvol), 2))
write.csv(annvol_tab, "results/tables/03_annual_volatility.csv")
print(annvol_tab, digits = 4)

saveRDS(mget(c("arch", "lbr", "ly", "P", "r", "ra_r", "rn", "y")),
        "results/rds/03_exploration.rds")
