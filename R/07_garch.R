###############################################################################
# EBC4008 Case 3 -- part 07: GARCH on r
# Reads : results/rds/01_data.rds (EPISODES)
#         results/rds/03_exploration.rds (r, rn)
# Writes: tables 07_garch_ic_grid, 07_persistence; figures 07_news_impact,
#         07_conditional_volatility, 07_std_residuals, 07_garch_forecast;
#         results/rds/07_garch.rds
###############################################################################

library(xts)       # xts(), index(); also attaches zoo
library(rugarch)   # ugarchspec(), ugarchfit(), ugarchforecast() + their diagnostics

d01 <- readRDS("results/rds/01_data.rds")
EPISODES <- d01$EPISODES
d03 <- readRDS("results/rds/03_exploration.rds")
r <- d03$r; rn <- d03$rn

# sGARCH(1,1), ML: sigma_t^2 = omega + alpha eps_{t-1}^2 + beta sigma_{t-1}^2

## --- Model selection grid
# 4 models x 4 dists -> combos; ugarchspec() GARCH(1,1), AR(1) mean -> mk_spec(m, d)
combos <- expand.grid(model = c("sGARCH", "eGARCH", "gjrGARCH", "apARCH"),
                      dist  = c("norm", "std", "sstd", "jsu"),
                      stringsAsFactors = FALSE)
mk_spec <- function(m, d) ugarchspec(
  variance.model     = list(model = m, garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = d)

set.seed(1)   # solver "hybrid" can fall back to the random-start gosolnp

# ugarchfit() all 16 (ML, solver "hybrid"); failed -> NA row, same names (else mapply -> list)
# -> ic_grid: 16 x (4 ICs per obs | npar | llh)
ic_grid <- t(mapply(function(m, d) {
  f <- tryCatch(ugarchfit(mk_spec(m, d), r, solver = "hybrid"),
                error = function(e) NULL)
  if (is.null(f) || convergence(f) != 0) {
    c(Akaike = NA_real_, Bayes = NA_real_, Shibata = NA_real_,
      `Hannan-Quinn` = NA_real_, npar = NA_real_, llh = NA_real_)
  } else {
    c(infocriteria(f)[, 1], npar = length(coef(f)), llh = likelihood(f))
  }
}, combos$model, combos$dist))

# rank by BIC; BIC_gap = total-BIC points behind best -> ic_tab, 07_garch_ic_grid
ic_tab <- cbind(combos, as.data.frame(ic_grid))
ic_tab <- ic_tab[order(ic_tab$Bayes), ]
ic_tab$BIC_total <- ic_tab$Bayes * nrow(r)
ic_tab$BIC_gap   <- (ic_tab$Bayes - min(ic_tab$Bayes, na.rm = TRUE)) * nrow(r)
write.csv(ic_tab, "results/tables/07_garch_ic_grid.csv")
print(ic_tab, digits = 5)

best_row <- ic_tab[1, ]   # BIC winner

# eGARCH-std by hand (BIC ~ ties apARCH-std; no positivity constraints); GJR-sstd = benchmark
SEL_MODEL <- "eGARCH"; SEL_DIST <- "std"
spec_main <- mk_spec(SEL_MODEL, SEL_DIST)
spec_gjr  <- mk_spec("gjrGARCH", "sstd")

# fit selected + benchmark -> fit_g, fit_gjr
fit_g   <- ugarchfit(spec_main, data = r, solver = "hybrid")
fit_gjr <- ugarchfit(spec_gjr,  data = r, solver = "hybrid")

stopifnot(convergence(fit_g) == 0)

## --- Diagnostics
# symmetric sGARCH-std -> fit_sym; Engle-Ng signbias() of 3 fits -> sign_bias (t | p tables)
fit_sym   <- ugarchfit(mk_spec("sGARCH", SEL_DIST), data = r, solver = "hybrid")
sign_bias <- list(selected = signbias(fit_g), gjr_sstd = signbias(fit_gjr),
                  sgarch = signbias(fit_sym))

# Nyblom-Hansen stability (H0: constant params) -> nyb
nyb <- nyblom(fit_g)

# persistence(), halflife() in days, 3 fits -> 07_persistence
pers <- data.frame(
  model = c(paste0(SEL_MODEL, "-", SEL_DIST), "gjrGARCH-sstd", paste0("sGARCH-", SEL_DIST)),
  persistence = c(persistence(fit_g), persistence(fit_gjr), persistence(fit_sym)),
  halflife_days = c(halflife(fit_g), halflife(fit_gjr), halflife(fit_sym)))
write.csv(pers, "results/tables/07_persistence.csv")
print(pers, digits = 4)

# newsimpact(): sigma_t^2 vs eps_{t-1}, 3 models -> 07_news_impact
png("results/figures/07_news_impact.png", width = 1500, height = 950, res = 130)
par(mar = c(4.2, 4.4, 2.5, 1))
nsel <- newsimpact(fit_g); ngjr <- newsimpact(fit_gjr); nsym <- newsimpact(fit_sym)
yl <- range(c(nsel$zy, ngjr$zy, nsym$zy))
plot(nsel$zx, nsel$zy, type = "l", lwd = 2.4, col = "firebrick", ylim = yl,
     xlab = expression(epsilon[t - 1]), ylab = expression(sigma[t]^2),
     main = "News impact curves: asymmetric models versus symmetric sGARCH")
lines(ngjr$zx, ngjr$zy, lwd = 2.2, col = "darkorange", lty = 4)
lines(nsym$zx, nsym$zy, lwd = 2.4, col = "steelblue",  lty = 2)
abline(v = 0, lty = 3); grid(col = "grey88")
legend("top", c(paste0(SEL_MODEL, " (selected)"), "gjrGARCH", "sGARCH"),
       col = c("firebrick", "darkorange", "steelblue"),
       lwd = 2.3, lty = c(1, 4, 2), bty = "n")
dev.off()

## --- Conditional variance, residuals
# in-sample sigma_t^2 of selected fit -> sigma_garch
sigma_garch <- sigma(fit_g)^2
colnames(sigma_garch) <- "garch_var"
stopifnot(inherits(sigma_garch, "xts"), nrow(sigma_garch) == nrow(r))

# annualised vol + episodes; red = sample sd, not uncvariance() -> 07_conditional_volatility
png("results/figures/07_conditional_volatility.png", width = 1500, height = 950, res = 130)
par(mar = c(3, 4.2, 2.5, 1))
av <- sqrt(252 * sigma_garch) * 100
plot(zoo::as.zoo(av), type = "n", xlab = "", ylab = "annualised sigma (%)",
     main = paste0(SEL_MODEL, "-", SEL_DIST,
                   " conditional volatility, annualised, with crisis episodes"))
usr <- par("usr")   # shade first: rect() paints over lines
for (w in EPISODES) rect(as.Date(w[1]), usr[3], as.Date(w[2]), usr[4],
                         col = adjustcolor("grey60", 0.28), border = NA)
lines(zoo::as.zoo(av), lwd = 1.1, col = "grey15")
abline(h = sd(rn) * sqrt(252) * 100, col = "red", lty = 2)
grid(col = "grey88")
legend("topright", c("conditional", "unconditional sd"), lwd = c(1.1, 1),
       lty = c(1, 2), col = c("grey15", "red"), bty = "n")
dev.off()

# z = eps / sigma: series, ACF z, ACF z^2, normal Q-Q -> 07_std_residuals
png("results/figures/07_std_residuals.png", width = 1500, height = 1050, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 3.8, 1))   # room for acf() titles
z <- as.numeric(residuals(fit_g, standardize = TRUE))
plot(zoo::as.zoo(xts(z, index(r))), xlab = "", ylab = "z",
     main = "Standardized residuals")
acf(z,   lag.max = 40, main = "ACF of z  (mean equation OK?)")
acf(z^2, lag.max = 40, main = "ACF of z^2  (variance equation OK?)")
qqnorm(z, pch = 16, cex = 0.35, col = "grey35", main = "Normal Q-Q of z")
qqline(z, col = "red", lwd = 2)
dev.off()

## --- Forecast
# ugarchforecast() 22 days; uncvariance() line above y-range (legend only) -> 07_garch_forecast
png("results/figures/07_garch_forecast.png", width = 1500, height = 950, res = 130)
par(mar = c(4, 4.2, 2.5, 1))
fcg <- ugarchforecast(fit_g, n.ahead = 22)
s   <- as.numeric(sigma(fcg)) * sqrt(252) * 100
plot(1:22, s, type = "b", pch = 16, lwd = 2, col = "firebrick",
     xlab = "days ahead", ylab = "annualised sigma (%)",
     main = "GARCH variance forecast reverts to the unconditional level")
abline(h = sqrt(uncvariance(fit_g) * 252) * 100, lty = 2)
grid(col = "grey88")
legend("topright", c("forecast", "unconditional"), lty = c(1, 2),
       col = c("firebrick", "black"), lwd = c(2, 1), bty = "n")
dev.off()

saveRDS(mget(c("best_row", "fit_g", "fit_gjr", "fit_sym", "ic_tab", "nyb",
               "SEL_DIST", "SEL_MODEL", "sigma_garch", "sign_bias", "spec_gjr",
               "spec_main")), "results/rds/07_garch.rds")
