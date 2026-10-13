#Modeling volatility: GARCH models 

# We take the returns (r=xts, rn=without dates) from 03_exploration

#   (i)  introduce GARCH: which stylized facts do they capture, and how do
#        they differ conceptually from the ARIMA models?
#   (ii) find an appropriate R package, estimate appropriate (G)ARCH models for 
#        your return series, present and interpret them.
#        We used the rugarch package
# Package: rugarch (Ghalanos). ugarchspec() + ugarchfit() estimate by ML;
# sigma() extracts the fitted conditional sd, residuals(standardize = TRUE)
# the standardized residuals z_t = (r_t - mu_t) / sigma_t


library(rugarch)

d01 <- readRDS("results/rds/01_data.rds")
EPISODES <- d01$EPISODES
d03 <- readRDS("results/rds/03_exploration.rds")
arch <- d03$arch; r <- d03$r; rn <- d03$rn

PRIMARY <- "JPM"

say <- function(...) cat(..., "\n", sep = "")

tab <- function(name, x, digits = NULL, caption = NULL) {
  write.csv(x, paste0("results/tables/", name, ".csv"))
  cat(caption, "\n")
  print(x, digits = digits)
}

fig <- function(name, expr, height = 950) {
  png(paste0("results/figures/", name, ".png"), width = 1500, height = height, res = 130)
  on.exit(dev.off())
  expr
}

cache <- function(key, expr) expr

stars <- function(p) ifelse(p < .01, "***", ifelse(p < .05, "**", ifelse(p < .10, "*", "")))

shade_episodes <- function(episodes) {
  usr <- par("usr")
  for (w in episodes) rect(as.Date(w[1]), usr[3], as.Date(w[2]), usr[4],
                           col = adjustcolor("grey60", 0.28), border = NA)
}

robust_acf_se <- function(x, K) {
  x <- x - mean(x)
  n <- length(x)
  d <- mean(x^2)^2
  se <- sapply(1:K, function(k) sqrt(mean(x[(k + 1):n]^2 * x[1:(n - k)]^2) / d / n))
  data.frame(lag = 1:K, acf = as.numeric(acf(x, lag.max = K, plot = FALSE)$acf)[-1],
             se_iid = 1 / sqrt(n), se_robust = se, inflation = se * sqrt(n))
}

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

##Modeling volatility: GARCH models

say("7. GARCH MODELS FOR THE RETURNS")

##MEAN equation 
#in script 05 we show that the conditional mean of returns is close to zero and hard to predict
#Consequently, we fix a constant mean (armaOrder = c(0, 0)) for every GARCH spec. below:
#the object of interest here is the VARIANCE equation.
#We fit on r (the xts object, not rn) needed for the episode-shaded volatility plot below.

MEAN_SPEC <- list(armaOrder = c(0, 0), include.mean = TRUE)

## --- Step 1: order search (iterative), sGARCH-std as the workhorse spec
# Rather than assuming (1,1) outright, we iterate over a small grid of
# candidate orders -- (1,1), (1,2), (2,1), (2,2) -- on the plain sGARCH with
# Student-t errors, and let BIC pick among them (same logic as the AR/MA
# order search done by hand for r in 05, just applied to the VARIANCE
# equation here instead of the mean). We search one model family (sGARCH) to
# keep the grid small: (1,1) is already the near-universal finding in the
# GARCH literature for daily equity returns, so this is a check that our
# data do not call for something richer, not an open-ended search.
ORDER_GRID <- list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))

order_fits <- lapply(ORDER_GRID, function(o) {
  tag <- paste0(o[1], o[2])
  cache(paste0("garch_order_", tag), {
    ugarchfit(spec = ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = o),
      mean.model = MEAN_SPEC, distribution.model = "std"),
      data = r, solver = "hybrid")
  })
})
names(order_fits) <- sapply(ORDER_GRID, function(o) paste0("sGARCH(", o[1], ",", o[2], ")_std"))

garch_order_cmp <- do.call(rbind, lapply(names(order_fits), function(k) {
  f  <- order_fits[[k]]
  ic <- infocriteria(f)
  data.frame(spec = k, loglik = likelihood(f), AIC = ic["Akaike", ], BIC = ic["Bayes", ],
             row.names = NULL)
}))
tab("07_garch_order_cmp", garch_order_cmp, digits = 5,
    caption = "sGARCH-std order search: (p,q) combinations compared by AIC/BIC")

best_order_key <- garch_order_cmp$spec[which.min(garch_order_cmp$BIC)]
best_order     <- ORDER_GRID[[which.min(garch_order_cmp$BIC)]]
say("Order search (sGARCH-std) selected by BIC: ", best_order_key)

#VARIANCE EQUATION
##1.Standard GARCH model
# sGARCH-norm: Bollerslev (1986) GARCH(1,1), Gaussian errors, but in script03 theJB test already rejects normality
# sGARCH-std: Student-t errors-->lets the model account for the fat tails 

##Asymmetric extensions of GARCH models take into account the well-known LEVERAGE EFFECT
#the leverage effect basically means that bad news tends to raise volatility more than good news of the same size
##2. Exponential GARCH 
# eGARCH-std: adds an asymmetric response of sigma_t to the SIGN of the
#               past shock. 

# gjrGARCH-std: an alternative, additive way of modelling the same asymmetry
#               (an indicator for negative shocks), kept alongside eGARCH as
#               a robustness check on how the leverage term is parameterised.

#Model specification with ugarchspec
specs <- list(
  sGARCH_norm  = ugarchspec(variance.model = list(model = "sGARCH",   garchOrder = best_order),
                            mean.model = MEAN_SPEC, distribution.model = "norm"),
  sGARCH_std   = ugarchspec(variance.model = list(model = "sGARCH",   garchOrder = best_order),
                            mean.model = MEAN_SPEC, distribution.model = "std"),
  eGARCH_std   = ugarchspec(variance.model = list(model = "eGARCH",   garchOrder = best_order),
                            mean.model = MEAN_SPEC, distribution.model = "std"),
  gjrGARCH_std = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = best_order),
                            mean.model = MEAN_SPEC, distribution.model = "std"))

#Model estimation with ugarchfit
garch_fits <- lapply(names(specs), function(k) {
  cache(paste0("garch_", k), {
    ugarchfit(spec = specs[[k]], data = r, solver = "hybrid")
  })
})
names(garch_fits) <- names(specs)

##We compare the four models using information criterions using infocriteria() 

#persistence = how long a volatility shock takes to die out. 
#It is computed in different ways according to the model we are estimating:
#sGARCH / gjrGARCH: alpha1 + beta1 (+1/2 gamma for gjr's asymmetric term)
#eGARCH:            beta1  
#Values close to 1 say shocks decay slowly
persistence <- function(fit, model) {
  cf <- coef(fit)
  switch(model,
         sGARCH   = unname(cf["alpha1"] + cf["beta1"]),
         eGARCH   = unname(cf["beta1"]),
         gjrGARCH = unname(cf["alpha1"] + cf["beta1"] + cf["gamma1"] / 2))
}
MODEL_OF <- c(sGARCH_norm = "sGARCH", sGARCH_std = "sGARCH",
              eGARCH_std = "eGARCH", gjrGARCH_std = "gjrGARCH")

#comparison table according to AIC, BIC, persistence
garch_cmp <- do.call(rbind, lapply(names(garch_fits), function(k) {
  f  <- garch_fits[[k]]
  ic <- infocriteria(f)
  data.frame(spec = k, model = MODEL_OF[[k]],
             dist = if (grepl("std", k)) "std" else "norm",
             loglik = likelihood(f), n = length(rn),
             AIC = ic["Akaike", ], BIC = ic["Bayes", ],
             persistence = persistence(f, MODEL_OF[[k]]),
             row.names = NULL)
}))
tab("07_garch_cmp", garch_cmp, digits = 5,
    caption = paste0("GARCH(", best_order[1], ",", best_order[2],
                     ") family on r: information criteria and shock persistence"))

#BIC picks the parsimony-favouring model, consistent with 05's choice of BIC
#for r (returns are close to unpredictable, so charge heavily per parameter)
best_key   <- garch_cmp$spec[which.min(garch_cmp$BIC)] #model with the lowest BIC
garch_best <- garch_fits[[best_key]]
say("Selected by BIC: ", best_key)

##We test the statistic significance of the coefficients for the selected model
mc <- garch_best@fit$matcoef
garch_coefs <- data.frame(term = rownames(mc), estimate = mc[, 1], se = mc[, 2],
                          t = mc[, 3], p = mc[, 4], row.names = NULL)
garch_coefs$sig <- stars(garch_coefs$p)
tab("07_garch_coefs", garch_coefs, digits = 5,
    caption = paste("Estimated coefficients,", best_key))

##Conditional volatility: plain, and against the crisis episodes
sigma_best <- sigma(garch_best)         
fig("07_garch_vol", {
  par(mar = c(3, 4, 2.5, 1))
  plot(zoo::as.zoo(sigma_best), col = "blue", lwd = 1.1, xlab = "",
       ylab = "conditional sd", main = paste(PRIMARY, best_key, "-- fitted conditional volatility"))
  grid(col = "lightgrey")
})

#We highlights crisis periods with grey areas 
#the GARCH volatility spikes line up with the crisis periods
fig("07_garch_vol_episodes", {
  par(mar = c(3, 4, 2.5, 1))
  plot(zoo::as.zoo(sigma_best), type = "n", xlab = "", ylab = "conditional sd",
       main = paste(PRIMARY, best_key, "conditional volatility with crisis episodes"))
  shade_episodes(EPISODES)
  lines(zoo::as.zoo(sigma_best), lwd = 1.3, col = "blue")
  grid(col = "lightgrey")
  legend("topleft", bty = "n", cex = 0.85,
         legend = paste(names(EPISODES),
                        sapply(EPISODES, function(w) paste(substr(w, 1, 7), collapse = " to "))))
})

##We use the squared returns as a proxy for the latent volatility
#True volatility is never observed directly. 
fig("07_garch_vol_vs_sqreturns", {
  par(mar = c(3, 4, 2.5, 1))
  r2 <- zoo::as.zoo(r)^2
  plot(r2, col = "darkgrey", type = "h", xlab = "", ylab = expression(r[t]^2~"  vs  "~sigma[t]^2),
       main = paste(PRIMARY, best_key, "r_t^2 vs fitted conditional variance"))
  lines(zoo::as.zoo(sigma_best)^2, col = "blue", lwd = 1.3)
  grid(col = "lightgrey")
  legend("topleft", bty = "n", cex = 0.85, lty = 1, lwd = c(1, 1.3),
         col = c("darkgrey", "blue"),
         legend = c("r_t^2 (proxy)", paste0("sigma_t^2 (fitted, ", best_key, ")")))
})

##DIAGNOSTICS
##Did the variance equation actually remove the ARCH effect?

#We compute the standardized residuals:  z_t=(r_t-mu_t)/sigma_t. 
#We expected not to find structure in the residuals, since if 
#sigma_t is well-specified, z_t should follow a WN
z <- as.numeric(residuals(garch_best, standardize = TRUE))

garch_std_resid_tests <- rbind(
  data.frame(test = "Ljung-Box(10) on z",   lb_robust(z,   10)[c("LB_robust", "p_robust")]),
  data.frame(test = "Ljung-Box(10) on z^2", lb_robust(z^2, 10)[c("LB_robust", "p_robust")]))

arch_z <- FinTS::ArchTest(z, lags = 12, demean = TRUE)
garch_std_resid_tests <- rbind(garch_std_resid_tests,
  data.frame(test = "ARCH-LM(12) on z", LB_robust = unname(arch_z$statistic), p_robust = arch_z$p.value))
tab("07_garch_std_resid", garch_std_resid_tests, digits = 4,
    caption = paste(best_key, "standardized residuals",
                    "(compare to arch/lbr from 03, run on r itself before GARCH)"))

fig("07_garch_std_resid", {
  par(mfrow = c(2, 2), mar = c(4, 4, 3.5, 1))
  acf(z,   lag.max = 40, main = "ACF of standardized residuals z")
  acf(z^2, lag.max = 40, main = "ACF of z^2")
  qqnorm(z, pch = 16, cex = 0.35, col = "grey35", main = "Normal Q-Q of z"); qqline(z, col = "red", lwd = 2)
  plot(density(z), col = "grey30", lwd = 1.5, main = "Density of z vs normal", xlab = "z")
  curve(dnorm(x), add = TRUE, col = "red", lwd = 2)
}, height = 1100)

##News Impact Curves (Engle & Ng, 1993)
#To show that the standard GARCH doesn't capture the leverage effect, indeed
#it gives us the classic "smile", which is symmetric
#Asymmetric specifications instead (eGARCH, gjr) are able to capture the 
#asymmetric effect between good and bad news of the same magnitude on the volatility

nic_keys <- c("sGARCH_std", "eGARCH_std", "gjrGARCH_std")
nic <- lapply(nic_keys, function(k) newsimpact(object = garch_fits[[k]]))
names(nic) <- nic_keys

fig("07_garch_nic", {
  par(mar = c(4, 4, 2.5, 1))
  xr <- range(unlist(lapply(nic, function(x) x$zx)))
  yr <- range(unlist(lapply(nic, function(x) x$zy)))
  cols <- c(sGARCH_std = "red", eGARCH_std = "pink", gjrGARCH_std = "lightgreen")
  plot(nic$sGARCH_std$zx, nic$sGARCH_std$zy, type = "l", lwd = 2, col = cols[1],
       xlim = xr, ylim = yr, xlab = "shock (past return)", ylab = "implied variance",
       main = paste(PRIMARY, "-- News Impact Curves"))
  for (k in nic_keys[-1]) lines(nic[[k]]$zx, nic[[k]]$zy, lwd = 2, col = cols[k])
  grid(col = "grey88")
  legend("top", bty = "n", lwd = 2, col = cols,
         legend = c("sGARCH-std (symmetric)", "eGARCH-std", "gjrGARCH-std"))
})

saveRDS(mget(c("best_order", "garch_best", "garch_cmp", "garch_fits", "garch_order_cmp",
               "garch_std_resid_tests", "nic", "sigma_best")), "results/rds/07_garch.rds")

