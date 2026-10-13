###############################################################################
# EBC4008 Case 3 -- part 11: Diebold-Yilmaz volatility spillovers
# Reads : results/rds/01_data.rds (prices, TICKERS, EPISODES)
# Writes: tables 11_panel_stationarity, 11_var_lag_criteria, 11_var_residual_serial,
#         11_spillover_table, 11_directional_spillovers, 11_robustness_p_H,
#         11_p_sensitivity_by_H; figures 11_net_directional, 11_rolling_spillover,
#         11_rolling_net, 11_rolling_overall_pkg; results/rds/11_spillovers.rds
###############################################################################

library(quantmod)                 # Hi(), Lo(); also attaches xts and zoo
library(bootUR)                   # boot_sqt()
library(vars)                     # VARselect(), VAR(), serial.test()
library(frequencyConnectedness)   # spilloverDY12(), spilloverRollingDY12()

d01 <- readRDS("results/rds/01_data.rds")
prices <- d01$prices; TICKERS <- d01$TICKERS; EPISODES <- d01$EPISODES

# Diebold-Yilmaz (2012): generalized FEVD of VAR(p) on 6 banks' log RR, horizon H

## --- Panel
# Parkinson range (log H - log L)^2 / (4 log 2) per bank, merged, logged;
# n_before = rows before na.omit() -> panel
panel <- do.call(merge, lapply(prices[TICKERS], function(x) {
  rr <- (log(Hi(x)) - log(Lo(x)))^2 / (4 * log(2))
  rr[!is.finite(rr) | rr <= 0] <- NA       # H == L days: log(0) = -Inf would survive na.omit()
  rr
}))
colnames(panel) <- TICKERS
panel <- log(panel)
n_before <- nrow(panel); panel <- na.omit(panel)

# boot_sqt(): intercept, MAIC lags, AWB p; Romano-Wolf Step-M -> ur_panel, 11_panel_stationarity
maxlag <- trunc(12 * (nrow(panel) / 100)^(1 / 4))
set.seed(1)
ur_panel <- boot_sqt(panel, data_name = "log RR panel", steps = 0:ncol(panel),
                     bootstrap = "AWB", B = 1999, union = FALSE, deterministics = "intercept",
                     detrend = "OLS", max_lag = maxlag, criterion = "MAIC",
                     criterion_scale = FALSE, show_progress = FALSE)
panel_ur <- data.frame(series    = TICKERS,
                       lags_used = as.integer(ur_panel$details$`selected lags`[TICKERS, 1]),
                       tau       = ur_panel$details$`individual statistics`[TICKERS, 1],
                       p_boot    = ur_panel$details$`individual p-values`[TICKERS, 1],
                       reject_stepM = ur_panel$rejections[TICKERS, 1], row.names = NULL)
write.csv(panel_ur, "results/tables/11_panel_stationarity.csv")
print(panel_ur, digits = 3)
stopifnot(all(panel_ur$reject_stepM))

## --- VAR
# VARselect() to 20 lags, p = BIC (SC) -> sel; 11_var_lag_criteria: lag x AIC | HQ | SC | FPE
sel <- VARselect(panel, lag.max = 20, type = "const")
lag_crit <- t(round(sel$criteria, 4))
write.csv(as.data.frame(lag_crit), "results/tables/11_var_lag_criteria.csv")
print(lag_crit, digits = 4)
P_VAR <- unname(as.integer(sel$selection["SC(n)"]))

# VAR(P_VAR) + const, OLS -> est; max_root = largest companion-root modulus (< 1: stable)
est <- VAR(panel, p = P_VAR, type = "const")
max_root <- max(abs(vars::roots(est)))

# serial.test() PT(16), BG(5) -> pt, bg; PT(16) at p = 2, 4, ..., 10 -> 11_var_residual_serial
pt <- serial.test(est, lags.pt = 16, type = "PT.asymptotic")
bg <- serial.test(est, lags.bg = 5,  type = "BG")
serial_by_p <- do.call(rbind, lapply(c(2L, 4L, 6L, 8L, 10L), function(pp) {
  v <- VAR(panel, p = pp, type = "const")
  s <- serial.test(v, lags.pt = 16, type = "PT.asymptotic")
  data.frame(p = pp, PT16_stat = unname(s$serial$statistic), PT16_p = unname(s$serial$p.value))
}))
write.csv(serial_by_p, "results/tables/11_var_residual_serial.csv")
print(serial_by_p, digits = 4)

## --- Spillovers, full sample
# spilloverDY12(no.corr = FALSE) = generalized FEVD (vars::fevd: Cholesky only) -> sp
H_BASE <- 10L
sp <- spilloverDY12(est, n.ahead = H_BASE, no.corr = FALSE)

# DY table (%): FROM/TO = off-diag row/col sums, corner = TOTAL -> 11_spillover_table
S   <- 100 * sp$tables[[1]]
off <- S; diag(off) <- 0
spill_tab <- rbind(cbind(S, "FROM others" = rowSums(off)),
                   "TO others" = c(colSums(off), mean(rowSums(off))))
write.csv(as.data.frame(spill_tab), "results/tables/11_spillover_table.csv")
print(spill_tab, digits = 2)
stopifnot(abs(mean(rowSums(off)) - overall(sp)[[1]]) < 1e-8)

# to(), from(), net() (= DY table / N), sorted by net -> 11_directional_spillovers
dir_tab <- data.frame(series = TICKERS,
                      to   = as.numeric(to(sp)[[1]]),
                      from = as.numeric(from(sp)[[1]]),
                      net  = as.numeric(net(sp)[[1]]))
dir_tab <- dir_tab[order(-dir_tab$net), ]
write.csv(dir_tab, "results/tables/11_directional_spillovers.csv")
print(dir_tab, digits = 2)

# net per bank, red = transmitter, blue = receiver -> 11_net_directional
png("results/figures/11_net_directional.png", width = 1500, height = 950, res = 130)
par(mar = c(4, 4.5, 2.5, 1))
barplot(dir_tab$net, names.arg = dir_tab$series,
        col = ifelse(dir_tab$net > 0, "firebrick", "steelblue"), border = NA,
        ylab = "net spillover (percentage points)",
        main = paste0("Net directional volatility spillovers, full sample (H = ", H_BASE, ")"))
abline(h = 0)
dev.off()

## --- Robustness
# total index over p x H grid; [[1]] on overall(), not on c(...) (would return p)
# -> 11_robustness_p_H: p | H=4 | H=10 | ... | H=200
ph_rob <- expand.grid(p = c(2L, 4L, 6L, 8L), H = c(4L, 10L, 20L, 50L, 100L, 200L))
ph_rob$total <- mapply(function(pp, hh)
  overall(spilloverDY12(VAR(panel, p = pp, type = "const"),
                        n.ahead = hh, no.corr = FALSE))[[1]], ph_rob$p, ph_rob$H)
ph_wide <- reshape(ph_rob, idvar = "p", timevar = "H", direction = "wide")
names(ph_wide) <- sub("^total\\.", "H=", names(ph_wide))
write.csv(ph_wide, "results/tables/11_robustness_p_H.csv")
print(ph_wide, digits = 2)

# index range across p per H -> 11_p_sensitivity_by_H; h100_200 = max |H=200 - H=100|
band <- tapply(ph_rob$total, ph_rob$H, function(z) max(z) - min(z))
p_sens <- data.frame(H = as.integer(names(band)), p_range_pp = as.numeric(band))
write.csv(p_sens, "results/tables/11_p_sensitivity_by_H.csv")
print(p_sens, digits = 3)
h100_200 <- max(abs(ph_wide[["H=200"]] - ph_wide[["H=100"]]))

# order check: banks shuffled (seed 7) must give same index -> perm, sp_perm
set.seed(7)
perm <- sample(TICKERS)
sp_perm <- spilloverDY12(VAR(panel[, perm], p = P_VAR, type = "const"),
                         n.ahead = H_BASE, no.corr = FALSE)
stopifnot(abs(overall(sp)[[1]] - overall(sp_perm)[[1]]) < 1e-6)

## --- Rolling index
# VAR(P_VAR) + DY index on every 200-day window, on all cores but 2 -> roll
# do.call(): workers get values, not promises (else "object 'P_VAR' not found")
ROLL_WINDOW <- 200L
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 2))
parallel::clusterEvalQ(cl, library(frequencyConnectedness))
roll <- do.call(spilloverRollingDY12,
                list(data = panel, n.ahead = H_BASE, no.corr = FALSE, func_est = "VAR",
                     params_est = list(p = P_VAR, type = "const"), window = ROLL_WINDOW,
                     cluster = cl))
parallel::stopCluster(cl)

# overall(roll) is POSIXct-indexed; rect() needs Date -> tci
tci <- overall(roll)[[1]]
tci <- xts(coredata(tci), order.by = as.Date(index(tci)))
colnames(tci) <- "total"
stopifnot(inherits(index(tci), "Date"))

# rolling index, episodes shaded, red dashed = full-sample index -> 11_rolling_spillover
png("results/figures/11_rolling_spillover.png", width = 1500, height = 950, res = 130)
par(mar = c(3, 4.5, 2.5, 1))
plot(zoo::as.zoo(tci), type = "n", xlab = "",
     ylab = "total spillover index (%)",
     main = paste0("Rolling total volatility spillover index (", ROLL_WINDOW,
                   "-day window, VAR(", P_VAR, "), H = ", H_BASE, ")"))
usr <- par("usr")
for (w in EPISODES) rect(as.Date(w[1]), usr[3], as.Date(w[2]), usr[4],
                         col = adjustcolor("grey60", 0.28), border = NA)
lines(zoo::as.zoo(tci), lwd = 1.6, col = "grey10")
abline(h = overall(sp)[[1]], col = "red", lty = 2)
grid(col = "grey90"); box()      # rect() paints over frame
legend("bottomleft", c("rolling index", "full-sample index"), lwd = c(1.6, 1),
       lty = c(1, 2), col = c("grey10", "red"), bty = "n", cex = 0.85)
dev.off()
# episodes must lie inside tci, else band not drawn
stopifnot(all(sapply(EPISODES, function(w)
  as.Date(w[1]) >= start(tci) & as.Date(w[2]) <= end(tci))))

# package plots: rolling net, rolling total -> 11_rolling_net, 11_rolling_overall_pkg
png("results/figures/11_rolling_net.png", width = 1500, height = 1150, res = 130)
plotNet(roll)
dev.off()
png("results/figures/11_rolling_overall_pkg.png", width = 1500, height = 950, res = 130)
plotOverall(roll)
dev.off()

saveRDS(mget(c("bg", "H_BASE", "h100_200", "max_root", "n_before", "P_VAR", "panel",
               "perm", "pt", "ROLL_WINDOW", "sel", "sp", "sp_perm", "tci", "ur_panel")),
        "results/rds/11_spillovers.rds")
