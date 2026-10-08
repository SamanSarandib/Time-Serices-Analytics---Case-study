###############################################################################
# EBC4008 Case 3 -- part 02: data-quality screen
# Reads : results/rds/01_data.rds (prices, jpm, TICKERS, EPISODES)
# Writes: tables 02_quality_screen, 02_largest_ranges_JPM, 02_range_ratio_subsamples,
#         02_return_screen, 02_series_summary, 02_series_screen;
#         results/rds/02_quality_screen.rds
###############################################################################

library(quantmod)   # Op(), Hi(), Lo(), Cl(), Vo(), Ad(); also attaches xts and zoo

d01 <- readRDS("results/rds/01_data.rds")
prices <- d01$prices; jpm <- d01$jpm; TICKERS <- d01$TICKERS; EPISODES <- d01$EPISODES

# per bank: OHLC integrity counts, range ratios, variance shares -> qual, 02_quality_screen
qual <- do.call(rbind, lapply(prices[TICKERS], function(x) {
  H <- as.numeric(Hi(x)); L <- as.numeric(Lo(x))
  C <- as.numeric(Cl(x)); O <- as.numeric(Op(x))
  V <- as.numeric(Vo(x)); A <- as.numeric(Ad(x))
  n    <- length(C)
  r    <- diff(log(A))                  # adjusted log return
  lhl  <- log(H / L)                    # daily log range
  rcc  <- log(C[-1] / C[-n])            # RAW close-to-close = overnight + open-to-close
  onr  <- log(O[-1] / C[-n])            # overnight
  intr <- log(C[-1] / O[-1])            # open-to-close
  vt   <- var(rcc,  na.rm = TRUE)
  vi   <- var(intr, na.rm = TRUE)
  data.frame(
    n            = n,
    na_any       = sum(!complete.cases(cbind(O, H, L, C, V, A))),
    na_HL        = sum(is.na(H) | is.na(L)),
    zero_range   = sum(H == L, na.rm = TRUE),          # -> log(H/L) = 0
    H_lt_L       = sum(H < L, na.rm = TRUE),
    C_outside_HL = sum(C > H | C < L, na.rm = TRUE),   # after-hours prints
    O_outside_HL = sum(O > H | O < L, na.rm = TRUE),
    zero_volume  = sum(V == 0, na.rm = TRUE),
    nonpos_price = sum(c(O, H, L, C, A) <= 0, na.rm = TRUE),
    max_logHL    = round(max(lhl, na.rm = TRUE), 4),
    max_absret   = round(max(abs(r), na.rm = TRUE), 4),
    range_ratio_cc = round(mean(lhl, na.rm = TRUE) / mean(abs(r), na.rm = TRUE), 3),
    range_ratio_oc = round(mean(lhl, na.rm = TRUE) / mean(abs(log(C / O)), na.rm = TRUE), 3),
    # variance shares of raw close-to-close: intraday + overnight + 2cov = 1
    intraday_share  = round(vi / vt, 3),
    overnight_share = round(var(onr, na.rm = TRUE) / vt, 3),
    cov_share       = round(2 * cov(onr, intr, use = "complete.obs") / vt, 3),
    pred_garch_over_RR   = round(1 / (vi / vt), 3),
    RR_over_intraday_var = round(mean(lhl^2 / (4 * log(2)), na.rm = TRUE) /
                                 mean(log(C / O)^2, na.rm = TRUE), 3)
  )
}))
write.csv(qual, "results/tables/02_quality_screen.csv")
print(qual, digits = 3)

# 12 largest JPM log(H/L) days, ranked -- should be known event dates -> 02_largest_ranges_JPM
lr <- as.numeric(log(Hi(jpm) / Lo(jpm)))
o  <- order(lr, decreasing = TRUE)[1:12]
largest_ranges <- data.frame(rank = 1:12, date = format(index(jpm)[o]), log_HL = lr[o])
write.csv(largest_ranges, "results/tables/02_largest_ranges_JPM.csv")
print(largest_ranges, digits = 4)

# range ratio, overnight share, crisis vs calm; returns split by fl[-1]
# -> 02_range_ratio_subsamples
sub_ratio <- do.call(rbind, lapply(TICKERS, function(tk) {
  x <- prices[[tk]]
  H <- as.numeric(Hi(x)); L <- as.numeric(Lo(x))
  C <- as.numeric(Cl(x)); O <- as.numeric(Op(x))
  n <- length(C)
  rcc <- log(C[-1] / C[-n]); onr <- log(O[-1] / C[-n])   # close-close, overnight
  fl  <- Reduce(`|`, lapply(EPISODES, function(w)
    index(x) >= as.Date(w[1]) & index(x) <= as.Date(w[2]))) # TRUE = in a crisis window
  ratio <- function(k) mean(log(H[k] / L[k])) / mean(abs(log(C[k] / O[k])))
  onsh  <- function(k) var(onr[k]) / var(rcc[k])
  data.frame(ticker = tk, ratio_oc_crisis = ratio(fl), ratio_oc_calm = ratio(!fl),
             ratio_diff = ratio(fl) - ratio(!fl),
             overnight_crisis = onsh(fl[-1]), overnight_calm = onsh(!fl[-1]),
             row.names = NULL)
}))
write.csv(sub_ratio, "results/tables/02_range_ratio_subsamples.csv")
print(sub_ratio, digits = 3)

# stale prices and implausible jumps in daily Ad() log returns; zero runs via rle()
# -> 02_return_screen
ret_screen <- do.call(rbind, lapply(TICKERS, function(tk) {
  rr <- diff(log(as.numeric(Ad(prices[[tk]])))); rr <- rr[is.finite(rr)]
  rl <- rle(rr == 0)
  data.frame(ticker = tk, n = length(rr),
             exact_zero = sum(rr == 0),
             abs_gt_25pct = sum(abs(rr) > 0.25),
             abs_gt_50pct = sum(abs(rr) > 0.50),
             max_abs_pct = round(max(abs(rr)) * 100, 2),
             longest_zero_run = max(0L, rl$lengths[rl$values]),   # 0 if no zero run
             row.names = NULL)
}))
write.csv(ret_screen, "results/tables/02_return_screen.csv")
print(ret_screen, digits = 4)

# JPM y = Ad(), log(y), r = diff(log(y)) -> ser
y   <- Ad(jpm)
ser <- list(y = y, "log(y)" = log(y), r = diff(log(y))[-1])   # [-1]: drop first-day NA

# size, range and extremes (with dates) of each ser series -> series_summary, 02_series_summary
series_summary <- do.call(rbind, lapply(names(ser), function(nm) {
  x <- ser[[nm]]; v <- as.numeric(x); d <- format(index(x))
  data.frame(series = nm, n = length(v), missing = sum(is.na(v)),
             first = d[1], last = d[length(d)],
             mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
             min = min(v, na.rm = TRUE), date_min = d[which.min(v)],
             max = max(v, na.rm = TRUE), date_max = d[which.max(v)],
             row.names = NULL)
}))
write.csv(series_summary, "results/tables/02_series_summary.csv")
print(series_summary, digits = 4)

# stale/jump screen on the daily change dx of each ser series; NA = check n/a
# -> series_screen, 02_series_screen
series_screen <- do.call(rbind, lapply(names(ser), function(nm) {
  x  <- ser[[nm]]; v <- as.numeric(x)
  dx <- diff(v)
  dd <- format(index(x))[-1]
  rl <- rle(dx == 0)
  data.frame(series = nm, n_changes = length(dx),
             nonpositive_price = if (nm == "y") sum(v <= 0) else NA_integer_,
             exact_zero_values = if (nm == "log(y)") NA_integer_ else sum(v == 0),
             zero_changes = sum(dx == 0),
             longest_zero_run = max(0L, rl$lengths[rl$values]),
             max_abs_change = max(abs(dx)), date_max_change = dd[which.max(abs(dx))],
             changes_gt_5sd = if (nm == "y") NA_integer_ else
                                sum(abs(dx - mean(dx)) > 5 * sd(dx)),
             row.names = NULL)
}))
write.csv(series_screen, "results/tables/02_series_screen.csv")
print(series_screen, digits = 4)

saveRDS(mget(c("qual", "series_screen", "series_summary")), "results/rds/02_quality_screen.rds")
