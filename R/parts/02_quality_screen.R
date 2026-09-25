###############################################################################
# EBC4008 Case 3 -- part 02: Data-quality screen
#
# Inputs come from need(),
# outputs go to results/parts/02_quality_screen.rds via part_save().
#
# Needs : jpm, prices <- 01_data
# Makes : adjfac, lr_adj, lr_raw, qual, series_screen, series_summary
# Writes: tables 02_quality_screen, 02_largest_ranges_JPM,
#         02_range_ratio_subsamples, 02_return_screen, 02_series_summary,
#         02_series_screen
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

need("01_data", c("jpm", "prices"))

part_begin("02_quality_screen")

## ---- 2. Data-quality screen -----------------------------------------------------

say("2. DATA-QUALITY SCREEN")

# The Parkinson range (section 8) is quadratic in log(H/L): one bad tick does a
# lot of damage and flows into HAR and the spillover VAR. Yahoo cleans nothing.

qual <- do.call(rbind, lapply(prices[TICKERS], screen_quality)) # do.call => stacks the 6 rows into 1 table
tab("02_quality_screen", qual, digits = 3,
    caption = "Integrity counts and range diagnostics, six-bank panel")
# What is measured by screen_quality:
# 1. a value is missing (no range)
# 2. High = Low         (log0 breaks) 
# 3. High < Low         (data error)
# 4. Close or Open outside [Low, High] (off-exchange trades mixed in)
# 5. no shares traded   (stale day)
# 6. a price <= 0       (data error)

# return is log(new price / old price) log used because they add up:
# log(C_t / C_{t-1})  =  log(O_t / C_{t-1})  +  log(C_t / O_t)

# log(H/L) => range within a day, measure of volatility

# Parkinson estimator for range 
# RR = (log(H/L))^2 / (4 * log(2) realized range (range used to est variance)
# variance measure behind HAR + spillover analysis

# Geometric Brownian motion (GBM)
# Textbook model of a stockprice, random + continuous + constant volatility + no jumps
# check the data against the model's predictions

# 1. range_ratio_oc = mean(log(H / L)) / mean(abs(log(C / O)))
# Numerator: the average daily range. Denominator: the average size of the open-to-close move.
# ratio is always >= 1  under GBM the benchmark is 2

# 2. range_ratio_cc = mean(log(H / L)) / mean(abs(r)) 
# denom is bigger so ratio must be below 2

# 3. pred_garch_over_RR:  var(cc) / var(intraday)
# GARCH is fitted close-to-close (24h), RR is H-L (trading hours)
# GARCH should be higher (delta is explained by the overnight gap) this says how much

# Parkinson's formula
# 4. RR_over_intraday_var: mean( (log(H/L))^2 / (4*log(2)) )  /  mean( log(C/O)^2 )
# Numerator: the Parkinson estimate from the range.
# Denominator: the simplest estimate there is, today's open-to-close return squared.
# tells you RR is sound on its own.


assert_quality(qual)  # Checks are binded

# Largest ranges should be known event dates. sort() on an xts drops the dates,
# and x[order(...), ] gets re-sorted by time -- top_by_value() avoids both.
# sort() on an xts object throws the dates away, so use custom
lr <- log(Hi(jpm) / Lo(jpm)); colnames(lr) <- "log_HL"
tab("02_largest_ranges_JPM", top_by_value(lr, n = 12, value_name = "log_HL"),
    digits = 4,
    caption = "Twelve largest daily log(H/L) for JPM, RANKED -- should be known event dates")

# Same diagnostics per regime. If the range ratio moves between crisis and calm,
# log(H/L) is also tracking a microstructure shift, and part of section 12's
# crisis effect would be the measure.
# Returns are built on the FULL series and split afterwards. Differencing inside
# each subsample would create fake returns across every crisis/calm seam.
sub_ratio <- do.call(rbind, lapply(TICKERS, function(tk) {
  x <- prices[[tk]]
  H <- as.numeric(Hi(x)); L <- as.numeric(Lo(x))
  C <- as.numeric(Cl(x)); O <- as.numeric(Op(x))
  n <- length(C)
  rcc <- log(C[-1] / C[-n]); onr <- log(O[-1] / C[-n])   # close-to-close and overnight returns
  fl  <- Reduce(`|`, lapply(EPISODES, function(w)
    index(x) >= as.Date(w[1]) & index(x) <= as.Date(w[2]))) # crisis flag | => combined as ORs
  ratio <- function(k) mean(log(H[k] / L[k])) / mean(abs(log(C[k] / O[k])))
  onsh  <- function(k) var(onr[k]) / var(rcc[k])          # exact overnight share
  data.frame(ticker = tk, ratio_oc_crisis = ratio(fl), ratio_oc_calm = ratio(!fl),
             ratio_diff = ratio(fl) - ratio(!fl),
             overnight_crisis = onsh(fl[-1]), overnight_calm = onsh(!fl[-1]),
             row.names = NULL)
}))
tab("02_range_ratio_subsamples", sub_ratio, digits = 3,
    caption = "Interval-matched range ratio and overnight share, crisis vs calm")

# Return-level screen: stale prices (exact zeros, runs of them) and absurd jumps.
ret_screen <- do.call(rbind, lapply(TICKERS, function(tk) {
  rr <- diff(log(as.numeric(Ad(prices[[tk]])))); rr <- rr[is.finite(rr)]  # daily log returns (log price - day before)
  data.frame(ticker = tk, n = length(rr),
             exact_zero = sum(rr == 0),
             abs_gt_25pct = sum(abs(rr) > 0.25),   # count days moving more than 25%
             abs_gt_50pct = sum(abs(rr) > 0.50),   # count days moving more than 50%
             max_abs_pct = round(max(abs(rr)) * 100, 2),  # max movement in a day % 
             longest_zero_run = {
               rl <- rle(rr == 0); if (any(rl$values)) max(rl$lengths[rl$values]) else 0L
             }, row.names = NULL) # run-length encoding to count the number of 0s
}))
tab("02_return_screen", ret_screen,
    caption = "Return-level screen: stale prices and implausible jumps")

# The three series the case works with, for JPM:
#   y = adjusted close (price), log(y), r = diff(log(y)) (log return)
# 03 builds the same three for the models; here they only get screened.
yx  <- Ad(jpm)
ser <- list(y = yx, "log(y)" = log(yx), r = diff(log(yx))[-1])   # [-1] drops only the first-day NA, so a real gap would still show

# What each series looks like: size, range, where the extremes fall
series_summary <- do.call(rbind, lapply(names(ser), function(nm) {
  x <- ser[[nm]]; v <- as.numeric(x); d <- format(index(x))
  data.frame(series = nm, n = length(v), missing = sum(is.na(v)),
             first = d[1], last = d[length(d)],
             mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
             min = min(v, na.rm = TRUE), date_min = d[which.min(v)],
             max = max(v, na.rm = TRUE), date_max = d[which.max(v)],
             row.names = NULL)
}))
tab("02_series_summary", series_summary, digits = 4,
    caption = "y (price), log(y) and r (log return): size, range and extremes")

# Same stale/jump screen on each series, on its own day-to-day change dx:
#   y      -> dx = price change in $ (stale price = dx of exactly 0)
#   log(y) -> dx = r (so this row repeats the return screen above for JPM)
#   r      -> dx = change in the return; a spike that reverses the next day
#             (bad tick) shows up here as a huge |dx|. The largest one
#             (2009-01-20 -23%, 2009-01-21 +22%) is real: both days are in the
#             top-12 log(H/L) table, and a bad close print wouldn't widen the
#             intraday range.
# changes_gt_5sd counts days with |dx| more than 5 sd from its mean. Not for y:
# $ changes scale with the price ($8 -> $365), so a pooled 5-sd count only finds
# 2020-26 days -- jumps are screened on log(y) and r.
# A price <= 0 is an error; a log price or a return below 0 isn't, and a log
# price of exactly 0 ($1) means nothing (NA there).
series_screen <- do.call(rbind, lapply(names(ser), function(nm) {
  v  <- as.numeric(ser[[nm]])
  dx <- diff(v)
  dd <- format(index(ser[[nm]]))[-1]
  rl <- rle(dx == 0)
  data.frame(series = nm, n_changes = length(dx),
             nonpositive_price = if (nm == "y") sum(v <= 0) else NA_integer_,
             exact_zero_values = if (nm == "log(y)") NA_integer_ else sum(v == 0),
             zero_changes = sum(dx == 0),
             longest_zero_run = if (any(rl$values)) max(rl$lengths[rl$values]) else 0L,
             max_abs_change = max(abs(dx)), date_max_change = dd[which.max(abs(dx))],
             changes_gt_5sd = if (nm == "y") NA_integer_ else sum(abs(dx - mean(dx)) > 5 * sd(dx)),
             row.names = NULL)
}))
tab("02_series_screen", series_screen, digits = 4,
    caption = "Stale/jump screen on the daily change of y, log(y) and r")

# log(H/L) doesn't care about split/dividend adjustment: H and L get the same
# factor. Check it instead of assuming it.
adjfac <- as.numeric(Ad(jpm)) / as.numeric(Cl(jpm))
lr_raw <- log(as.numeric(Hi(jpm)) / as.numeric(Lo(jpm)))
lr_adj <- log((as.numeric(Hi(jpm)) * adjfac) / (as.numeric(Lo(jpm)) * adjfac))
stopifnot(max(abs(lr_raw - lr_adj)) < 1e-12)

part_save(c("adjfac", "lr_adj", "lr_raw", "qual", "series_screen", "series_summary"))
