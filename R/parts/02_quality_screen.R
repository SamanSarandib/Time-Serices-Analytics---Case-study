###############################################################################
# EBC4008 Case 3 -- part 02: data-quality screen
# Needs : jpm, prices <- 01_data
# Makes : adjfac, lr_adj, lr_raw, qual, series_screen, series_summary
# Writes: tables 02_quality_screen, 02_largest_ranges_JPM,
#         02_range_ratio_subsamples, 02_return_screen, 02_series_summary,
#         02_series_screen
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

need("01_data", c("jpm", "prices"))

part_begin("02_quality_screen")

# screen_quality() per bank: OHLC integrity counts, range ratios, variance shares -> 02_quality_screen, rows = tickers:
#   n | na_any | na_HL | zero_range | H_lt_L | C_outside_HL | O_outside_HL | zero_volume | nonpos_price | max_logHL | max_absret |
#   range_ratio_cc | range_ratio_oc | intraday_share | overnight_share | cov_share | pred_garch_over_RR | RR_over_intraday_var
qual <- do.call(rbind, lapply(prices[TICKERS], screen_quality))
tab("02_quality_screen", qual, digits = 3,
    caption = "Integrity counts and range diagnostics, six-bank panel")

# top_by_value(): the 12 largest JPM log(H/L) days, with dates
# -> 02_largest_ranges_JPM: rank | date | log_HL
lr <- log(Hi(jpm) / Lo(jpm)); colnames(lr) <- "log_HL"
tab("02_largest_ranges_JPM", top_by_value(lr, n = 12, value_name = "log_HL"),
    digits = 4,
    caption = "Twelve largest daily log(H/L) for JPM, RANKED -- should be known event dates")

# range ratio mean(log(H/L)) / mean(|log(C/O)|) and overnight share var(overnight) / var(c2c), crisis (EPISODES) vs calm
# returns built on the full series, then split by fl[-1] (no fake returns across crisis/calm seams)
# -> 02_range_ratio_subsamples: ticker | ratio_oc_crisis | ratio_oc_calm | ratio_diff | overnight_crisis | overnight_calm
sub_ratio <- do.call(rbind, lapply(TICKERS, function(tk) {
  x <- prices[[tk]]
  H <- as.numeric(Hi(x)); L <- as.numeric(Lo(x))
  C <- as.numeric(Cl(x)); O <- as.numeric(Op(x))
  n <- length(C)
  rcc <- log(C[-1] / C[-n]); onr <- log(O[-1] / C[-n])   # close-to-close and overnight returns
  fl  <- Reduce(`|`, lapply(EPISODES, function(w)
    index(x) >= as.Date(w[1]) & index(x) <= as.Date(w[2]))) # TRUE = inside any crisis window
  ratio <- function(k) mean(log(H[k] / L[k])) / mean(abs(log(C[k] / O[k]))) # How large is the average high-low range relative to the average open-to-close movement?
  onsh  <- function(k) var(onr[k]) / var(rcc[k]) 
  data.frame(ticker = tk, ratio_oc_crisis = ratio(fl), ratio_oc_calm = ratio(!fl),
             ratio_diff = ratio(fl) - ratio(!fl),
             overnight_crisis = onsh(fl[-1]), overnight_calm = onsh(!fl[-1]),
             row.names = NULL)
}))
tab("02_range_ratio_subsamples", sub_ratio, digits = 3,
    caption = "Interval-matched range ratio and overnight share, crisis vs calm")

# stale/jump screen on daily log returns of Ad(); zero runs via rle()
# -> 02_return_screen: ticker | n | exact_zero | abs_gt_25pct | abs_gt_50pct | max_abs_pct | longest_zero_run
ret_screen <- do.call(rbind, lapply(TICKERS, function(tk) {
  rr <- diff(log(as.numeric(Ad(prices[[tk]])))); rr <- rr[is.finite(rr)]
  data.frame(ticker = tk, n = length(rr),
             exact_zero = sum(rr == 0),
             abs_gt_25pct = sum(abs(rr) > 0.25),
             abs_gt_50pct = sum(abs(rr) > 0.50),
             max_abs_pct = round(max(abs(rr)) * 100, 2),
             longest_zero_run = {
               rl <- rle(rr == 0); if (any(rl$values)) max(rl$lengths[rl$values]) else 0L
             }, row.names = NULL)
}))
tab("02_return_screen", ret_screen,
    caption = "Return-level screen: stale prices and implausible jumps")

# JPM series: y = Ad() price, log(y), r = diff(log(y)); list of 3 xts
yx  <- Ad(jpm)
ser <- list(y = yx, "log(y)" = log(yx), r = diff(log(yx))[-1])   # [-1]: drops the first-day NA only

# summary stats per series in ser, with the dates of min and max
# -> 02_series_summary: series | n | missing | first | last | mean | sd | min | date_min | max | date_max
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

# stale/jump screen on each series' daily change dx (rle() zero runs, |dx - mean| > 5 sd); NA = check doesn't apply
# -> 02_series_screen: series | n_changes | nonpositive_price | exact_zero_values | zero_changes |
#    longest_zero_run | max_abs_change | date_max_change | changes_gt_5sd
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

part_save(c("qual", "series_screen", "series_summary"))