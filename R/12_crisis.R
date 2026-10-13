###############################################################################
# Case 3 -- part 12: crisis vs non-crisis volatility spillovers
#
# Reuses the multivariate method from the supplied part 11 script:
# Parkinson realized range -> log range panel -> VAR selected by BIC ->
# 10-step generalized FEVD spillovers (Diebold-Yilmaz 2012).
#
# The crisis is the same GFC window used previously: 2007-08-01 to 2009-06-30.
# It is compared with the most recent equally long contiguous non-crisis window
# before the GFC. The full-sample estimate is included as context.
###############################################################################

library(xts)

d01 <- readRDS("results/rds/01_data.rds")
prices <- d01$prices; TICKERS <- d01$TICKERS; EPISODES <- d01$EPISODES

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

var_max_lag <- function(Tn, K, lag_max, min_obs_per_param) {
  p <- seq_len(lag_max)
  ok <- (Tn - p) >= min_obs_per_param * (K * p + 1)
  if (any(ok)) max(p[ok]) else 0L
}

say("12. CRISIS VS NON-CRISIS VOLATILITY SPILLOVERS")

FEVD_HORIZON <- 10L
VAR_LAG_MAX <- 20L
MIN_OBS_PER_PARAM <- 10L

## ---- Build the same log realized-range panel as part 11 ----------------------

range_by_bank <- lapply(TICKERS, function(ticker) {
  x <- prices[[ticker]]
  rr <- (log(quantmod::Hi(x)) - log(quantmod::Lo(x)))^2 / (4 * log(2))
  rr[!is.finite(rr) | rr <= 0] <- NA_real_
  colnames(rr) <- ticker
  rr
})

range_panel <- do.call(merge, range_by_bank)
colnames(range_panel) <- TICKERS
n_before <- nrow(range_panel)
log_range_panel <- na.omit(log(range_panel))
dates <- zoo::index(log_range_panel)
say("Balanced range panel: ", n_before, " rows before removing incomplete dates; ",
    nrow(log_range_panel), " complete trading dates remain.")

if (ncol(log_range_panel) != length(TICKERS))
  stop("The balanced realized-range panel does not contain all configured banks.")
if (nrow(log_range_panel) <= FEVD_HORIZON + length(TICKERS) * 10L)
  stop("Too few balanced observations for the spillover VAR.")

## ---- Define matched crisis and non-crisis periods ----------------------------

gfc_window <- EPISODES[["gfc"]]
gfc_start <- as.Date(gfc_window[1L])
gfc_end <- as.Date(gfc_window[2L])
gfc_idx <- which(dates >= gfc_start & dates <= gfc_end)
if (length(gfc_idx) < 60L)
  stop("Too few observations in the configured GFC window.")
if (any(diff(gfc_idx) != 1L))
  stop("The GFC observations are not a single contiguous block in the range panel.")

# Exclude every named crisis when looking for a clean comparison window.
in_any_episode <- rep(FALSE, length(dates))
for (window in EPISODES) {
  in_any_episode <- in_any_episode |
    (dates >= as.Date(window[1L]) & dates <= as.Date(window[2L]))
}
eligible_noncrisis <- dates < gfc_start & !in_any_episode

# Pick the latest contiguous clean block before the GFC, matching its length.
runs <- rle(eligible_noncrisis)
run_end <- cumsum(runs$lengths)
eligible_runs <- which(runs$values & runs$lengths >= length(gfc_idx))
if (!length(eligible_runs))
  stop("No contiguous pre-GFC non-crisis block is long enough for a matched comparison.")
chosen_run <- max(eligible_runs)
noncrisis_idx <- seq.int(run_end[chosen_run] - length(gfc_idx) + 1L,
                         run_end[chosen_run])

period_indices <- list(
  full_sample = seq_len(nrow(log_range_panel)),
  noncrisis_pre_gfc = noncrisis_idx,
  crisis_gfc = gfc_idx
)
period_labels <- c(full_sample = "Full sample",
                   noncrisis_pre_gfc = "Matched non-crisis",
                   crisis_gfc = "GFC crisis")

say("GFC window: ", format(min(dates[gfc_idx])), " to ", format(max(dates[gfc_idx])),
    " (", length(gfc_idx), " trading dates)")
say("Matched non-crisis window: ", format(min(dates[noncrisis_idx])), " to ",
    format(max(dates[noncrisis_idx])), " (", length(noncrisis_idx), " trading dates)")

## ---- Estimate VAR and generalized spillovers for each period -----------------

fit_spillovers <- function(label, idx) {
  Y <- as.matrix(log_range_panel[idx, , drop = FALSE])
  storage.mode(Y) <- "double"
  Tn <- nrow(Y)
  K <- ncol(Y)

  lag_cap <- var_max_lag(Tn, K, lag_max = VAR_LAG_MAX,
                         min_obs_per_param = MIN_OBS_PER_PARAM)
  if (lag_cap < 1L)
    stop("The ", label, " sample is too short for a stable VAR lag search.")

  selection <- vars::VARselect(Y, lag.max = lag_cap, type = "const")
  p <- unname(as.integer(selection$selection["SC(n)"]))
  if (length(p) != 1L || is.na(p) || p < 1L || p > lag_cap)
    stop("BIC lag selection failed for ", label, ".")

  fit <- vars::VAR(Y, p = p, type = "const")
  max_root <- max(abs(vars::roots(fit)))
  spill <- frequencyConnectedness::spilloverDY12(
    fit, n.ahead = FEVD_HORIZON, no.corr = FALSE)

  total <- as.numeric(frequencyConnectedness::overall(spill)[[1L]])
  to_values <- as.numeric(frequencyConnectedness::to(spill)[[1L]])
  from_values <- as.numeric(frequencyConnectedness::from(spill)[[1L]])
  net_values <- as.numeric(frequencyConnectedness::net(spill)[[1L]])

  if (any(lengths(list(to_values, from_values, net_values)) != K))
    stop("Unexpected directional spillover dimensions for ", label, ".")

  summary <- data.frame(
    period = label,
    start = as.Date(min(dates[idx])),
    end = as.Date(max(dates[idx])),
    n = Tn,
    VAR_lag_cap = lag_cap,
    VAR_lag_BIC = p,
    FEVD_horizon = FEVD_HORIZON,
    total_spillover_pct = total,
    max_root = max_root,
    VAR_stable = max_root < 1,
    row.names = NULL
  )

  directional <- data.frame(
    period = label,
    ticker = TICKERS,
    TO_others_pct = to_values,
    FROM_others_pct = from_values,
    net_TO_pct = net_values,
    row.names = NULL
  )

  lag_criteria <- data.frame(
    period = label,
    lag = seq_len(ncol(selection$criteria)),
    AIC = as.numeric(selection$criteria["AIC(n)", ]),
    HQ = as.numeric(selection$criteria["HQ(n)", ]),
    BIC = as.numeric(selection$criteria["SC(n)", ]),
    FPE = as.numeric(selection$criteria["FPE(n)", ]),
    row.names = NULL
  )

  list(summary = summary, directional = directional,
       lag_criteria = lag_criteria, fit = fit,
       spillover = spill, max_root = max_root)
}

crisis_spillover_results <- lapply(names(period_indices), function(label) {
  fit_spillovers(label, period_indices[[label]])
})
names(crisis_spillover_results) <- names(period_indices)

crisis_spillover_summary <- do.call(rbind, lapply(
  crisis_spillover_results, `[[`, "summary"))
crisis_spillover_directional <- do.call(rbind, lapply(
  crisis_spillover_results, `[[`, "directional"))
crisis_var_lag_criteria <- do.call(rbind, lapply(
  crisis_spillover_results, `[[`, "lag_criteria"))

tab("12_crisis_spillover_summary", crisis_spillover_summary, digits = 4,
    caption = "Full sample and matched non-crisis / GFC generalized spillover indices")
tab("12_crisis_spillover_directional", crisis_spillover_directional, digits = 4,
    caption = "Directional volatility spillovers by bank and sample")
tab("12_crisis_var_lag_criteria", crisis_var_lag_criteria, digits = 4,
    caption = "VAR lag criteria by sample; BIC (SC) selects the reported lag")

print(crisis_spillover_summary, digits = 4, row.names = FALSE)
print(crisis_spillover_directional, digits = 4, row.names = FALSE)

## ---- Figures -----------------------------------------------------------------

fig("12_crisis_spillover_index", {
  values <- crisis_spillover_summary$total_spillover_pct
  names(values) <- crisis_spillover_summary$period
  plot_values <- values[names(period_labels)]
  barplot(plot_values, names.arg = unname(period_labels),
          col = c("grey65", "steelblue", "firebrick"),
          ylim = c(0, max(plot_values, na.rm = TRUE) * 1.15),
          ylab = "total spillover index (%)",
          main = paste0("10-day generalized volatility spillovers (VAR lag by BIC)"))
  grid(nx = NA, ny = NULL, col = "grey88")
})

fig("12_crisis_net_directional", {
  comparison_periods <- c("noncrisis_pre_gfc", "crisis_gfc")
  net_by_period <- lapply(comparison_periods, function(period) {
    d <- crisis_spillover_directional[
      crisis_spillover_directional$period == period, , drop = FALSE]
    d$net_TO_pct[match(TICKERS, d$ticker)]
  })
  net_matrix <- do.call(rbind, net_by_period)
  rownames(net_matrix) <- c("Matched non-crisis", "GFC crisis")
  colnames(net_matrix) <- TICKERS
  barplot(net_matrix, beside = TRUE,
          col = c("steelblue", "firebrick"),
          ylab = "net directional spillover (percentage points)",
          main = "Net spillovers by bank: matched non-crisis vs GFC",
          legend.text = rownames(net_matrix), args.legend = list(x = "topright", bty = "n"))
  abline(h = 0, col = "grey30")
  grid(nx = NA, ny = NULL, col = "grey88")
})

saveRDS(mget(c("period_indices", "crisis_spillover_directional",
               "crisis_spillover_results", "crisis_spillover_summary",
               "crisis_var_lag_criteria", "log_range_panel", "range_panel")),
        "results/rds/12_crisis.rds")


