###############################################################################
# EBC4008 Case 3 -- part 01: data
# Reads : data/case3_prices.rds (downloaded from Yahoo first if missing)
# Writes: tables 01_history_depth; results/rds/01_data.rds
###############################################################################

library(quantmod)   # getSymbols(); also attaches xts and zoo

# folders for the data and every part's outputs
for (d in c("data", "results/tables", "results/figures", "results/rds"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

TICKERS <- c("JPM", "BAC", "C", "WFC", "GS", "MS")
BENCH   <- c("XLF", "^GSPC")   # sector ETF, S&P 500
START   <- "1999-05-04"        # GS IPO
END     <- Sys.Date()

# crisis windows, c(start, end)
EPISODES <- list(
  dotcom  = c("2000-03-24", "2002-10-09"),  # S&P 500 peak to trough
  gfc     = c("2007-08-01", "2009-06-30"),  # quant quake / BNP freeze
  euro    = c("2011-07-01", "2012-07-26"),
  covid   = c("2020-02-19", "2020-06-30"),  # pre-covid S&P peak
  banks23 = c("2023-02-15", "2023-06-30")   # SVB/Signature
)

# download from Yahoo once, then load from disk -> prices (list of OHLCVA xts)
if (!file.exists("data/case3_prices.rds")) {
  prices <- list()
  for (s in c(TICKERS, BENCH)) {
    prices[[s]] <- getSymbols(s, src = "yahoo", from = START, to = END, auto.assign = FALSE)
    Sys.sleep(1)
  }
  names(prices) <- c(TICKERS, "XLF", "GSPC")
  saveRDS(prices, "data/case3_prices.rds")
} else {
  prices <- readRDS("data/case3_prices.rds")   # delete the file to re-download
}

# start, end, n per series -> 01_history_depth
depth <- data.frame(
  series = names(prices),
  start  = sapply(prices, function(x) format(start(x))),
  end    = sapply(prices, function(x) format(end(x))),
  n      = sapply(prices, nrow),
  row.names = NULL)
write.csv(depth, "results/tables/01_history_depth.csv")
print(depth, digits = 4)

# balanced panel: all series start on START
stopifnot(all(depth$start == START))

jpm <- prices[["JPM"]]   # xts: OHLCV + Adjusted

saveRDS(mget(c("prices", "jpm", "TICKERS", "EPISODES")), "results/rds/01_data.rds")
