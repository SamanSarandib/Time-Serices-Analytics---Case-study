###############################################################################
# EBC4008 Case 3 -- part 01: Data
# Needs : nothing (first part)
# Makes : prices (list of 8 xts: Open|High|Low|Close|Volume|Adjusted), jpm (prices$JPM)
# Writes: tables 01_history_depth
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

part_begin("01_data")

# quantmod::getSymbols via safe_get()
note("Downloading (once) or loading from disk: ", DATAFILE)
if (!file.exists(DATAFILE)) { 
  prices <- lapply(c(TICKERS, BENCH), function(s) { Sys.sleep(1); safe_get(s, START, END) })
  names(prices) <- c(TICKERS, "XLF", "GSPC")  
  saveRDS(prices, DATAFILE)
  message("downloaded and saved")
} else {
  prices <- readRDS(DATAFILE)
  message("loaded from disk (delete the file to re-download)")
}

# -> 01_history_depth: series | start | end | n
depth <- data.frame(
  series = names(prices),
  start  = sapply(prices, function(x) format(start(x))),
  end    = sapply(prices, function(x) format(end(x))),
  n      = sapply(prices, nrow),
  row.names = NULL)
tab("01_history_depth", depth, caption = "History depth of every downloaded series")

# every series must start on START (balanced panel), else stop
stopifnot(all(sapply(prices, function(x) format(start(x))) == START))

jpm <- prices[[PRIMARY]] # xts, JPM OHLCV + Adjusted

part_save(c("jpm", "prices"))
