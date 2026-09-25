###############################################################################
# EBC4008 Case 3 -- part 01: Data
#
# Inputs come from need(),
# outputs go to results/parts/01_data.rds via part_save().
#
# Needs : nothing (first part)
# Makes : jpm (needed by 02,03,08,13), prices (needed by 02, 11)
# Writes: tables 01_history_depth
###############################################################################

# bootstrap: find R/part
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

part_begin("01_data")

## ---- 1. Data ------------------------------------------------------------------

say("1. DATA") # Print banner of = lines + section title

# Why JPM: one ticker, no splices, since before 1999. Very liquid, so High/Low
# reflect volatility rather than bid-ask bounce (matters for the range in
# section 8). A bank, so 2007-09 is a real regime shift inside the series. And
# five big-bank peers give a natural panel for the spillover part.
# Why this sample: GS listed on 1999-05-04, so starting there makes the six-bank
# panel balanced from day one (same start as Diebold-Yilmaz). ~6,900 days, five
# crises, enough data for GARCH/HAR and 200-day rolling VARs.

note("Downloading (once) or loading from disk: ", DATAFILE)
if (!file.exists(DATAFILE)) { 
  # download once, then always work from the .rds (Yahoo rate-limits)
  prices <- lapply(c(TICKERS, BENCH), function(s) { Sys.sleep(1); safe_get(s, START, END) })  # combines the six bank tickers with the two benchmark
  names(prices) <- c(TICKERS, "XLF", "GSPC")   # getSymbols strips the ^ from ^GSPC
  saveRDS(prices, DATAFILE)
  message("  downloaded and saved")
} else {
  prices <- readRDS(DATAFILE)
  message("  loaded from disk (delete the file to re-download)")
}

depth <- data.frame(
  series = names(prices),
  start  = sapply(prices, function(x) format(start(x))),
  end    = sapply(prices, function(x) format(end(x))),
  n      = sapply(prices, nrow),
  row.names = NULL)
tab("01_history_depth", depth, caption = "History depth of every downloaded series")
# series      start          end      n
# JPM       1999-05-04   2026-09-18   6887
# ................

# balanced from day one, or na.omit() hides the early years later
# All must have the same start else throw error
stopifnot(all(sapply(prices, function(x) format(start(x))) == START))

jpm <- prices[[PRIMARY]] # returns the primary as the element itsself (xts)

part_save(c("jpm", "prices"))
