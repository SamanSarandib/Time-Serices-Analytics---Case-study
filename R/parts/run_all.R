###############################################################################
# EBC4008 Case 3 -- run the parts
#
#   Rscript R/parts/run_all.R            every part, 01 -> 13 (no 06: price and
#                                        return ARIMA are both in 05)
#   Rscript R/parts/run_all.R 7          07 and everything that depends on it
#                                        (08, 09, 10, 13) -- use after editing 07
#   Rscript R/parts/run_all.R 7 --only   exactly 07 (later parts then say they
#                                        are out of date when they load it)
# or Source this file in RStudio (runs every part).
#
# Each part runs via run_part() in its own environment, a child of the global
# one: its objects stay there and don't reach other parts. Under Rscript that
# is all a part sees besides 00_setup.R and what it need()s; Sourced in
# RStudio it can also see whatever is in your workspace.
# R/case3_run.R is the OLDER one-file version (y = log price, price ARIMA in 5,
# returns in 6); its outputs no longer match the parts. The parts write
# results/tables, results/figures, results/case3_results.rds,
# results/parts/<part>.rds and one log per part in results/log/parts/.
# After a data refresh or a model change, also delete results/cache or bump
# CACHE_TAG in 00_setup.R: cached fits don't know the data changed.
###############################################################################

## ---- bootstrap: find R/parts/, load the shared setup -------------------------
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

local({                                  # keep run_all's own variables out of the workspace
  args <- commandArgs(TRUE)
  only <- "--only" %in% args
  args <- setdiff(args, "--only")
  todo <- if (length(args)) {
    ids <- as.integer(substr(PARTS, 1, 2))
    n   <- suppressWarnings(as.numeric(args))
    bad <- args[is.na(n) | n != round(n) | !(n %in% ids)]
    if (length(bad)) stop("no such part: ", paste(bad, collapse = ", "),
                          " -- parts are ", paste(ids, collapse = ", "),
                          " (05 covers price and return ARIMA; there is no 06), e.g.  Rscript R/parts/run_all.R 7")
    picked <- PARTS[ids %in% n]
    if (only) picked else part_downstream(picked)
  } else PARTS
  message("  [run_all] running: ", paste(todo, collapse = " "))
  t0 <- Sys.time()
  for (p in todo) run_part(p)
  message(sprintf("  [run_all] %d part(s) in %.1f min", length(todo),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
})
