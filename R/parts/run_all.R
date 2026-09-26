###############################################################################
# EBC4008 Case 3 -- run the parts
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
