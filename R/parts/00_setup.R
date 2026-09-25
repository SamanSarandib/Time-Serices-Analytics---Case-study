###############################################################################
# EBC4008 Case 3 -- parts: shared setup
#
# Every part sources this first. Packages, paths, constants, and the toolkit
#   need(part, objs)  load objects an earlier part saved. Runs that part first
#                     if it has no saved output; says so if the output is out of
#                     date. Call it BEFORE part_begin().
#   part_begin(part)  start a part: drop its old output, seed, own log, tables
#   part_save(objs)   save what later parts use to results/parts/<part>.rds,
#                     with checksums of what it was built from
#   run_part(part)    run one part in a fresh environment (run_all.R uses it)
# Safe to source again and again
# handy for poking around:  source("R/parts/00_setup.R"); need("07_garch", "fit_g")
###############################################################################

# where R/parts/ is: this file's own location wins, then the session's last one
# R reads files relative to the working directory
# local {} => all variables vanish only normalizePath remains
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

# Switches
AUTO_INSTALL <- TRUE      # install missing packages
CACHE_ON     <- TRUE      # memoise slow fits in results/cache
CACHE_TAG    <- "v1"      # bump to invalidate every cached object
SAVE_FIGS    <- TRUE      # FALSE = draw on screen instead of png

pkgs <- c(
  "quantmod", "xts", "zoo",                 # data
  "urca", "forecast", "tseries",            # unit roots, ARIMA, Jarque-Bera
  "FinTS", "rugarch",                       # ARCH-LM, GARCH
  "sandwich", "lmtest",                     # HAC
  "vars", "frequencyConnectedness",         # VAR, Diebold-Yilmaz
  "moments"                                 # skewness, kurtosis
)
# Check which packages is missing in installed (True => missing)
# Then either install automatically or manually if any is missing
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs)) {
  if (AUTO_INSTALL) {
    message("Installing: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")   # Rscript has no mirror set
  } else {
    stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
         "\nSet AUTO_INSTALL <- TRUE or install them yourself.")
  }
}
invisible(lapply(pkgs, library, character.only = TRUE))

# case3/ is two levels above R/parts/
ROOT <- normalizePath(file.path(PARTS_DIR, "..", ".."), winslash = "/")
if (!file.exists(file.path(ROOT, "R", "case3_functions.R")))
  stop("expected R/case3_functions.R in ", ROOT)
source(file.path(ROOT, "R", "case3_functions.R"))

# every table that goes through tab() is also kept in TABLES (emptied per part)
# Upgraded tab
if (!exists("TABLES")) TABLES <- list()
tab <- local({
  write_tab <- tab
  function(name, x, ...) { TABLES[[name]] <<- as.data.frame(x); write_tab(name, x, ...) }
})

CFG <- case3_config(root = ROOT, cache_on = CACHE_ON, cache_tag = CACHE_TAG,
                    save_figs = SAVE_FIGS, seed = 1)
CFG$partdir    <- file.path(CFG$resdir, "parts")     # each part's saved outputs
CFG$partlogdir <- file.path(CFG$logdir, "parts")     # each part's console log
for (d in c(CFG$partdir, CFG$partlogdir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# Stock-exchange codes
TICKERS  <- c("JPM", "BAC", "C", "WFC", "GS", "MS")
# XLF is an exchange-traded fund holding the whole US financial sector; ^GSPC is the S&P 500 index. Downloaded, but only used in the history table
BENCH    <- c("XLF", "^GSPC")         
# Stock studied in detail
PRIMARY  <- "JPM"
# Goldman Sachs' IPO
START    <- "1999-05-04"   
# today
END      <- Sys.Date()
DATAFILE <- file.path(CFG$datadir, "case3_prices.rds")
# crisis dates
EPISODES <- case3_episodes()           

# HAC SEs at several truncations as one table (parts 08 and 09)
hac_se_by_lag <- function(fit, lags) do.call(rbind, lapply(lags, function(L) {
  V <- NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE)
  data.frame(hac_lag = L, t(sqrt(diag(V))), check.names = FALSE)
}))


## ---- the parts toolkit -------------------------------------------------------

PARTS <- c("01_data", "02_quality_screen", "03_exploration", "04_unit_roots",
           "05_arima", "07_garch", "08_realized_range_har",
           "09_monthly_rv", "10_out_of_sample", "11_spillovers", "12_crisis",
           "13_gather")

part_out    <- function(part) file.path(CFG$partdir, paste0(part, ".rds"))  # Saved outputs of the parts
part_script <- function(part) file.path(PARTS_DIR, paste0(part, ".R"))      # Code

# a checkout can hold only some of the parts (e.g. 01-05): keep the ones present
PARTS <- PARTS[file.exists(part_script(PARTS))]

# Reads the script and determines what the part needs from need() through RE
part_deps <- function(part) {
  src <- readLines(part_script(part), warn = FALSE)
  unique(regmatches(src, regexpr('(?<=need\\(")[0-9]{2}_[a-z_]+', src, perl = TRUE)))
}
# a part plus every part that depends on it, directly or not, in run order
part_downstream <- function(parts) {
  todo <- parts
  repeat {
    more <- setdiff(PARTS[vapply(PARTS, function(p) any(part_deps(p) %in% todo), logical(1))], todo)
    if (!length(more)) break
    todo <- c(todo, more)
  }
  PARTS[PARTS %in% todo]
}

# checksums of everything a part's output depends on
# makes a fingerprint of everything a part's output needs, part_save() stores the list, while need() makes a fresh list to compare checksums
provenance <- function(part) {
  deps <- part_deps(part)
  f <- c(script = part_script(part), setup = file.path(PARTS_DIR, "00_setup.R"),
         functions = file.path(ROOT, "R", "case3_functions.R"),
         if (part == "01_data") c(data = DATAFILE),
         if (length(deps)) setNames(part_out(deps), deps))   # paste0(character(0), ...) isn't empty
  setNames(unname(tools::md5sum(f)), names(f))
}
# what changed since a part's output was saved (empty = up to date)
part_stale <- function(part, saved) {
  old <- saved$.provenance
  if (is.null(old)) return("no record of its inputs")
  now <- provenance(part)
  same <- vapply(names(old), function(n) identical(unname(now[n]), unname(old[n])), logical(1))
  ch   <- union(names(old)[!same], setdiff(names(now), names(old)))
  lab <- c(script = paste0(part, ".R"), setup = "00_setup.R", functions = "case3_functions.R",
           data = "the data file")
  ifelse(ch %in% names(lab), lab[ch], paste0("output of ", ch))
}

# Run one part in a fresh environment (a child of the global one): its objects
# stay there, and other parts' objects are not in it. The global workspace is
# still visible, so the clean check that a need() list is complete is a fresh
# Rscript run.
run_part <- function(part) {
  if (!part %in% PARTS) stop("unknown part: ", part)
  sys.source(part_script(part), envir = new.env(parent = globalenv()))
  invisible(part)
}

# Load `objs` saved by `part` into the caller's environment.
need <- function(part, objs, envir = parent.frame()) {
  out <- part_out(part)
  if (!file.exists(out)) {
    message("  [need] ", part, " has no saved output (never run, or its last run failed) -- running it first")
    run_part(part)
  }
  saved <- readRDS(out)
  stale <- part_stale(part, saved)
  if (length(stale))
    message("  [need] NOTE: ", part, " is out of date (changed since: ", paste(stale, collapse = ", "),
            "). Rerun it and what depends on it:  Rscript R/parts/run_all.R ", as.integer(substr(part, 1, 2)))
  miss <- setdiff(objs, names(saved))
  if (length(miss)) stop(part, " doesn't save: ", paste(miss, collapse = ", "))
  list2env(saved[objs], envir = envir)
  invisible(objs)
}

# the tables a part wrote (13_gather rebuilds TABLES from these). Never runs a
# part: 13_gather need()s every part before its part_begin().
part_tables <- function(part) {
  if (!file.exists(part_out(part))) stop(part, " has no saved output -- run it first")
  readRDS(part_out(part))$.tables
}

if (!exists(".part_state")) .part_state <- new.env()

# Resets the env (encapsulated in .part_state)
# deletes old files
# opens and empties the log file
# clears the TABLES
part_begin <- function(part) {
  while (sink.number() > 0) sink()          # leftover from a crashed part
  if (!is.null(.part_state$log)) try(close(.part_state$log), silent = TRUE)
  unlink(part_out(part))                    # a failed run must not leave the old output behind
  .part_state$part <- part
  .part_state$log  <- file(file.path(CFG$partlogdir, paste0(part, ".txt")), open = "wt")
  sink(.part_state$log, split = TRUE)
  TABLES <<- list()
  set.seed(1)                               # every part starts from the same seed
  invisible(part)
}

# urca tests keep the function frame they ran in (lag matrix, helper fits),
# through the terms of their lm() fit and of its model frame. Summaries,
# statistics and critical values don't use it, so drop it before saving.
slim <- function(o) {
  if (isS4(o) && "testreg" %in% slotNames(o) && !is.null(o@testreg$terms)) {
    attr(o@testreg$terms, ".Environment") <- globalenv()
    if (!is.null(o@testreg$model))
      attr(attr(o@testreg$model, "terms"), ".Environment") <- globalenv()
  }
  o
}

part_save <- function(objs, envir = parent.frame(), tables = TRUE) {
  part <- .part_state$part
  miss <- objs[!vapply(objs, exists, logical(1), envir = envir, inherits = FALSE)]
  if (length(miss)) stop(part, ": outputs not created: ", paste(miss, collapse = ", "))
  saved <- lapply(mget(objs, envir = envir), slim)
  saved$.tables     <- if (tables) TABLES else list()
  saved$.provenance <- provenance(part)
  saveRDS(saved, part_out(part))
  message(sprintf("  [part] %s saved: %d objects, %d tables", part, length(objs), length(saved$.tables)))
  sink(); close(.part_state$log); .part_state$log <- NULL
  invisible(objs)
}
