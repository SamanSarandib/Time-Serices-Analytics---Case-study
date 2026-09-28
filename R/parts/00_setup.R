###############################################################################
# EBC4008 Case 3 -- shared setup (every part sources it first; safe to re-source)
# Needs:  R/case3_functions.R
# Makes:  packages loaded, PARTS_DIR, ROOT, CFG, PARTS, constants (TICKERS ... EPISODES), tab() + TABLES,
#         hac_se_by_lag(), toolkit: need(), part_begin(), part_save(), run_part(), part_tables(), part_downstream()
# Writes: folders only: data/, results/{figures,tables,cache,log}/, results/parts/, results/log/parts/
###############################################################################

# find R/parts/ (Rscript, source() or RStudio Ctrl+Enter) -> PARTS_DIR
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
# install missing packages (or stop if AUTO_INSTALL is FALSE), then load them all
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

# ROOT = case3/ (two levels above R/parts/); load the function library
ROOT <- normalizePath(file.path(PARTS_DIR, "..", ".."), winslash = "/")
if (!file.exists(file.path(ROOT, "R", "case3_functions.R")))
  stop("expected R/case3_functions.R in ", ROOT)
source(file.path(ROOT, "R", "case3_functions.R"))

# wrap tab(): also keep each table as TABLES[[name]] (data.frame; emptied per part)
if (!exists("TABLES")) TABLES <- list()
tab <- local({
  write_tab <- tab
  function(name, x, ...) { TABLES[[name]] <<- as.data.frame(x); write_tab(name, x, ...) }
})

# CFG: list of paths (datadir, resdir, figdir, tabdir, cachedir, logdir) + switches; makes the folders
CFG <- case3_config(root = ROOT, cache_on = CACHE_ON, cache_tag = CACHE_TAG,
                    save_figs = SAVE_FIGS, seed = 1)
CFG$partdir    <- file.path(CFG$resdir, "parts")     # each part's saved outputs
CFG$partlogdir <- file.path(CFG$logdir, "parts")     # each part's console log
for (d in c(CFG$partdir, CFG$partlogdir)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# the six bank tickers
TICKERS  <- c("JPM", "BAC", "C", "WFC", "GS", "MS")
# XLF = financial-sector ETF, ^GSPC = S&P 500: downloaded, only used in 01_history_depth
BENCH    <- c("XLF", "^GSPC")         
# the bank studied in detail
PRIMARY  <- "JPM"
# sample: Goldman Sachs' IPO to today
START    <- "1999-05-04"   
END      <- Sys.Date()
DATAFILE <- file.path(CFG$datadir, "case3_prices.rds")
# crisis windows -> named list (dotcom, gfc, euro, covid, banks23), each c(start, end)
EPISODES <- case3_episodes()           

# Newey-West HAC SEs (no prewhitening) of an lm `fit` at each lag in `lags`; used by parts 08, 09
# -> data.frame(hac_lag, one SE column per coefficient), one row per lag
hac_se_by_lag <- function(fit, lags) do.call(rbind, lapply(lags, function(L) {
  V <- NeweyWest(fit, lag = L, prewhite = FALSE, adjust = TRUE)
  data.frame(hac_lag = L, t(sqrt(diag(V))), check.names = FALSE)
}))


## ---- the parts toolkit -------------------------------------------------------

# all parts, in run order
PARTS <- c("01_data", "02_quality_screen", "03_exploration", "04_unit_roots",
           "05_arima", "07_garch", "08_realized_range_har",
           "09_monthly_rv", "10_out_of_sample", "11_spillovers", "12_crisis",
           "13_gather")

part_out    <- function(part) file.path(CFG$partdir, paste0(part, ".rds"))  # -> results/parts/<part>.rds
part_script <- function(part) file.path(PARTS_DIR, paste0(part, ".R"))      # -> R/parts/<part>.R

# keep only the parts whose script exists
PARTS <- PARTS[file.exists(part_script(PARTS))]

# regex the script for need("<part>" calls -> character vector of the parts it depends on
part_deps <- function(part) {
  src <- readLines(part_script(part), warn = FALSE)
  unique(regmatches(src, regexpr('(?<=need\\(")[0-9]{2}_[a-z_]+', src, perl = TRUE)))
}
# repeat part_deps() until nothing new -> `parts` + every part depending on them, in run order
part_downstream <- function(parts) {
  todo <- parts
  repeat {
    more <- setdiff(PARTS[vapply(PARTS, function(p) any(part_deps(p) %in% todo), logical(1))], todo)
    if (!length(more)) break
    todo <- c(todo, more)
  }
  PARTS[PARTS %in% todo]
}

# md5sum of the part's script, 00_setup.R, case3_functions.R, the data file (01 only), deps' .rds
# -> named character vector (script, setup, functions, [data], <dep parts>); part_save() stores it
provenance <- function(part) {
  deps <- part_deps(part)
  f <- c(script = part_script(part), setup = file.path(PARTS_DIR, "00_setup.R"),
         functions = file.path(ROOT, "R", "case3_functions.R"),
         if (part == "01_data") c(data = DATAFILE),
         if (length(deps)) setNames(part_out(deps), deps))   # guard: part_out(character(0)) isn't empty
  setNames(unname(tools::md5sum(f)), names(f))
}
# saved vs current provenance() -> labels of the changed inputs (length 0 = up to date)
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

# sys.source() the part's script in a fresh env (child of global); used by need() and run_all.R
run_part <- function(part) {
  if (!part %in% PARTS) stop("unknown part: ", part)
  sys.source(part_script(part), envir = new.env(parent = globalenv()))
  invisible(part)
}

# load `objs` from results/parts/<part>.rds into the caller's env, e.g. need("07_garch", "fit_g")
# runs the part first if it has no .rds, notes if stale; call it before part_begin()
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

# -> named list of the data.frames `part` wrote with tab() (used by 13_gather); never runs the part
part_tables <- function(part) {
  if (!file.exists(part_out(part))) stop(part, " has no saved output -- run it first")
  readRDS(part_out(part))$.tables
}

# .part_state: env holding the current part name and its log connection
if (!exists(".part_state")) .part_state <- new.env()

# start a part: delete its old .rds, sink the console to results/log/parts/<part>.txt,
# empty TABLES, set.seed(1)
part_begin <- function(part) {
  while (sink.number() > 0) sink()          # leftover from a crashed part
  if (!is.null(.part_state$log)) try(close(.part_state$log), silent = TRUE)
  unlink(part_out(part))                    # so a failed run leaves no old output
  .part_state$part <- part
  .part_state$log  <- file(file.path(CFG$partlogdir, paste0(part, ".txt")), open = "wt")
  sink(.part_state$log, split = TRUE)
  TABLES <<- list()
  set.seed(1)
  invisible(part)
}

# urca test objects: point their lm terms' .Environment at globalenv (drops the test's function frame)
# -> same object, smaller .rds; anything else returned as is
slim <- function(o) {
  if (isS4(o) && "testreg" %in% slotNames(o) && !is.null(o@testreg$terms)) {
    attr(o@testreg$terms, ".Environment") <- globalenv()
    if (!is.null(o@testreg$model))
      attr(attr(o@testreg$model, "terms"), ".Environment") <- globalenv()
  }
  o
}

# end a part: saveRDS(list(<objs> via slim(), .tables, .provenance)) -> results/parts/<part>.rds,
# stops if an obj is missing, closes the log
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
