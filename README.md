# EBC4008 Time Series Analytics: Case 3

Modeling stock prices, returns and volatility: JPMorgan Chase (JPM), daily data from Yahoo Finance, 1999-05-04 to 2026-09-18. The brief is in `case3_assignment.pdf`.

This repository holds parts 01–05 of the analysis:

| Part | Script | What it does |
|---|---|---|
| 01 | `R/parts/01_data.R` | Loads the six-bank panel (JPM, BAC, C, WFC, GS, MS) plus XLF and the S&P 500 |
| 02 | `R/parts/02_quality_screen.R` | Checks the High/Low data; summary and stale/jump screen for y, log(y) and r |
| 03 | `R/parts/03_exploration.R` | Builds y (price), log(y) and r (log return); stylized facts, ACFs, boxplots, rolling statistics |
| 04 | `R/parts/04_unit_roots.R` | Unit-root tests in Pantula order (ADF, cross-checked with KPSS) on y, log(y) and r |
| 05 | `R/parts/05_arima.R` | ARIMA by hand from the correlogram and by `auto.arima`, for y, log(y) and r; compares the price fit with the return fit |

## Run it

From the `case3` folder:

```
Rscript R/parts/run_all.R        # every part, in order
Rscript R/parts/04_unit_roots.R  # one part; it runs the parts it needs first
```

You can also source a part in RStudio. It needs R 4.6. The packages are listed in `R/parts/00_setup.R` and missing ones install on first run. The first run takes about 5 minutes. After that the slow fits come from `results/cache/`.

- `data/case3_prices.rds` is the frozen Yahoo download. Delete it to re-download, but the numbers will then change.
- `results/tables/` and `results/figures/` hold the output of parts 01–05.
