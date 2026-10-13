# EBC4008 Time Series Analytics: Case 3

Modeling stock prices, returns and volatility: JPMorgan Chase (JPM), daily data from Yahoo Finance, 1999-05-04 to 2026-09-18. The brief is in `case3_assignment.pdf`.

Same analysis as `case3-bootur`, same tables and figures, but every part is a plain standalone R script: no shared setup file, no helper library, no caching. Parts pass results on through `results/rds/`.

| Part | Script | Reads | What it does |
|---|---|---|---|
| 01 | `R/01_data.R` | `data/case3_prices.rds` | Loads (or downloads) the six-bank panel (JPM, BAC, C, WFC, GS, MS) plus XLF and the S&P 500; defines the tickers and crisis episodes |
| 02 | `R/02_quality_screen.R` | 01 | Checks the High/Low data; summary and stale/jump screen for y, log(y) and r |
| 03 | `R/03_exploration.R` | 01 | Builds y (price), log(y) and r (log return); stylized facts, ACFs, boxplots, rolling statistics |
| 04 | `R/04_unit_roots.R` | 03 | Unit-root tests in Pantula order (`bootUR` wild-bootstrap ADF, cross-checked with KPSS from `tseries`) on y, log(y) and r |
| 05 | `R/05_arima.R` | 03 | ARIMA by hand from the correlogram and by `auto.arima`, for y, log(y) and r; compares the price fit with the return fit |
| 07 | `R/07_garch.R` | 01, 03 | GARCH on r: 4 variance models x 4 error densities ranked by BIC, eGARCH-std selected; sign bias, stability, persistence, news impact, 22-day forecast |
| 08 | `R/08_realized_range_har.R` | 01, 02, 07 | Parkinson realized range of JPM, `bootUR` ADF on log(RR); HAR-RR (1/5/22 days) with Newey-West se, residual checks and other horizon sets; informal comparison with the 07 GARCH |
| 09 | `R/09_monthly_rv.R` | 03, 08 | Monthly realized variance from daily returns (plain and French-Schwert-Stambaugh), `bootUR` ADF on log(RV); monthly HAR (1/3/12 months); daily vs monthly HAR |
| 11 | `R/11_spillovers.R` | 01 | Diebold-Yilmaz spillovers: VAR on the six banks' log ranges (`bootUR` Step-M test that all six are I(0)), spillover table, directional and net spillovers, robustness over p and H, 200-day rolling index |

05 also covers the returns, so there is no 06.

## Run it

Open `case3.Rproj` in RStudio (this makes the project folder the working directory), then run a script top to bottom. A part needs the parts in its "Reads" column to have run once. To run everything in order:

```
Rscript R/run_all.R
```

Install the packages once (R 4.6):

```
install.packages(c("quantmod", "bootUR", "tseries", "forecast", "FinTS", "rugarch",
                   "sandwich", "lmtest", "moments", "vars", "frequencyConnectedness"))
```

Nothing is cached, so the slow steps run every time. A full run takes about 17 minutes on a 16-core PC: 04 about 9 min (15 bootstrap ADF tests), 11 about 5 min, 05 about 2 min, 07 about 1 min; the rest take seconds. bootUR and 11's rolling index use all cores but two, so expect it to be slower on fewer cores.

- `data/case3_prices.rds` is the frozen Yahoo download; delete it to re-download (the numbers will then change).
- `results/tables/` (csv) and `results/figures/` (png) hold the output; `results/rds/` holds the objects each part passes on.
