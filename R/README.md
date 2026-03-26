# Bypassing the Truncation Problem of Truncated Lévy Flights

Replication code for the paper:

> **Bypassing the Truncation Problem of Truncated Lévy Flights**  
> Mateus Hiro Nagata  
> *Related to: "On the Boundedness of an Exchange Rate"* (USD/CHF intraday, 2000–2015)

Published papers from this project available on [Google Scholar](https://scholar.google.com).

---

## Overview

This project models intraday USD/CHF exchange rate returns using **truncated Lévy flights** — a class of heavy-tailed distributions that are stable over intermediate time scales but truncated at the extremes. The key contribution is a method to estimate the Lévy stable index `α` from the log-price process without being biased by the truncation, by exploiting the scaling of `P(0|Δt)` across multiple time scales.

### Mathematical core

Let `X_t` be the USD/CHF bid price at minute `t`. Define log-returns at time scale Δt:

```
w(Δt) = log X_{t+Δt} - log X_t
```

Under a (truncated) Lévy flight, the characteristic function of `w(Δt)` scales as:

```
φ(k, Δt) = exp( i·ω·Δt^β·k − γ·Δt·|k|^α )
```

where `α ∈ (0,2]` is the **Lévy stable index**, `β` is the **drift exponent**, and `γ` is the scale parameter.

**Key identification strategy:** The density at zero satisfies:

```
P(0 | Δt) ∝ Δt^{−1/α}
```

So a log-log regression of `P̂(0|Δt)` on `Δt` gives a clean estimate of `α`, completely bypassing the truncation bias in tail estimators like the Hill estimator.

The **Proposition 1** file additionally verifies variance properties of the truncated Cauchy and truncated Normal distributions analytically.

---

## Repo structure

```
truncated-levy-flights/
├── README.md
├── R/
│   ├── 00_load_data.R       # Data loading utility (all years, one call)
│   ├── 01_log_returns.R     # Lévy scaling: log-returns at multiple Δt
│   ├── 02_main_analysis.R   # Main pipeline: α, γ, Fréchet EVT fit
│   ├── 03_proposition1.R    # Analytical: truncated Cauchy & Normal variance
│   ├── 04_figures.R         # Publication figures (ggplot2)
│   └── 05_animations.R      # Animated figures (gganimate)
├── data/
│   └── README.md            # Data source and format instructions
└── output/
    └── figures/             # PNG outputs (git-ignored)
```

---

## Data

Data is **1-minute OHLC bid/ask quotes** for USD/CHF from **TickData Inc.**, covering January 2000 – June 2015. It is not redistributed here due to licensing.

See `data/README.md` for the expected file naming convention and column layout.

---

## Requirements

```r
install.packages(c(
  "ggplot2", "extrafont", "grid", "scales",
  "stable", "rmutil", "moments", "pracma",
  "gganimate", "magick"
))
```

> **Note on fonts:** Scripts use Times New Roman. On Windows, run `extrafont::font_import()` once. On Mac/Linux, the code falls back gracefully to the default font.

---

## How to run

1. Clone this repo.
2. Place the TickData CSV files in `data/` (see `data/README.md`).
3. Run scripts in order: `00_` → `01_` → `02_` → `03_` → `04_` → `05_`.

Each script is self-contained and sources `00_load_data.R` at the top.
