# Bypassing the Truncation Problem of Truncated Lévy Flights

**Physica A: Statistical Mechanics and its Applications, 559 (2020) 125035**  
Matsushita, Da Silva, Da Fonseca, **Nagata**  
[📄 Paper](https://doi.org/10.1016/j.physa.2020.125035) · [👤 Author](https://hiromn2.github.io)

---

## The problem

Truncated Lévy Flights (TLFs) were introduced to model leptokurtic financial returns without the infinite-variance pathology of the full Lévy stable distribution. But they face a structural impasse: the truncation cutoff ℓ must be set by hand, and every new extreme event makes the current cutoff obsolete. Each modeling cycle incorporates the last crisis and is immediately invalidated by the next one. This is the problem of induction applied to tail risk.

## The solution

We show that the truncation cutoff ℓ and the return standard deviation σ satisfy a stable **power law**:

$$\ell_t = \zeta \cdot \sigma_t^\beta$$

Both ℓ and σ are time-varying and estimable from intraday data. The ratio ζ = ℓ / σ^β is time-invariant by construction, and β ∈ [0,1] is a sensitivity index. The paper proves this relationship holds for any symmetric distribution with finite density at zero (Proposition 1), and that β = 1 when σ/ℓ → c > 0 as ℓ → 0 (Proposition 2).

The normalized bound $\hat{\zeta}_t = \hat{\ell}_t / \hat{\sigma}_t^{\hat{\beta}}$ then follows a **Fréchet (Type II extreme value) distribution**:

$$F(z) = \exp\!\left[-\gamma\, z^{-\xi}\right], \quad z > \zeta_0$$

which yields closed-form tail probability estimates beyond any previously observed threshold.

---

## The return definition

Standard log returns ∇ ln X_t are not used here. Instead, the paper defines the **return over the mean**:

$$r_t = \frac{X_t}{m_t} - 1$$

where m_t = E[X_t] is the time-varying daily mean price, estimated from intraday data as the sample average. This makes r_t first-order stationary with zero mean, while its variance σ²_t is left free to vary over time. The normalized return is then:

$$Z_t = \frac{r_t}{\sigma_t^\beta}$$

After normalization, Z_t has a time-invariant bound ζ despite the heavy tails in r_t.

---

## Data

| | |
|---|---|
| Asset | USD/CHF bid price |
| Period | 30 May 2000 – 12 June 2015 |
| Observations | 4,842,688 ticks across 4,588 trading days |
| Source | Tick Data, LLC |

The dataset spans the January 15, 2015 Swiss National Bank shock, when CHF/USD moved from 1.02 to 0.84 within a single trading day — one of the largest single-day FX moves in recent history.

---

## Key results

**Theoretical benchmark.** For both the truncated Cauchy (α=1) and truncated Normal (α=2), Propositions 1–2 give β = 1 and ℓ ≈ √3 σ analytically.

**Empirical finding.** On the CHF/USD data:

| Parameter | Symbol | Estimate |
|-----------|--------|----------|
| Sensitivity index (log-log slope) | β̂ | 0.85 |
| Fréchet tail index | ξ̂ | 8.52 (NLS) / 8.43 (CDF fit) |
| Fréchet scale parameter | γ̂ | 0.44 |

**β̂ < 1** is the central empirical result. Because both the Cauchy and Gaussian TLFs give β = 1 in theory, the finding β̂ = 0.85 < 1 indicates that the CHF/USD series cannot be described by a single sharp truncation — the mixture structure (Example 2.3 in the paper) applies. The power law nonetheless holds, providing a stable regularity that absorbs all future extreme events without requiring the cutoff to be re-estimated.

**Flash-crash robustness.** The intraday return series r_t shows a clear outlier on January 15, 2015 (r_t = 0.133). The normalized return series Z_t shows no outlier on that day. The model absorbed the shock.

**Tail probability estimates** beyond the observed range (Table 1 in the paper):

| ζ | Empirical CDF | Fréchet model |
|---|--------------|---------------|
| 1.0 | 0.3376 | 0.3582 |
| 1.5 | 0.0118 | 0.0144 |
| 2.0 | 0.0011 | 0.0013 |
| 2.5 | 0.0002 | 0.0002 |
| 3.0 | 0 (never observed) | 0.0000418 |
| 4.0 | 0 | 0.0000037 |
| 5.0 | 0 | 0.0000006 |

---

## Repo structure

```
/
├── R/
│   ├── 00_load_data.R       # Loads and combines all yearly CSV files
│   ├── 01_log_returns.R     # Lévy scaling: estimates Lévy index from P(0|Δt)
│   ├── 02_main_analysis.R   # Main pipeline: β̂, ζ̂_t, Fréchet fit, NLS tail index
│   ├── 03_proposition1.R    # Proposition 1–2 verification: truncated Cauchy and Normal
│   ├── 04_figures.R         # All publication figures (ggplot2)
│   └── 05_animations.R      # Animated figures (gganimate)
├── data/
│   └── README.md            # Expected file names and column layout
└── output/
    └── figures/             # Generated PNGs and GIFs (git-ignored)
```

---

## Quick start

**1. Install dependencies**

```r
install.packages(c(
  "ggplot2", "extrafont", "grid", "scales",
  "stable", "rmutil", "moments", "pracma",
  "gganimate", "magick"
))

# Windows only — run once after installing extrafont:
# extrafont::font_import()
```

**2. Add data**

Place the TickData CSV files in `data/`. See `data/README.md` for the expected filenames and column layout. The raw tick data is proprietary and not redistributed here.

**3. Run in order**

```r
source("R/00_load_data.R")     # loads data.1 into memory
source("R/01_log_returns.R")   # estimates Lévy stable index
source("R/02_main_analysis.R") # estimates β, fits Fréchet, NLS tail index ξ̂
source("R/03_proposition1.R")  # Propositions 1–2 (standalone, no data needed)
source("R/04_figures.R")       # saves PNGs to output/figures/
source("R/05_animations.R")    # saves GIFs to output/figures/
```

> Scripts `04` and `05` use objects created by `02`. Run them in the same R session, or uncomment the `source()` calls at the top of each file.

---

## Related work

- [Retrodicting with the truncated Lévy flight](https://doi.org/10.1016/j.cnsns.2022.106721) — *Communications in Nonlinear Science and Numerical Simulation, 2022*
- [The duration of historical pandemics](https://doi.org/10.1016/j.cnsns.2022.106461) — EVT applied to pandemic duration data

---

## Citation

```bibtex
@article{matsushita2020truncated,
  title   = {Bypassing the truncation problem of truncated {L\'evy} flights},
  author  = {Matsushita, Raul and Da Silva, Sergio and Da Fonseca, Regina and Nagata, Mateus},
  journal = {Physica A: Statistical Mechanics and its Applications},
  volume  = {559},
  pages   = {125035},
  year    = {2020},
  doi     = {10.1016/j.physa.2020.125035}
}
```
