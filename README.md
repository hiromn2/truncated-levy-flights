# truncated-levy-flights

# Bypassing the Truncation Problem of Truncated Lévy Flights

**Physica A: Statistical Mechanics and its Applications, 559 (2020) 125035**  
Matsushita, Da Silva, Da Fonseca, **Nagata**  
[📄 Read the paper](https://doi.org/10.1016/j.physa.2020.125035)

---

## What this is

Standard financial models assume Gaussian returns — but real markets don't behave that way. 
Flash crashes, currency crises, and tail events happen far more often than a Gaussian would predict.

Truncated Lévy Flights (TLFs) were proposed to fix this, but they face a fundamental problem: 
every new extreme event makes the current truncation cutoff obsolete. You model the last crisis, 
not the next one.

This paper proposes a solution: instead of chasing an ideal cutoff, we find a **power law** 
between the truncation cutoff ℓ and the return standard deviation σ:
```
ℓₜ = ζ · σₜᵝ
```

Both ℓ and σ are time-varying, β is a fixed sensitivity index, and the distribution of returns 
is left unknown. This gives a model that is robust to future extreme events by design.

---

## Why it matters for risk modeling

- **VaR and CVaR**: Standard models underestimate tail probabilities. Our Fréchet-based tail 
  estimator provides sharper extreme quantile estimates directly applicable to risk limits and 
  margin calculations.
  
- **Flash crash robustness**: The model was fitted through the January 15, 2015 Swiss franc 
  crisis — when CHF/USD moved ~15% in minutes. The standardized returns show no outlier on 
  that day, meaning the model absorbed the shock rather than breaking.

- **Distribution-free**: No parametric assumption on the return distribution. The power law 
  holds for Cauchy, Gaussian, and mixture distributions alike.

---

## Data

4,842,688 tick-by-tick bid prices of CHF/USD  
**Period:** May 30, 2000 – June 12, 2015 (4,588 trading days)  
**Source:** Tick Data, LLC  

The dataset covers the full 2015 Swiss National Bank announcement shock.

---

## Key results

| Parameter | Estimate |
|-----------|----------|
| Sensitivity index β̂ | 0.85 |
| Fréchet tail index ξ̂ | ~8.5 |
| Scale parameter γ̂ | 0.44 |

β < 1 indicates the CHF/USD series has **multiple truncation regimes** — the single-cutoff 
TLF assumption fails on real data with extreme moves.

The Fréchet fit allows tail probability estimation beyond observed data. For example, a 
normalized threshold ζ = 3 (never observed in 15 years) is estimated at P ≈ 0.0000418.

---

## Repo contents
```
/
├── data/           # Summary statistics (raw tick data not redistributable)
├── estimation.R    # Power law estimation: β̂ via log-log regression
├── tail_model.R    # Fréchet tail fitting and probability estimation
├── figures.R       # Replication of paper figures
└── README.md
```

---

## Related work

This is part of a series on heavy-tailed financial modeling:

- [Retrodicting with the truncated Lévy flight](https://doi.org/10.1016/j.cnsns.2022.106721) — *Nonlinear Science, 2022*  
- [The duration of historical pandemics](https://doi.org/10.1016/j.cnsns.2022.106461) — EVT applied to pandemic duration

---

## Citation
```bibtex
@article{matsushita2020truncated,
  title   = {Bypassing the truncation problem of truncated Lévy flights},
  author  = {Matsushita, Raul and Da Silva, Sergio and Da Fonseca, Regina and Nagata, Mateus},
  journal = {Physica A: Statistical Mechanics and its Applications},
  volume  = {559},
  pages   = {125035},
  year    = {2020},
  doi     = {10.1016/j.physa.2020.125035}
}
```