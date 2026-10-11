# Phase 2b: score reliability and power of the primary dispersion call

Run `outputs/phase2b_power` (scripts/phase2b_power.py, config/phase2b.json). NTC cells only; no target data were used. Recomputed NTC scores match Phase 2 (max abs diff 2.3e-07).

## Split-half reliability (NTC cells)
Program genes split at random into halves (200 seeded splits); each half is scored with the same expression-matched background procedure. Reliability = Spearman-Brown 2r/(1+r) of the half-score correlation: the share of score variance reproducible across gene halves (biology plus any shared cell-level technical factor), the rest being gene-sampling noise.

| Readout | Role | Genes | Condition | Reliability median [min, max] | Raw (no background) median |
|---|---|---:|---|---|---:|
| Tsoi_undifferentiated | secondary | 118 | DMSO | 0.39 [0.30, 0.44] | 0.47 |
| Tsoi_undifferentiated | secondary | 118 | Vem | 0.32 [0.23, 0.39] | 0.32 |
| Tsoi_neural_crest_like | secondary | 66 | DMSO | 0.27 [0.04, 0.39] | 0.46 |
| Tsoi_neural_crest_like | secondary | 66 | Vem | 0.35 [0.22, 0.42] | 0.52 |
| Tsoi_transitory | secondary | 29 | DMSO | -0.02 [-0.13, 0.11] | 0.07 |
| Tsoi_transitory | secondary | 29 | Vem | 0.14 [-0.04, 0.24] | 0.23 |
| Tsoi_melanocytic | secondary | 62 | DMSO | 0.09 [0.01, 0.17] | 0.18 |
| Tsoi_melanocytic | secondary | 62 | Vem | 0.11 [0.01, 0.19] | 0.20 |
| Kanai2024_YAP_TAZ_up | secondary | 80 | DMSO | 0.61 [0.53, 0.65] | 0.69 |
| Kanai2024_YAP_TAZ_up | secondary | 80 | Vem | 0.32 [0.24, 0.40] | 0.51 |
| Kanai2024_YAP_TAZ_down | exploratory | 52 | DMSO | 0.33 [0.20, 0.40] | 0.42 |
| Kanai2024_YAP_TAZ_down | exploratory | 52 | Vem | 0.19 [0.08, 0.27] | 0.40 |
| Fufa2019_AZD_up_SOX10_bound_pigmentation | primary | 45 | DMSO | 0.01 [-0.09, 0.10] | 0.31 |
| Fufa2019_AZD_up_SOX10_bound_pigmentation | primary | 45 | Vem | 0.12 [-0.02, 0.23] | 0.37 |

**Fufa (primary).** Background-adjusted reliability is 0.01 (DMSO) and 0.12 (Vem); without the background it is 0.31 (DMSO) and 0.37 (Vem). The two gene halves of the raw score share a cell-level component, but the expression-matched background removes it, and what remains does not replicate across gene halves. In NTC cells, nearly all between-cell variance of the Fufa score as defined is therefore gene-sampling noise. Population means still carry signal (Phase 2 mean contrasts); single-cell dispersion of this score carries little.

## Power of the primary dispersion call (Fufa)
Pseudo-targets of n NTC cells vs the remaining NTC cells; 200 replicates per scenario, B = 5000 bootstrap draws (Phase 2: 5000). Test: log variance ratio bootstrap p < alpha ("lvr"), and the full primary call that also requires the mean CI to include 0 ("call"). alpha 0.05 is an optimistic bound for BH across 56 targets; 0.05/56 = 8.9e-04 (Bonferroni) is the bound for a single true hit.

Null rejection rates (lvr test, alpha 0.05): DMSO n=90: 0.070, DMSO n=209: 0.025, DMSO n=396: 0.035, Vem n=98: 0.065, Vem n=170: 0.080, Vem n=238: 0.060 (mean 0.056; Monte Carlo SE about 0.015 each). The percentile bootstrap test is close to nominal.

### Minimal detectable true log variance ratio at 80% power (scale model)
Linear interpolation on the scale-model grid (0, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5). Signal scale: the implied log variance ratio of the reproducible signal, log(1 + (e^MDE - 1)/reliability), assuming equal noise variance in target and NTC cells; median reliability, with the range at the split p95 and p05 reliabilities. It is not estimable when the split p05 reliability is <= 0.

| Condition | n (quantile) | MDE, alpha 0.05 | MDE, Bonferroni | Full call, Bonferroni | Bonferroni MDE, signal scale |
|---|---|---:|---:|---:|---|
| DMSO | 90 (p10) | 0.48 | 0.74 | 0.74 | not estimable (reliability ~ 0) |
| DMSO | 209 (p50) | 0.30 | 0.48 | 0.48 | not estimable (reliability ~ 0) |
| DMSO | 396 (p90) | 0.25 | 0.42 | 0.42 | not estimable (reliability ~ 0) |
| Vem | 98 (p10) | 0.46 | 0.71 | 0.71 | 2.28 [1.87, 3.83] |
| Vem | 170 (p50) | 0.38 | 0.57 | 0.58 | 2.02 [1.62, 3.54] |
| Vem | 238 (p90) | 0.30 | 0.46 | 0.46 | 1.79 [1.41, 3.28] |

### Mixture scenarios on the same scale (median n)
A fraction pi of cells shifted by +d and pi by -d (d in NTC SD units); true log variance ratio log(1 + 2 pi d^2). Cells: true lvr, then power of the lvr test at alpha 0.05 / Bonferroni.

| Condition | pi | d = 0.25 | d = 0.5 | d = 1 | d = 1.5 |
|---|---:|---|---|---|---|
| DMSO | 0.1 | 0.01: 0.04 / 0.01 | 0.05: 0.07 / 0.01 | 0.18: 0.41 / 0.02 | 0.37: 0.93 / 0.47 |
| DMSO | 0.2 | 0.02: 0.04 / 0.00 | 0.10: 0.12 / 0.00 | 0.34: 0.90 / 0.37 | 0.64: 1.00 / 1.00 |
| DMSO | 0.3 | 0.04: 0.06 / 0.00 | 0.14: 0.28 / 0.00 | 0.47: 1.00 / 0.85 | 0.85: 1.00 / 1.00 |
| DMSO | 0.5 | 0.06: 0.08 / 0.00 | 0.22: 0.55 / 0.07 | 0.69: 1.00 / 1.00 | 1.18: 1.00 / 1.00 |
| Vem | 0.1 | 0.01: 0.07 / 0.00 | 0.05: 0.09 / 0.00 | 0.18: 0.39 / 0.03 | 0.37: 0.85 / 0.33 |
| Vem | 0.2 | 0.02: 0.07 / 0.01 | 0.10: 0.11 / 0.00 | 0.34: 0.83 / 0.25 | 0.64: 1.00 / 0.96 |
| Vem | 0.3 | 0.04: 0.07 / 0.00 | 0.14: 0.23 / 0.01 | 0.47: 0.98 / 0.69 | 0.85: 1.00 / 1.00 |
| Vem | 0.5 | 0.06: 0.07 / 0.00 | 0.22: 0.48 / 0.05 | 0.69: 1.00 / 1.00 | 1.18: 1.00 / 1.00 |

## How to report the Phase 2 null
Largest observed |Fufa log variance ratio| in Phase 2: DMSO 0.23, Vem 0.32.

Suggested wording: "No knockdown changed the between-cell dispersion of the Fufa score without a mean shift. At the median target size the design had 80% power, after multiplicity control, only for log variance ratios of about 0.48 (DMSO) / 0.57 (Vem), and the score's split-half reliability in NTC cells is only 0.01 (DMSO) / 0.12 (Vem). The null result therefore means that heterogeneity of the Fufa response is not measurable with this single-cell score; it is not evidence that the response is homogeneous."

Figure: `reports/figures/phase2b_power_curves.png`. Tables: `phase2b_reliability.csv`, `phase2b_power.csv`.
