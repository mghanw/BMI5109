# Held-out evaluation of the CSCN-SEM knockdown predictions

Run: 2026-10-21T18:00:06, code at commit 11985fa; rules in `config/evaluation_run.json` (frozen before the run), endpoints in `config/evaluation_endpoints.json` (authoritative).

## Deviations

- None: the evaluation was run once with the committed code.

## Design

- Observed response: target - NTC mean of log1p CP10k (whole-transcriptome totals) over the 249 panel genes other than the perturbed gene, on the evaluation cells: the prediction calibrated on calib_fold c (delta from the held-out target's own gene on calib_fold c) is evaluated on calib_fold 1 - c, for both c.
- Predicted means average over all NTC reference cells of the condition (both calib folds), while the delta uses only calib_fold c cells; NTC reference cells carry no target response.
- Per target: Pearson r and RMSE per condition and direction; directions and then conditions combined on the Fisher z scale (atanh, |r| clipped at 0.999999999999); RMSE averaged.
- Amendment (2026-10-21, before model fitting): in the primary decision test and every paired test on z, a constant prediction gets r = 0 (z = 0), so all 56 targets enter; the pre-registered complete-pairs analysis is supporting. Wilcoxon zero differences are dropped (zero_method 'wilcox') and counted.
- Sign convention: differences are cell-specific minus reference; a positive z difference or a negative RMSE difference favours the cell-specific model. Hodges-Lehmann (HL) estimates include zero differences; CIs are percentile bootstraps over targets (5,000 draws, seed 20261011).

## Primary decision test (cell-specific vs shared CSCN-SEM, amended rule)

- n = 56, NaN pairs dropped 0, zero differences dropped 0, n_eff = 56, approx; W = 600.0, p = 0.106; HL -0.009 on the z scale (tanh -0.009; 95% bootstrap CI [-0.019, 0.002]).
- Decision at alpha 0.05: H0 (no difference in held-out Fisher z) is not rejected.
- Median back-transformed r over the paired targets: cell-specific -0.012, shared -0.000.

## Supporting analyses

### Coverage (share of targets with a non-constant prediction) and accuracy (complete pairs)

| Model | Non-constant DMSO | Non-constant Vem | Non-constant in both | NaN pooled z (share) | > 10% flag |
|---|---|---|---|---|---|
| Cell-specific CSCN-SEM | 56/56 | 54/56 | 54/56 | 2 (0.04) | no |
| Shared CSCN-SEM | 33/56 | 29/56 | 16/56 | 40 (0.71) | yes |
| Cluster-level SEM | 53/56 | 55/56 | 52/56 | 4 (0.07) | no |
| Local association SEM | 56/56 | 56/56 | 56/56 | 0 (0.00) | no |
| Unweighted propagation | 56/56 | 54/56 | 54/56 | 2 (0.04) | no |
| Zero downstream | 0/56 | 0/56 | 0/56 | r not computed | - |

- Complete pairs (pre-registered NaN rule, supporting): n = 56, NaN pairs dropped 40, zero differences dropped 0, n_eff = 16, exact; W = 13.0, p = 0.003; HL -0.029 on the z scale (tanh -0.029; 95% bootstrap CI [-0.049, -0.012]).
- Coverage and accuracy are separate questions: the amended primary test credits a model for predicting a non-constant change for more targets (constant predictions score r = 0), whereas the complete-pairs test compares accuracy only on targets where both models predict a non-constant change.
- NaN flag (> 10% of targets with undefined pooled z under the pre-registered rule, any model with r): FLAGGED.

### Per-condition and RMSE tests

- DMSO (amended z): n = 56, NaN pairs dropped 0, zero differences dropped 0, n_eff = 56, approx; W = 638.0, p = 0.192; HL -0.012 on the z scale (tanh -0.012; 95% bootstrap CI [-0.028, 0.005]).
- Vem (amended z): n = 56, NaN pairs dropped 0, zero differences dropped 2, n_eff = 54, approx; W = 633.0, p = 0.346; HL -0.004 on the z scale (tanh -0.004; 95% bootstrap CI [-0.013, 0.006]).
- RMSE (mean over conditions): n = 56, NaN pairs dropped 0, zero differences dropped 0, n_eff = 56, approx; W = 566.0, p = 0.058; HL -4.69e-06 (95% bootstrap CI [-1.51e-05, -3.43e-08]).

### Medians over targets

| Model | Condition | Median r (NaN excluded, n) | Median r (r = 0, n) | Median RMSE |
|---|---|---|---|---|
| Cell-specific CSCN-SEM | DMSO | -0.019 (56) | -0.019 (56) | 0.0820 |
| Shared CSCN-SEM | DMSO | 0.002 (33) | 0.000 (56) | 0.0820 |
| Cluster-level SEM | DMSO | -0.002 (53) | -0.000 (56) | 0.0821 |
| Local association SEM | DMSO | -0.006 (56) | -0.006 (56) | 0.0820 |
| Unweighted propagation | DMSO | -0.018 (56) | -0.018 (56) | 0.0821 |
| Zero downstream | DMSO | not computed | not computed | 0.0820 |
| Cell-specific CSCN-SEM | Vem | -0.014 (54) | -0.011 (56) | 0.0774 |
| Shared CSCN-SEM | Vem | -0.008 (29) | 0.000 (56) | 0.0774 |
| Cluster-level SEM | Vem | -0.001 (55) | -0.000 (56) | 0.0775 |
| Local association SEM | Vem | -0.023 (56) | -0.023 (56) | 0.0774 |
| Unweighted propagation | Vem | -0.012 (54) | -0.010 (56) | 0.0773 |
| Zero downstream | Vem | not computed | not computed | 0.0774 |
| Cell-specific CSCN-SEM | pooled | -0.012 (54) | -0.012 (56) | 0.0804 |
| Shared CSCN-SEM | pooled | -0.019 (16) | -0.000 (56) | 0.0804 |
| Cluster-level SEM | pooled | -0.007 (52) | -0.004 (56) | 0.0804 |
| Local association SEM | pooled | -0.011 (56) | -0.011 (56) | 0.0804 |
| Unweighted propagation | pooled | -0.009 (54) | -0.008 (56) | 0.0804 |
| Zero downstream | pooled | not computed | not computed | 0.0804 |

### Per-fold medians (pooled over conditions)

| Fold | Model | Condition | Median r (NaN excluded, n) | Median r (r = 0, n) | Median RMSE |
|---|---|---|---|---|---|
| 0 | Cell-specific CSCN-SEM | pooled | -0.033 (12) | -0.033 (12) | 0.0914 |
| 0 | Shared CSCN-SEM | pooled | -0.017 (3) | -0.010 (12) | 0.0916 |
| 0 | Cluster-level SEM | pooled | -0.001 (11) | 0.001 (12) | 0.0915 |
| 0 | Local association SEM | pooled | -0.026 (12) | -0.026 (12) | 0.0914 |
| 0 | Unweighted propagation | pooled | -0.034 (12) | -0.034 (12) | 0.0917 |
| 0 | Zero downstream | pooled | not computed | not computed | 0.0914 |
| 1 | Cell-specific CSCN-SEM | pooled | 0.009 (10) | 0.010 (11) | 0.0795 |
| 1 | Shared CSCN-SEM | pooled | -0.026 (1) | 0.000 (11) | 0.0795 |
| 1 | Cluster-level SEM | pooled | -0.002 (10) | -0.001 (11) | 0.0795 |
| 1 | Local association SEM | pooled | -0.008 (11) | -0.008 (11) | 0.0795 |
| 1 | Unweighted propagation | pooled | 0.007 (10) | 0.009 (11) | 0.0794 |
| 1 | Zero downstream | pooled | not computed | not computed | 0.0795 |
| 2 | Cell-specific CSCN-SEM | pooled | 0.006 (11) | 0.006 (11) | 0.0812 |
| 2 | Shared CSCN-SEM | pooled | 0.015 (5) | 0.000 (11) | 0.0811 |
| 2 | Cluster-level SEM | pooled | -0.014 (11) | -0.014 (11) | 0.0811 |
| 2 | Local association SEM | pooled | -0.012 (11) | -0.012 (11) | 0.0812 |
| 2 | Unweighted propagation | pooled | -0.002 (11) | -0.002 (11) | 0.0812 |
| 2 | Zero downstream | pooled | not computed | not computed | 0.0812 |
| 3 | Cell-specific CSCN-SEM | pooled | -0.008 (11) | -0.008 (11) | 0.0781 |
| 3 | Shared CSCN-SEM | pooled | -0.067 (2) | -0.003 (11) | 0.0781 |
| 3 | Cluster-level SEM | pooled | -0.001 (11) | -0.001 (11) | 0.0781 |
| 3 | Local association SEM | pooled | -0.011 (11) | -0.011 (11) | 0.0781 |
| 3 | Unweighted propagation | pooled | -0.005 (11) | -0.005 (11) | 0.0781 |
| 3 | Zero downstream | pooled | not computed | not computed | 0.0781 |
| 4 | Cell-specific CSCN-SEM | pooled | 0.001 (10) | -0.012 (11) | 0.0783 |
| 4 | Shared CSCN-SEM | pooled | 0.021 (5) | 0.007 (11) | 0.0783 |
| 4 | Cluster-level SEM | pooled | -0.027 (9) | -0.015 (11) | 0.0783 |
| 4 | Local association SEM | pooled | -0.004 (11) | -0.004 (11) | 0.0783 |
| 4 | Unweighted propagation | pooled | 0.004 (10) | -0.007 (11) | 0.0783 |
| 4 | Zero downstream | pooled | not computed | not computed | 0.0783 |

### By connectivity stratum (descriptive; connected = an edge in >= 10% of NTC cells in both conditions)

| Stratum | Model | Condition | Median r (NaN excluded, n) | Median r (r = 0, n) | Median RMSE |
|---|---|---|---|---|---|
| connected | Cell-specific CSCN-SEM | pooled | -0.027 (29) | -0.027 (29) | 0.0830 |
| connected | Shared CSCN-SEM | pooled | -0.010 (12) | -0.001 (29) | 0.0830 |
| connected | Cluster-level SEM | pooled | -0.013 (27) | -0.013 (29) | 0.0829 |
| connected | Local association SEM | pooled | -0.019 (29) | -0.019 (29) | 0.0829 |
| connected | Unweighted propagation | pooled | -0.028 (29) | -0.028 (29) | 0.0834 |
| connected | Zero downstream | pooled | not computed | not computed | 0.0829 |
| weak | Cell-specific CSCN-SEM | pooled | 0.004 (25) | 0.004 (27) | 0.0730 |
| weak | Shared CSCN-SEM | pooled | -0.028 (4) | 0.000 (27) | 0.0730 |
| weak | Cluster-level SEM | pooled | -0.001 (25) | -0.001 (27) | 0.0730 |
| weak | Local association SEM | pooled | -0.009 (27) | -0.009 (27) | 0.0730 |
| weak | Unweighted propagation | pooled | -0.002 (25) | -0.002 (27) | 0.0730 |
| weak | Zero downstream | pooled | not computed | not computed | 0.0730 |

### Split-half noise ceiling (context only; not part of any test)

- DMSO: median noise ceiling 0.224 (IQR 0.076 to 0.491); targets with a ceiling < 0: 5 of 56.
  - Descriptive, not part of any test: targets whose held-out r exceeds their noise ceiling, among targets with a non-constant prediction: Cell-specific CSCN-SEM 7 of 56; Shared CSCN-SEM 3 of 33.
- Vem: median noise ceiling 0.119 (IQR 0.002 to 0.222); targets with a ceiling < 0: 13 of 56.
  - Descriptive, not part of any test: targets whose held-out r exceeds their noise ceiling, among targets with a non-constant prediction: Cell-specific CSCN-SEM 12 of 54; Shared CSCN-SEM 5 of 29.

## Secondary comparisons (cell-specific vs each comparator; descriptive, Holm within each metric family)

| Comparator | Metric | Version | n_eff (zeros dropped) | Method | p | Holm p | HL [95% CI] |
|---|---|---|---|---|---|---|---|
| Cluster-level SEM | z | amended (descriptive) | 56 (0) | approx | 0.184 | 0.467 | -0.009 [-0.025, 0.005] (tanh) |
| Local association SEM | z | amended (descriptive) | 56 (0) | approx | 0.271 | 0.467 | 0.006 [-0.005, 0.016] (tanh) |
| Unweighted propagation | z | amended (descriptive) | 56 (0) | approx | 0.156 | 0.467 | -0.002 [-0.005, 0.001] (tanh) |
| Cluster-level SEM | z | complete_pairs (supporting) | 50 (0) | exact | 0.152 | 0.455 | -0.011 [-0.029, 0.004] (tanh) |
| Local association SEM | z | complete_pairs (supporting) | 54 (0) | approx | 0.236 | 0.455 | 0.007 [-0.005, 0.017] (tanh) |
| Unweighted propagation | z | complete_pairs (supporting) | 54 (0) | approx | 0.172 | 0.455 | -0.002 [-0.006, 0.001] (tanh) |
| Zero downstream | rmse | all (descriptive) | 56 (0) | approx | 0.065 | 0.196 | 1.77e-06 [-4.28e-08, 5.95e-06] |
| Cluster-level SEM | rmse | all (descriptive) | 56 (0) | approx | 0.562 | 0.562 | -1.98e-06 [-9.39e-06, 5.11e-06] |
| Local association SEM | rmse | all (descriptive) | 56 (0) | approx | 0.181 | 0.362 | 9.48e-07 [-2.48e-07, 4.93e-06] |
| Unweighted propagation | rmse | all (descriptive) | 56 (0) | approx | 0.002 | 0.007 | -1.98e-05 [-5.56e-05, -5.06e-06] |

## Secondary endpoint: Fufa score mean-effect error

Observed: Aim 1 per-cell scores, target - NTC on the evaluation fold; predicted: mean stored delta-scaled Fufa shift over included reference cells (genes outside the panel change by 0 for every model). Error |predicted - observed| per direction, averaged over directions.

| Contrast | Cell-specific CSCN-SEM | Shared CSCN-SEM | Cluster-level SEM | Local association SEM | Unweighted propagation | Zero downstream |
|---|---|---|---|---|---|---|
| KD_minus_NTC_DMSO | 0.098 | 0.098 | 0.098 | 0.098 | 0.100 | 0.098 |
| KD_minus_NTC_Vem | 0.118 | 0.119 | 0.119 | 0.119 | 0.114 | 0.119 |
| interaction_Vem_minus_DMSO | 0.156 | 0.156 | 0.156 | 0.156 | 0.157 | 0.156 |

- KD_minus_NTC_DMSO, cell-specific - shared absolute error: n = 56, NaN pairs dropped 0, zero differences dropped 8, n_eff = 48, approx; W = 524.0, p = 0.512; HL 2.02e-06 (95% bootstrap CI [-9.36e-06, 3.22e-05]) (no multiplicity correction pre-specified).
- KD_minus_NTC_Vem, cell-specific - shared absolute error: n = 56, NaN pairs dropped 0, zero differences dropped 7, n_eff = 49, approx; W = 384.0, p = 0.023; HL -1.38e-05 (95% bootstrap CI [-7.46e-05, -9.28e-09]) (no multiplicity correction pre-specified).
- interaction_Vem_minus_DMSO, cell-specific - shared absolute error: n = 56, NaN pairs dropped 0, zero differences dropped 1, n_eff = 55, approx; W = 493.0, p = 0.020; HL -2.58e-05 (95% bootstrap CI [-1.71e-04, -1.08e-06]) (no multiplicity correction pre-specified).

## Tertiary endpoints (expected low power)

Aim 1 found 0 primary dispersion calls and the Fufa score has low single-cell reliability, so these endpoints have low power (results/summary/power_summary.md). Predicted distribution: evaluation-fold NTC scores plus each cell's own stored shift (shared and zero-downstream models: a constant shift, so their predicted dispersion change is exactly 0); the predicted statistic compares the shifted with the same unshifted cells.

| Model | Median abs. error: log var ratio | IQR difference | Centred W1 | Over-dispersion rate | Coverage |
|---|---|---|---|---|---|
| Cell-specific CSCN-SEM | 0.124 | 0.140 | 0.109 | 0.080 | 0.915 |
| Shared CSCN-SEM | 0.124 | 0.140 | 0.110 | 0.080 | 0.915 |
| Cluster-level SEM | 0.124 | 0.140 | 0.110 | 0.080 | 0.915 |
| Local association SEM | 0.124 | 0.139 | 0.110 | 0.080 | 0.915 |
| Unweighted propagation | 0.124 | 0.139 | 0.108 | 0.080 | 0.915 |
| Zero downstream | 0.124 | 0.140 | 0.110 | 0.080 | 0.915 |

- Over-dispersion: predicted log variance ratio above the upper bound of the observed 95% bootstrap interval (evaluation cells), averaged over directions per target x condition row (112 rows); nominal <= 0.025. Aim 1 all-cell intervals are listed alongside in evaluation_distribution.csv (columns only).
- Evaluation NTC cells without a defined shift (dropped, summed over targets): Cell-specific CSCN-SEM DMSO c=0: 1736; Unweighted propagation DMSO c=0: 1736; Cell-specific CSCN-SEM DMSO c=1: 2184; Unweighted propagation DMSO c=1: 2184

## Post hoc description (not pre-registered)

Descriptive only; no test. Computed from the stored predictions and the observed change vectors over the 249 evaluation genes: for each model, target, condition and direction, the root mean square (RMS) of the predicted downstream change, relative to the knockdown size |delta| and to the RMS of the observed downstream change; medians over the 224 target x condition x direction units. For scale, the median |delta| is 0.167 and the median observed RMS is 0.080 (log1p CP10k); the median absolute observed Fufa KD - NTC shift is 0.103 (DMSO NTC SD).

| Model | Median predicted RMS | Ratio to abs(delta) | Ratio to observed RMS | Median abs. predicted Fufa KD - NTC shift |
|---|---|---|---|---|
| Cell-specific CSCN-SEM | 7.55e-05 | 4.02e-04 | 9.34e-04 | 4.05e-05 |
| Shared CSCN-SEM | 2.09e-04 | 1.78e-03 | 2.60e-03 | 0.00e+00 |
| Cluster-level SEM | 5.13e-04 | 2.87e-03 | 6.22e-03 | 2.41e-04 |
| Local association SEM | 3.27e-05 | 2.11e-04 | 3.94e-04 | 4.74e-05 |
| Unweighted propagation | 5.04e-04 | 2.91e-03 | 6.58e-03 | 3.36e-04 |
| Zero downstream | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |

- The predicted downstream changes of every model are much smaller than the observed downstream changes (median ratio to the observed RMS at most 0.7%), and the predicted Fufa shifts are much smaller than the observed ones. Each model's RMSE is therefore essentially the RMS of the observed change, which is the RMSE of the zero-downstream baseline, and each model's Fufa error is essentially the observed Fufa shift. This explains why every model is close to the zero baseline on RMSE and on the Fufa endpoint. Pearson r does not depend on the scale of the prediction, so the small magnitude does not by itself explain the near-zero correlations.

## Interpretation notes and limitations

- The cell-specific vs shared contrast is mainly structural: the NTC-only cross-validation selected lambda 10 (DMSO) and 10 (Vem), the largest grid value, and pure pooled coefficients (lambda = infinity, reference only) fitted marginally better, so local coefficients add little beyond per-cell graph structure.
- The shared and cluster graphs use a population Spearman partial-correlation PC search, whereas the cell-specific graphs use the local CSCN count test: the comparison differs in CI-test family as well as in locality. Coverage matching equalises the number of targets with an edge, not the number of edges.
- Weak target connectivity: connectivity follows NTC detection of the target gene; many targets have few or no downstream edges, and 40 targets are constant under the shared model in at least one condition.
- Low score reliability (Fufa split-half 0.01 / 0.12) limits the distribution-level endpoints.
- 70 DMSO cells (3.5%) are excluded from the cell-specific and ablation pseudobulk means because their stored graphs contain a cycle among fixed edges.
- Fallback completions (relaxed rule) may introduce v-structures not supported by the data; the local association model has no separating sets and no collider step.
