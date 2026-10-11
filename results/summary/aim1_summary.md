# Aim 1: mean and distribution contrasts of program scores

Run `outputs/aim1` (scripts/aim1.py, config/aim1.json): 53,213 primary-rule singlets,
56 eligible targets, both conditions, 5,000 stratified cell-bootstrap draws, 2,000 permutations.
All tables carry the target's `fold` from `config/target_split.json`.

**Leakage rule:** these descriptive results may be reported for all targets. Later analyses may use
them to choose panel genes, hyperparameters or models only through training folds.

## Held-out split (frozen before any score was computed)
- Seeded 5-fold split over the 56 eligible targets (fold sizes 12/11/11/11/11), stratified by
  knockdown-strength tertile only.
- The commit 'Freeze the held-out target split' (2026-10-11 12:39 +08:00) precedes the first
  scoring commit. `tests/test_target_split.py` checks that the split is deterministic, disjoint
  and complete, and that it regenerates exactly.

## Program gene sets
| Program | Role | Source genes | Mapped to GENCODE v38 |
|---|---|---:|---:|
| Fufa 2019 AZD-up SOX10-bound pigmentation (Table S2) | primary | 45 | 45 (44 exact, ADRBK2→GRK3) |
| Kanai 2024 YAP/TAZ up (Table S3, tab 6) | secondary | 80 | 80 |
| Kanai 2024 YAP/TAZ down | exploratory | 52 | 52 |
| Tsoi 2018 undifferentiated / neural-crest-like / transitory / melanocytic (Table S3) | secondary | 118 / 66 / 29 / 62 | all |

- **Fufa mapping:** the earlier report mapped 44 of 45 Fufa genes. The difference is GRK3, which
  is now resolved through its HGNC previous symbol.
- **Tsoi definition:** each Tsoi program uses only the genes labelled with that subtype. Labels
  shared by adjacent subtypes are not used. This reproduces the earlier report's Tsoi numbers;
  the inclusive definition did not (|z| up to 2.7).
- **Download sources:** PMC serves supplements only behind a browser proof-of-work check. The
  identical publisher / Europe PMC copies were fetched instead (`data/ref/SOURCES.md`, checksums
  in `config/programs.json`).
- **Scoring:** `icpdag/program_scores.py` was used unchanged, with an expression-matched
  background on a frozen DMSO NTC scale. All 143 library target genes are excluded from the
  background controls.
- **Perturbed gene inside a program:** only one pair, PTEN in Fufa. Dropping PTEN from the
  program changes the PTEN contrasts by less than 0.1 SD (e.g. DMSO mean −0.19 → −0.10, CI
  includes 0 in both cases); see `aim1_target_in_program_sensitivity.csv`.

## Main results (Fufa, primary readout)
- **Mean shifts are rare.**
  - Knockdown − NTC is significant (BH q < 0.05) for 3 targets in DMSO and 4 in Vem.
  - The gene × drug interaction is significant for 4: MED12 −0.76 [−0.98, −0.53], MED15 −0.72,
    PTPN14 −0.31, KEAP1 −0.30.
  - The drug response (Vem − DMSO) is positive in all 57 groups.
- **No pre-specified dispersion call.** No target reaches log variance ratio q < 0.05 in either
  condition, so there are **0 primary highlight calls**.
  - The smallest q is MED12 in Vem: log variance ratio −0.29 [−0.51, −0.11], p = 0.002,
    q = 0.13. This is a narrower distribution together with a large mean drop (−0.67), not a
    distribution change without a mean shift.
  - IQR difference and centered W1 are also null after BH.
  - W1 and KS reach q < 0.05 only for targets with clear mean shifts (MED12, MED15, SUPT20H,
    CAND1, KEAP1). These statistics react to location as well as shape.
- **Vemurafenib widens the Fufa distribution.**
  - NTC log variance ratio Vem vs DMSO is +0.38 [0.29, 0.47], and 29 of 57 groups are
    significant.
  - No target changes that drug effect on dispersion (0 significant interactions for log variance
    ratio or IQR).

## Secondary and exploratory readouts
- **Kanai up and Tsoi undifferentiated (DMSO):** many targets have significant log variance
  ratios (22 and 25). Almost all of these targets also have a significant mean shift (20 and 23).
  Across targets, the mean shift and the log variance ratio correlate strongly (r = +0.78 and
  +0.81): dispersion grows with the mean, as expected for scores built from sparse counts.
- **Exploratory OR rule:** 15 (target, readout, condition) rows meet "mean CI includes 0 and (log
  variance ratio, IQR or centered W1 q < 0.05)". All are in DMSO, in Kanai or Tsoi readouts, and
  11 of them survive the depth adjustment. Examples:
  - Tsoi neural-crest-like: TADA2B, KEAP1, PTPN14.
  - Tsoi undifferentiated: CDH13, CCDC117.
  - Kanai up: CST3.
  
  These are exploratory leads only (`aim1_highlight_targets.csv`).

## Technical sensitivity
- **Depth.** The log total UMI explains little of each score in DMSO NTC: R² 0.002 for Fufa and
  at most 0.07 (Tsoi undifferentiated).
  - After residualizing on depth (fit on DMSO NTC only), no Fufa log variance ratio, IQR or
    centered W1 result changes significance.
  - In the secondary readouts, 14 results cross q = 0.05 in either direction
    (`aim1_dispersion_depth_sensitivity.csv`).
- **Well clustering** (Fufa; ligation-1 well bootstrap, 48 wells per condition, shared draws).
  - Cluster intervals have about the same width as cell-bootstrap intervals (median width ratio
    0.98–0.99). Whether the CI excludes 0 agrees for 278 of 280 contrasts; the exceptions are
    MED10 Vem mean (cell CI excludes 0, cluster CI does not) and KIRREL Vem log variance ratio
    (the reverse).
  - The MED12 interaction keeps a CI that excludes 0: cluster [−0.98, −0.52].
  - There were no primary highlight calls to re-check.

## Expression contrasts
- **Genes tested:** 7,963 genes detected in ≥ 5% of NTC cells in either condition. The perturbed
  gene is excluded. The full matrix is `outputs/aim1/expression_contrasts.csv.gz` (not
  committed).
- **Genes changed per target** (|Δ mean log1p CP10k| ≥ 0.1 and BH q < 0.05):
  - Median 33 (DMSO) and 10.5 (Vem). Four targets per condition have none.
  - The largest are NF2 DMSO (614), MED12 DMSO (578), MED12 Vem (555) and TP53 DMSO (493).
- **Tables:** `aim1_expression_summary.csv` and `aim1_expression_top_genes.csv`.

## Reproduction of the earlier report (7 targets)
- **Cells:** the cell sets are identical in every contrast.
- **Fit:** all 105 contrasts (gene, drug, interaction and pair contrasts for Fufa, Kanai up and
  Tsoi undifferentiated) fall inside the report's 95% CI. The largest gap is |z| = 1.86 (in units
  of the report's SE).
  - MED12 Fufa interaction: −0.76 [−0.98, −0.53] vs −0.835 [−1.06, −0.61].
  - Drug response, NTC: +1.06 vs +1.13. MED12: +0.30 vs +0.29. NF2: +0.82 vs +0.85. MED13: +0.82
    vs +0.84.
  - MED12 − NF2 drug-response difference: −0.51 vs −0.55. MED12 − MED13: −0.51 vs −0.55.
- **Why they differ:** the remaining gaps come from which background genes are drawn. The earlier
  run excluded its unsaved 1,998 NMF network genes from the background pool; this run excludes the
  library target genes instead.
  - Reseeding the background draw 20 times moves the Fufa contrasts by SD 0.008–0.035.
  - The NTC drug response ranges 1.01–1.14 across seeds, which covers the report's 1.13.
  - 95 of 105 report values lie inside the reseeded range. The other 10 lie just outside it, with
    gaps of 0.02–0.10 (at most 0.82 report SE). The largest are the MED12 − NF2 DMSO pair
    contrast for Tsoi undifferentiated (0.10) and the MED12 Fufa interaction (0.08). For these the
    excluded background pool, not the random draw, is the likely source.

## Caveats
- **Interval scope:** cell-bootstrap intervals describe this cell pool. The well cluster bootstrap
  suggests little extra between-well variance, but it is not independent screen replication. The
  Tri-omics replicate is reserved for external validation.
- **Interactions for W1 and KS:** none are defined. Drug-effect W1 is descriptive only.
- **Figures:** `reports/figures/aim1_forest_fufa`, `aim1_mean_vs_dispersion`,
  `aim1_ecdf_fufa`. The ECDF figure shows the three smallest-q targets per condition,
  because there were no primary calls.
