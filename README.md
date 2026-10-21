# Cell-specific causal networks for perturbation response prediction

We asked whether cell-specific causal network models predict held-out CRISPRi knockdown responses better
than a shared network. Using the PerturbFate A375 melanoma screen (GEO GSE291147, Dual-omics RNA) under
DMSO and vemurafenib, we learned cell-specific causal networks (CSCN) from control cells, turned them into
linear structural equation models (CSCN-SEM), and predicted the downstream expression changes of
knockdown targets that were held out from all model fitting. The comparison with a shared-network
CSCN-SEM was pre-registered. The primary test did not detect a difference in held-out accuracy between
cell-specific and shared CSCN-SEM (two-sided Wilcoxon p = 0.106, 56 targets). Every model's
median held-out correlation with the observed changes is close to 0 (between −0.02 and 0), well below the
split-half noise ceiling (median 0.224 in DMSO and 0.119 in vemurafenib). A post hoc description shows that
every model predicts downstream changes far smaller than those observed (median at most 0.7% of the
observed magnitude). In the pre-registered complete-pairs analysis, restricted to the 16 targets for which
both models predict a change, the shared model is more accurate (p = 0.003).

## Team

BMI5109 Group 8, National University of Singapore: Yao Lihong, Wang Menghan, Lai Yingtong,
Chapala Dharani.

CSCN was developed by Wang Menghan and colleagues (Bioinformatics 2026;42(7):btag480). The copy in
`icpdag/_vendor/cscn/` is the unmodified official implementation; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Study design

- **Data.** PerturbFate A375 CRISPRi screen, Dual-omics RNA modality, under DMSO and vemurafenib (Vem).
  The authors' sgRNA assignment rule gives 53,213 singlets (29,038 DMSO, 24,175 Vem), including 1,989 and
  1,864 non-targeting control (NTC) cells. Of 143 library targets, 56 are eligible (more than 50 singlets
  and adequate knockdown in both conditions).
- **Held-out targets.** A seeded 5-fold split over the 56 eligible targets (fold sizes 12/11/11/11/11),
  stratified only by knockdown strength and frozen before any response statistic was computed
  ([config/target_split.json](config/target_split.json)). Network structure, coefficients, gene
  selection and settings use control cells and training-fold targets only.
- **Cross-fitting.** Replicate labels cannot be recovered within the Dual-omics data, so cells are split by
  whole well into two calibration folds. A held-out target's knockdown size (delta, its own gene only) is
  calibrated on one fold and its downstream response is evaluated on the other, in both directions.
- **Pre-registered endpoints.** [config/evaluation_endpoints.json](config/evaluation_endpoints.json):
  Pearson r (Fisher z) and RMSE across the 249 panel genes other than the perturbed gene; one primary
  decision test, a two-sided Wilcoxon signed-rank test over held-out targets on the paired z difference
  between cell-specific and shared CSCN-SEM, with a Hodges-Lehmann effect size and a bootstrap CI;
  secondary comparisons with Holm correction; a secondary Fufa score endpoint; tertiary
  distribution endpoints. One dated amendment (2026-10-21): a constant prediction counts as r = 0 in the
  primary test, so all 56 targets enter, and the complete-pairs analysis becomes supporting.
- **Pre-registration record.** The endpoints were frozen on 2026-10-11 before any modelling code; the
  amendment was committed in `81a3bd1` before any model fitting; the evaluation rules and code
  ([config/evaluation_run.json](config/evaluation_run.json)) were committed in `11985fa` before any
  held-out response was read.
- **Six pre-specified models** ([config/sem_prediction.json](config/sem_prediction.json)):
  1. cell-specific CSCN-SEM (primary): per-cell graphs with local coefficients;
  2. shared CSCN-SEM (primary reference): one graph per condition with pooled coefficients;
  3. cluster-level SEM: six clusters of control cells, each with its own graph;
  4. local association SEM: per-cell correlation graphs matched in density to the CSCN graphs;
  5. unweighted propagation: the cell-specific graphs with equal, unfitted coefficients;
  6. zero-downstream baseline: no downstream change.

## Results

### Aim 1: Few knockdowns change the Fufa score response

Primary readout: the Fufa 2019 SOX10-bound pigmentation program score, in DMSO NTC standard deviations.
Details: [aim1_summary.md](results/summary/aim1_summary.md) and
[power_summary.md](results/summary/power_summary.md).

- Vemurafenib raises the score in NTC cells by +1.06 [0.99, 1.13].
- A knockdown shifts the mean score (knockdown − NTC, BH q < 0.05) for 3 targets in DMSO and 4 in
  vemurafenib. Four targets have a significant gene × drug interaction; the largest is MED12, −0.76
  [−0.98, −0.53].
- No knockdown changes the dispersion of the score without a mean shift (0 primary dispersion calls; the
  smallest q is 0.13, MED12 in vemurafenib).
- This null does not show homogeneity: the score's split-half reliability in NTC cells is 0.01 (DMSO) and
  0.12 (vemurafenib), and at the median target size 80% power after Bonferroni correction needs a log
  variance ratio of about 0.48 (DMSO) and 0.57 (vemurafenib).
- The earlier mechanism report is reproduced: all 105 of its contrasts fall inside its 95% CIs.

Figures: [mean effects](reports/figures/aim1_forest_fufa.png),
[mean vs dispersion](reports/figures/aim1_mean_vs_dispersion.png),
[score distributions](reports/figures/aim1_ecdf_fufa.png), [power curves](reports/figures/power_curves.png).

### Aim 2: Cell-specific networks and knockdown predictions

Details: [cscn_panel_summary.md](results/summary/cscn_panel_summary.md) and
[sem_prediction_summary.md](results/summary/sem_prediction_summary.md).

- **Panel.** 250 genes chosen from NTC cells only: the 56 eligible targets and the 45 mapped Fufa genes
  (100 forced genes; PTEN is both), plus 150 fill genes ranked by standardised dispersion.
- **Graphs.** One CSCN graph per NTC cell, with a median of 64 (DMSO) and 61 (vemurafenib) edges per cell.
- **Target connectivity is weak.** The share of cells in which a target has an edge tracks the target's
  detection in NTC cells (Spearman 0.74 in DMSO, 0.80 in vemurafenib). 34 of 56 targets have no edge in
  more than half of the cells in both conditions; 29 targets are connected (an edge in at least 10% of
  cells in both conditions); 2 targets (ATP5G3, VBP1) have no edge in any vemurafenib cell.
- **Shared graphs.** 181 (DMSO) and 97 (vemurafenib) edges, matched so that as many targets have an edge
  (35 and 32 of 56) as have one in at least 10% of the cells' CSCN graphs. They use a population
  partial-correlation test, not the local CSCN test.
- **Coefficients.** Control-only cross-validation selected the largest shrinkage towards pooled
  coefficients (lambda = 10, the grid boundary) in both conditions, and fully pooled coefficients fit
  marginally better (MSE 0.49290 vs 0.49293 in DMSO, 0.47972 vs 0.47981 in vemurafenib). The
  cell-specific versus shared contrast is therefore mainly structural.
- **Constant predictions.** Targets whose predicted downstream change is zero for every gene in at least
  one condition: cell-specific 2, shared 40, cluster-level 4, local association 0, unweighted propagation 2
  of 56 targets.

Figure: [graph overview](reports/figures/cscn_panel_graph_overview.png).

### Aim 3: Held-out evaluation

Details: [evaluation_summary.md](results/summary/evaluation_summary.md); all evaluation tables are
`results/summary/evaluation_*.csv`.

- **Primary decision test (as pre-registered, with the amendment).** Cell-specific − shared paired Fisher z
  over 56 targets: Wilcoxon p = 0.106; Hodges-Lehmann −0.009 (back-transformed), 95% bootstrap CI
  [−0.019, 0.002]. The null hypothesis of no difference is not rejected at alpha 0.05.
- **Complete pairs (pre-registered rule, supporting).** On the 16 targets for which both models predict a
  non-constant change, the shared model is more accurate: p = 0.003, Hodges-Lehmann −0.029 [−0.049, −0.012].
- **Coverage versus accuracy.** The cell-specific model predicts a non-constant change for 54 of 56
  targets in both conditions, the shared model for 16 of 56. The amended primary test therefore mixes
  coverage (constant predictions count as r = 0) with accuracy; the complete-pairs test measures accuracy
  only.
- **Accuracy is low for every model.** Median held-out r over targets lies between −0.02 and 0 for every
  model, against a split-half noise ceiling with median 0.224 (DMSO) and 0.119 (vemurafenib). Per-condition
  tests (p = 0.192 in DMSO, 0.346 in vemurafenib) and RMSE (p = 0.058, Hodges-Lehmann −4.69e-06) show no
  meaningful difference.
- **Post hoc description (not pre-registered).** The predicted downstream changes are far smaller than the
  observed ones (median at most 0.7% of the observed root mean square for every model), so every model
  is close to the zero-downstream baseline on RMSE and on the Fufa score endpoint.
- **Secondary and tertiary endpoints (descriptive).** Against the cluster-level, local association and
  unweighted models, the paired z differences of the cell-specific model lie between −0.009 and 0.006
  (back-transformed Hodges-Lehmann; Holm p = 0.467 for each). RMSE differences are of order 1e-6 to 1e-5;
  the only one with Holm p < 0.05 is against unweighted propagation (−1.98e-05, Holm p = 0.007). Fufa
  score errors of the cell-specific and shared models differ by at most order 1e-5 (Hodges-Lehmann). The
  distribution endpoints are almost identical across models (over-dispersion rate 0.080 and coverage 0.915
  for every model).
- **Limitations.** The shared and cluster graphs use a different conditional-independence test family
  from the cell-specific graphs; target connectivity is weak; the Fufa score has low single-cell
  reliability; and 70 DMSO cells (3.5%) are excluded from the cell-specific and unweighted-propagation
  pseudobulk means because their stored graphs contain a cycle among fixed edges.

Figures: [paired Fisher z](reports/figures/evaluation_paired_z.png),
[r vs noise ceiling](reports/figures/evaluation_r_vs_noise_ceiling.png).

An overview of all analyses is in [project_summary.md](results/summary/project_summary.md).

## Repository structure

| Path | Contents |
|---|---|
| `icpdag/` | Python package: vendored CSCN, network panel and cell-specific graphs, SEM orientation, coefficients and knockdown simulation, held-out evaluation statistics, Aim 1 scores and contrasts |
| `scripts/` | Command-line steps of the analysis, from data download to the held-out evaluation |
| `config/` | Frozen target split, program gene sets, network panel, model definitions, pre-registered endpoints and evaluation run rules |
| `results/summary/` | Aggregate result tables and summaries (no cell-level data) |
| `reports/` | Figures and the earlier mechanism report |
| `data/` | Input data; see [data/README.md](data/README.md) |
| `docs/` | [Reproduction guide](docs/reproduce.md) |
| `tests/` | `pytest` test suite |

## Reproducing the analysis

The environment and the full ordered sequence of commands, from the GEO files to the held-out
evaluation, are in [docs/reproduce.md](docs/reproduce.md).

## Data source and citation

Xu Z, Lu Z, Ugurbil A, et al. Mapping convergent regulators of melanoma drug resistance by
PerturbFate. *Nature* 2026;654(8117):261-71.
[doi:10.1038/s41586-026-10367-0](https://www.nature.com/articles/s41586-026-10367-0).
Data: GEO [GSE291147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291147); see
[data/README.md](data/README.md).

Third-party code: [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Earlier analysis

The [mechanism report](reports/mechanism-report.zh.md) (in Chinese) describes an earlier analysis of the
same A375 RNA data that combined perturbation information, cell-specific networks and average treatment
effects. It found that NF2, MED12 and MED13 knockdowns share enhanced YAP-associated outputs without drug,
whereas MED12 shows a weaker drug-induced differentiation response. Its summary tables (average effects,
CKM state and factor annotations, network directionality counts) are in `results/summary/`.
