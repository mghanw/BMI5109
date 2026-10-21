# Project summary: cell-specific CSCN-SEM models of CRISPRi perturbation responses

Status on 2026-10-21. Data preparation, Aim 1, the reliability and power analysis, the pre-registered
evaluation design and both parts of Aim 2 (network panel with cell-specific CSCN graphs; SEM knockdown
predictions) are complete. **The pre-registered held-out evaluation has not been run, and no held-out
performance metric (correlation, error or noise ceiling) has been computed.** Everything below is built
from the committed tables in `results/summary/` and the configuration files in `config/`.

## 1. Data and design
- **Data.** PerturbFate A375 melanoma CRISPRi screen, GEO GSE291147, Dual-omics RNA modality, under DMSO and
  vemurafenib (Vem).
- **Perturbation assignment.** The authors' rule is primary: total sgRNA UMI > 5, max/total >= 0.5 and
  second/max <= 0.3. Sensitivity rules are reported alongside it (`data_prep_assignment_summary.csv`,
  `data_prep_assignment_note.md`).
- **Cells.** The primary rule gives 53,213 singlets: 29,038 in DMSO and 24,175 in Vem.
  - Non-targeting control (NTC) cells: 1,989 (DMSO) and 1,864 (Vem).
  - Library targets: 143. Of these, 56 are eligible, meaning more than 50 singlets and adequate knockdown
    in both conditions (`data_prep_target_qc.csv`, `data_prep_eligibility_by_rule.csv`).
- **Held-out split.** A seeded 5-fold split over the 56 eligible targets (fold sizes 12/11/11/11/11),
  stratified only by knockdown-strength tertile. It was frozen before any response statistic was computed
  (`config/target_split.json`).
- **Calibration and evaluation cells.** Replicate labels cannot be recovered within the Dual-omics data,
  so the whole-well `calib_fold` split is used. A held-out target's knockdown size (delta) is calibrated
  on one fold and evaluated on the other, in both directions (`data_prep_replicate_evidence.md`).
- **Leakage rules.**
  - Structure, coefficients, gene selection and hyperparameters use training data only.
  - Delta uses only the held-out target's own perturbed gene on its calibration fold.
  - The perturbed gene is excluded from evaluation.
  - All models share features, calibration information and splits.

## 2. Aim 1: mean and distribution contrasts (unchanged since data-v1)
Full details are in `aim1_summary.md` and `power_summary.md`.

The primary readout is the Fufa 2019 AZD-up SOX10-bound pigmentation program score.
- **Drug effect.** Vemurafenib raises it in NTC cells by +1.06 [0.99, 1.13] DMSO-NTC SD
  (`aim1_mean_contrasts.csv`).
- **Mean shifts by knockdown are rare.** The largest gene × drug interactions are MED12 −0.76
  [−0.98, −0.53] and MED15 −0.72 [−1.03, −0.42].
- **No dispersion call.** No knockdown changed the dispersion of the score without a mean shift
  (0 primary dispersion calls).
- **This null is uninformative.** The score's split-half reliability in NTC cells is 0.01 (DMSO) and
  0.12 (Vem) (`score_reliability.csv`). At the median target size, 80% power after Bonferroni correction
  requires a log variance ratio of about 0.48 (DMSO) and 0.57 (Vem) (`power_simulation.csv`).

## 3. Pre-registered evaluation
- **Endpoints.** The held-out endpoints were pre-registered in `config/evaluation_endpoints.json`
  (frozen 2026-10-11) before any modelling code.
  - Primary endpoint: Pearson r across panel genes between predicted and observed downstream changes,
    together with RMSE.
  - Decision test: a single two-sided Wilcoxon signed-rank test over held-out targets on the paired
    Fisher z difference between cell-specific and shared CSCN-SEM.
  - Effect size: Hodges-Lehmann, with a bootstrap 95% CI.
  - Further comparisons: secondary comparisons with Holm correction, plus secondary and tertiary
    endpoints.
- **Amendment of 2026-10-21.** The file carries one dated `amendments` entry. It was committed in
  `81a3bd1` ("Freeze CSCN-SEM knockdown prediction model definitions"), together with the frozen model
  definitions, before any SEM coefficient, prediction or held-out response was computed.
  - It changes two rules:
    - Constant predictions get Pearson r = 0 (z = 0) in the primary decision test, so all 56 targets
      enter. The pre-registered complete-pairs analysis becomes supporting.
    - Wilcoxon zero differences use `zero_method='wilcox'`: pairs with a zero difference are dropped
      and their number is reported.
  - Reason: structure-only checks showed that under the shared model 40 of 56 targets have a constant
    prediction in at least one condition. Complete pairs would therefore leave 16 targets, selected by
    shared-graph connectivity (`sem_prediction_structure_targets.csv`).
  - Every other field is unchanged.

## 4. Aim 2, part 1: network panel and cell-specific CSCN graphs
Details are in `cscn_panel_summary.md`, with configuration in `config/cscn_panel.json` and
`config/cscn_networks.json`.
- **Panel (250 genes, NTC cells only).**
  - Forced genes: the 56 eligible targets and the 45 mapped Fufa genes (PTEN is both), 100 in total.
    17 Fufa genes have NTC detection below 5% in both conditions; they are flagged and kept.
  - Fill genes: 150, ranked by a Seurat-style mean.var.plot standardisation on log1p CP10k with
    equal-frequency mean bins. The universe is the 3,253 genes detected in at least 10% of NTC cells in
    each condition, excluding MT- and ribosomal genes (`cscn_panel_genes.csv`).
  - The panel and graphs are fold-independent, which is stricter than the pre-registered "training
    targets plus NTC" rule.
- **Hyperparameters.** These were chosen from graph diagnostics and runtime only, never from prediction
  performance.
  - Panel size 250 and maximum conditioning set 4 (`cscn_panel_timing.csv`).
  - CSCN window 0.35, fixed a priori.
  - CI alpha 0.10. Against 0.05 it roughly doubles the edges per cell (42 → 78) and raises target
    connectivity (0.18 → 0.29) (`cscn_panel_sensitivity.csv`).
- **Graphs.** One post-collider PDAG per NTC cell with the pinned official CSCN, with 0 CI-test errors
  (`cscn_panel_graph_summary.csv`).
  - Median edges per cell: 64 [IQR 49, 77] in DMSO and 61 [36, 81] in Vem; about 26% of edges are
    undirected.
  - Run time: 24.9 and 19.0 CPU hours.
- **Target connectivity is weak and tracks NTC detection.** The Spearman correlation between the share
  of cells with an edge (`cscn_panel_target_connectivity.csv`) and NTC detection
  (`cscn_panel_genes.csv`) is 0.74 (DMSO) and 0.80 (Vem).
  - 34 of 56 targets have no edge in more than half of the cells in both conditions.
  - 29 targets are "connected", meaning they have an edge in at least 10% of cells in both conditions.

## 5. Aim 2, part 2: SEM knockdown predictions
Details are in `sem_prediction_summary.md`; the definitions are frozen in `config/sem_prediction.json`.
- **Models.**
  1. `cell_specific_cscn_sem` (primary): per-cell DAGs with local coefficients.
  2. `shared_cscn_sem` (primary reference): one graph per condition with pooled coefficients.
  3. `cluster_level_sem`: K = 6 clusters of NTC cells, each with its own graph and pooled coefficients.
  4. `local_association_network_sem`: per cell, the most strongly correlated gene pairs within the
     cell's neighbourhood, density-matched to that cell's CSCN graph.
  5. `unweighted_propagation_ablation`: the model-1 DAGs with b = 1/|pa|.
  6. `zero_downstream_baseline`.
- **Orientation (fold-specific).** Seeds come from training-fold targets only: Welch tests of the
  target − NTC response, BH within fold and condition, q < 0.05. That gives 501–547 seed proposals per
  fold in DMSO and 147–203 in Vem. The pipeline is: unchanged seeding → official Meek rules → official
  DAG completion (`sem_prediction_orientation.csv`).
- **Local IDA for the held-out target.**
  - Enumerate the parent sets among its undirected neighbours that create no new unshielded collider.
  - The cell prediction is the equal-weight mean over the valid choices.
  - The limit is |U| <= 10; the largest observed |U| is 9 (`sem_prediction_structure.csv`).
- **Shared and cluster graphs.**
  - The CSCN test cannot use a single all-cell neighbourhood, because the graph is empty.
  - These graphs therefore use the vendored PC-stable search with a population Spearman
    partial-correlation test.
  - The threshold is target-coverage matched on structure only: the sparsest graph in which at least as
    many targets have an edge as have one in at least 10% of the same cells' CSCN graphs.
  - Result: 181 (DMSO) and 97 (Vem) edges in the shared graphs (`sem_prediction_shared_graphs.csv`).
- **Unified completion fallback.**
  - If the collider-preserving completion fails, a relaxed completion is accepted: acyclic, keeping
    every fixed edge.
  - This affected 2.6% of DMSO cell-specific fold DAGs and about 49% (DMSO) and 15% (Vem) of the local
    association fold DAGs.
- **Coefficients.** Ridge shrinkage towards pooled OLS on the 0.35·n nearest NTC cells.
  - The NTC-only cross-validation selected the boundary value λ = 10 in both conditions.
  - Pure pooled coefficients (λ = ∞, a reference row that cannot be selected) fit marginally better:
    MSE 0.49290 vs 0.49293 (DMSO) and 0.47972 vs 0.47981 (Vem) (`sem_prediction_lambda_cv.csv`).
  - So local coefficients give no detectable gain in control cells, and the cell-specific versus shared
    contrast is mainly structural.
- **Predictions.** Per-unit and delta-scaled pseudobulk changes for every model, target, condition and
  cross-fit direction, plus per-cell Fufa score shifts (outputs only, not committed).
  - Targets with a constant prediction in at least one condition (`sem_prediction_model_summary.csv`):

    | model | targets |
    |---|---|
    | cell-specific | 2 |
    | shared | 40 |
    | cluster | 4 |
    | local association | 0 |
    | unweighted ablation | 2 |
    | zero baseline | 56 (by definition) |
  - Exclusions: 3.5% of DMSO cell-specific units, from 70 of 1,989 DMSO cells whose stored PDAG already
    contains a directed cycle among its fixed edges. This is below the 5% exclusion flag; no (model,
    target, condition) row is flagged.
  - Runtime: about 1.6 h wall with 14 workers; 4.9 CPU hours for the cell-specific model and 16.2 for
    the local association model.

## 6. Limitations
- **Different test family.** The shared and cluster graphs use a population partial-correlation test,
  whereas the cell-specific graphs use the local CSCN count test. The primary contrast therefore differs
  in test family as well as in locality. Coverage matching equalises the number of targets with an edge,
  not the number of edges.
- **Weak target connectivity.** Connectivity follows NTC detection of the target gene. Many targets have
  few or no downstream edges, and 40 targets are constant under the shared model.
- **Low score reliability.** The Fufa score has low single-cell reliability (0.01 / 0.12), which limits
  the distribution-level (tertiary) endpoints.
- **Excluded cells.** 70 DMSO cells (3.5%) are excluded from the cell-specific and ablation pseudobulk
  means because their stored graphs contain a cycle among fixed edges.
- **No evaluation result yet.** No held-out evaluation result exists.

## 7. Next step
Run the pre-registered held-out evaluation in `config/evaluation_endpoints.json`, including its
2026-10-21 amendment, on the frozen predictions in `outputs/sem_prediction/predictions/`:
- observed target − NTC changes on the evaluation calibration fold;
- Pearson r and RMSE per target, model, condition and direction;
- the single primary Wilcoxon test (cell-specific vs shared), with the supporting and secondary analyses.

No held-out performance metric has been computed so far.
