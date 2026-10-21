# CSCN-SEM knockdown predictions

Predicted pseudobulk expression changes (log1p CP10k) for every pre-registered model, held-out target,
condition and cross-fit direction (config/sem_prediction.json, frozen before the first prediction run).
This summary contains no evaluation metric: no observed response on evaluation cells, correlation or error
was computed. The evaluation follows config/evaluation_endpoints.json including its 2026-10-21 amendment
(constant predictions get r = 0 in the primary decision test; complete pairs are supporting).

## Graphs (structure only)
- DMSO shared graph: 181 edges, 35/56 targets with an incident edge (required 35), log10 alpha -4.7.
- Vem shared graph: 97 edges, 32/56 targets with an incident edge (required 32), log10 alpha -5.7.
- DMSO cluster graphs (K = 6): 66-138 edges, 30-37 targets with an edge.
- Vem cluster graphs (K = 6): 65-115 edges, 28-36 targets with an edge.
- Shared and cluster graphs use a population Spearman partial-correlation test inside the vendored PC-stable
  search (target-coverage matched); the cell-specific graphs use the local CSCN test. This difference in test
  family is a stated limitation.

## Orientation and local IDA (structure + training-fold seeds)
| graph_type | condition | max_U | fallback | no_choice |
|---|---|---|---|---|
| cell_specific | DMSO | 3 | 0.026 | 0.035 |
| cell_specific | Vem | 3 | 0.001 | 0.000 |
| cluster | DMSO | 1 | 0.000 | 0.000 |
| cluster | Vem | 1 | 0.000 | 0.000 |
| local_association | DMSO | 8 | 0.492 | 0.000 |
| local_association | Vem | 9 | 0.150 | 0.000 |
| shared | DMSO | 1 | 0.000 | 0.000 |
| shared | Vem | 1 | 0.000 | 0.000 |

- fallback: share of fold DAGs completed only under the relaxed rule (colliders may change); no_choice: share
  of (reference cell, held-out target) units without any valid completion.
- The remaining 3.5% of DMSO cell-specific units without a valid choice come from 70 of 1,989 DMSO cells whose
  stored post-collider PDAG already contains a directed cycle among its fixed (collider-oriented) edges, from
  conflicting collider orientations in the CSCN step. No completion can keep every fixed edge, so these cells
  are excluded from the cell-specific and unweighted-ablation pseudobulk means. This is below the 5% exclusion
  flag, so no rule was changed. Vem has no such cells.
- Seed proposals (eligible training-target responses, q < 0.05) per fold: DMSO 511, 522, 514, 547, 501; Vem 159, 147, 203, 192, 185.

## Shrinkage lambda (NTC-only cross-validation)
- DMSO: selected lambda 10.0 (grid boundary: True); MSE 0.49293 vs 0.49290 for lambda = infinity (reference only).
- Vem: selected lambda 10.0 (grid boundary: True); MSE 0.47981 vs 0.47972 for lambda = infinity (reference only).
- The NTC-only cross-validation selected the boundary value 10 in both conditions, and pure pooled coefficients
  (lambda = infinity, reference only, not selectable) fit marginally better: there is no detectable gain from
  local coefficients in control cells. The cell-specific model therefore differs from the shared model mainly
  through per-cell structure. Lambda and the grid were not changed.

## Predictions per model
| model | condition | n_constant_prediction | genes_nonzero_median | excluded_share | n_targets_exclusion_flag | mean_undirected_per_unit | fallback_share_of_valid_choices |
|---|---|---|---|---|---|---|---|
| cell_specific_cscn_sem | DMSO | 0 | 24 | 0.035 | 0 | 0.120 | 0.027 |
| shared_cscn_sem | DMSO | 23 | 1 | 0.000 | 0 | 0.214 | 0.000 |
| cluster_level_sem | DMSO | 3 | 4 | 0.000 | 0 | 0.264 | 0.000 |
| local_association_network_sem | DMSO | 0 | 180 | 0.000 | 0 | 0.206 | 0.497 |
| unweighted_propagation_ablation | DMSO | 0 | 24 | 0.035 | 0 | 0.120 | 0.027 |
| zero_downstream_baseline | DMSO | 56 | 0 | 0.000 | 0 | 0.000 | nan |
| cell_specific_cscn_sem | Vem | 2 | 27 | 0.000 | 0 | 0.113 | 0.001 |
| shared_cscn_sem | Vem | 27 | 1 | 0.000 | 0 | 0.286 | 0.000 |
| cluster_level_sem | Vem | 1 | 4 | 0.000 | 0 | 0.253 | 0.000 |
| local_association_network_sem | Vem | 0 | 238 | 0.000 | 0 | 0.394 | 0.167 |
| unweighted_propagation_ablation | Vem | 2 | 27 | 0.000 | 0 | 0.113 | 0.001 |
| zero_downstream_baseline | Vem | 56 | 0 | 0.000 | 0 | 0.000 | nan |

- n_constant_prediction: held-out targets whose predicted change is 0 for every gene other than the
  perturbed gene (undefined Pearson r; r = 0 under the amended primary rule). Targets constant in at least one
  condition (NaN under the supporting complete-pairs rule): cell_specific_cscn_sem 2, shared_cscn_sem 40, cluster_level_sem 4, local_association_network_sem 0, unweighted_propagation_ablation 2, zero_downstream_baseline 56.
- zero_downstream_baseline is constant by definition (RMSE only).
- Exclusion flag: (model, target, condition) rows with > 5% of reference cells excluded.
