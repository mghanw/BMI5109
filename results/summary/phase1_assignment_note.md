# Phase 1b: perturbation assignment rule

**Decision (Phase 1 review):** the authors' sgRNA rule is the primary assignment rule. The
presence-threshold rule is kept only as a sensitivity analysis. The authors' rule reproduces the
earlier report's singlet counts exactly. Source tables: `phase1_assignment_summary.csv`,
`phase1_report_comparison.csv` and `phase1_eligibility_by_rule.csv`.

## Rules
| Label | Role | Definition |
|---|---|---|
| `paper_gt5_max0.5_sec0.3` | primary | total sgRNA UMI > 5, max/total >= 0.5, second/max <= 0.3 |
| `paper_ge5_max0.5_sec0.3` | sensitivity | as above, with total sgRNA UMI >= 5 |
| `threshold_T5_auto` (= `threshold_T5`) | sensitivity | exactly one construct with >= T UMI; T = 5 chosen at the density valley |
| `threshold_T3`, `threshold_T10` | sensitivity | the same, with a fixed T |

## Deviation from the authors' code
The authors' code (`perturbfate_repo_bafdfff/Dual_omics_cell_calling_dual_omics_combine_multi_modal_cells_functions.R`,
line 160) filters `total_sgRNA_UMI >= 5`.
- GEO's sgRNA count matrix has no cells below 5 UMI.
- 2,962 cells (DMSO 1,460, Vem 1,502) have exactly 5 UMI. Under `> 5` they are all
  `unassigned_low`.
- Under `>= 5`, 97 of them would become singlets (DMSO 19, Vem 78). The rest would be multiplets.
- We keep `> 5` because only that version reproduces the report's counts. With `>= 5`, NTC
  would be 1,990 / 1,869 against the report's 1,989 / 1,864.

## Counts per condition
| Rule | Condition | Singlet | Multiplet | Unassigned (low) | No sgRNA |
|---|---|---:|---:|---:|---:|
| paper `> 5` (primary) | DMSO | 29,038 | 12,908 | 1,460 | 9,545 |
| | Vem | 24,175 | 9,427 | 1,502 | 12,535 |
| paper `>= 5` | DMSO | 29,057 | 14,349 | 0 | 9,545 |
| | Vem | 24,253 | 10,851 | 0 | 12,535 |
| threshold T = 5 | DMSO | 29,364 | 5,467 | 8,575 | 9,545 |
| | Vem | 24,696 | 3,648 | 6,760 | 12,535 |

The two rule families agree on the number of singlets to within about 2%, but they label
low-purity cells differently:
- The threshold rule labels a cell with one construct at >= 5 UMI and background reads below 5 a
  singlet, even when the cell's purity is low.
- The authors' rule calls such a cell a multiplet when max/total < 0.5 or second/max > 0.3.

`rna_singlets.h5ad` now holds 53,213 primary-rule singlets (previously 54,060).

## Eligible targets
Eligibility: > 50 singlets and adequate knockdown (remaining ratio <= 0.5, BH q < 0.05) in both
conditions.

| Rule | Eligible |
|---|---:|
| paper `> 5` (primary) | 56 |
| paper `>= 5` | 56 (same set as primary) |
| threshold T = 5 | 55 (same set as Phase 1) |

- **Gained under the primary rule: CCDC117.** Under the threshold rule, its DMSO remaining ratio
  was 0.508 (208 cells), just above the 0.5 cut, so its status was `inadequate`. Under the primary
  rule it is 0.496 (211 cells, q = 0.019). In Vem it passes under both rules (0.478 vs 0.447). It
  is a borderline target, and its eligibility depends on which cells are counted.
- **Lost:** none.

## Calibration folds
`calib_fold` depends only on the first-ligation well, the condition and the seed, not on the
assignment rule.
- Old vs new cell tables: 0 of 100,590 cells and 0 of 96 wells changed fold.
- Each condition keeps 24 wells per fold.

Primary-rule singlets per condition and fold:

| Condition | Fold 0 | Fold 1 |
|---|---:|---:|
| DMSO | 14,380 | 14,658 |
| Vem | 12,555 | 11,620 |
