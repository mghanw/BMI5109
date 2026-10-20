# Cell-specific CSCN graphs on the network panel

Structure is learned from NTC cells only (both conditions), so the graphs are fold-independent and one
graph set serves all five held-out folds. Stored graph per cell: the PDAG after unshielded-collider
orientation; intervention-based orientation is fold-specific and done when SEMs are fitted.

## Panel
- 250 genes: 100 forced (56 eligible targets + Fufa program genes; PTEN is
  both) and 150 fill genes from a universe of 3,253 genes detected in >= 10% of NTC cells in each condition.
- Fill ranking: Seurat-style mean.var.plot standardization, computed on log1p CP10k with equal-frequency
  mean bins (mean of the two within-condition ranks).
- 17 forced genes (all Fufa) have NTC detection < 5% in both conditions; they are kept and flagged (cscn_panel_genes.csv).

## Hyperparameters (graph diagnostics and runtime only, never prediction performance)
- window 0.35, CI alpha 0.1, max conditioning 4, stable PC, bitmap queries (official CSCN, pinned).
- Window kept at the a priori 0.35: window changes (0.5, 0.65) gave no consistent benefit, and 0.65 disconnected the most-detected targets. ci_alpha 0.10 roughly doubles edges per cell (median 42 -> 78) and raises target connectivity (mean share of target nodes with >= 1 edge 0.18 -> 0.29), while edges stay far below the number expected if the local CI test were calibrated, i.e. the test is conservative on zero-inflated counts. max_conditioning 2-4 changed runtime and edges by < 3% (sparse skeletons); 4 is used.
- Grids: cscn_panel_timing.csv (12 configurations of panel size x max conditioning) and
  cscn_panel_sensitivity.csv (6 configurations of window x CI alpha at P 250).

## Graphs
| condition | cells | CI errors | edges per cell, median [IQR] | undirected share | target nodes with an edge | unseeded DAG completion valid | CPU hours |
|---|---|---|---|---|---|---|---|
| DMSO | 1,989 | 0 | 64 [49, 77] | 0.26 | 0.27 | 0.939 | 24.9 |
| Vem | 1,864 | 0 | 61 [36, 81] | 0.26 | 0.27 | 0.999 | 19.0 |

- Edge frequency: 1,582 (condition, gene pair) rows present in >= 1% of cells (cscn_panel_edge_frequency.csv), for the shared and cluster-level comparators.

## Target connectivity
- Targets isolated in a majority of cells (no edge in > 50% of cells): DMSO 41, Vem 42 of 56; in both conditions 34 of 56.
- Connectivity tracks NTC detection of the target: in a cell where the target is not detected, its
  neighbourhood is the zero tie block and the local CI tests return independence.
- Targets with no edge in any cell: DMSO: none; Vem: ATP5G3, VBP1.
- Only targets without an edge in any cell get constant downstream predictions (undefined Pearson r).
  Under the pre-registered rule (NaN if either condition is undefined) the expected NaN share is 2/56 (3.6%),
  below the 10% flag.
- Weakly connected (an edge in < 1% of cells in at least one condition; 17 targets): AHNAK (DMSO 5.3% / Vem 0.9%), ARRDC3 (DMSO 0.3% / Vem 1.1%), ATP5G3 (DMSO 19.0% / Vem 0.0%), CCDC117 (DMSO 0.1% / Vem 8.7%), CSK (DMSO 15.0% / Vem 0.6%), DGKD (DMSO 0.4% / Vem 4.3%), DPH2 (DMSO 1.1% / Vem 0.3%), FHOD1 (DMSO 0.6% / Vem 20.8%), HPRT1 (DMSO 8.9% / Vem 0.4%), KEAP1 (DMSO 5.0% / Vem 0.2%), MED10 (DMSO 1.1% / Vem 0.1%), NRP1 (DMSO 73.0% / Vem 0.7%), PGRMC2 (DMSO 0.1% / Vem 42.6%), RNF7 (DMSO 1.0% / Vem 5.2%), SBF1 (DMSO 0.9% / Vem 0.4%), SUPT7L (DMSO 0.5% / Vem 4.1%), VBP1 (DMSO 3.7% / Vem 0.0%).
  Weakly connected targets get predictions from very few cells.
- Isolated in a majority of cells in both conditions: AHNAK, ARRDC3, ASAH1, ATP5G3, CAND1, CCDC117, CLTA, CSK, CUL5, DGKD, DPH2, DUSP6, FHOD1, GFM1, HPRT1, KCTD10, KCTD5, KEAP1, LATS2, MED10, MED12, MED15, PDCD10, PGRMC2, PLXND1, RNF7, SBF1, SUPT20H, SUPT7L, TADA2B, TADA3, TP53, VBP1, ZNF592.
