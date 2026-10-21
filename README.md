# Cell-specific causal networks for perturbation response prediction

This project extends cell-specific causal network (CSCN) analysis to a structural equation model
(CSCN-SEM) and asks whether cell-specific causal network models predict held-out perturbation
responses better than shared-network models. We use the PerturbFate A375 melanoma CRISPRi screen
(GEO GSE291147, Dual-omics RNA) under DMSO and vemurafenib. Held-out knockdown targets will be
predicted from networks fitted on training targets only, and the predictions will be compared with the
observed downstream expression changes and with the mean and distribution of a fixed
differentiation program score.

## Team

BMI5109 Group 8, National University of Singapore: Yao Lihong, Wang Menghan, Lai Yingtong,
Chapala Dharani.

CSCN was developed by Wang Menghan and colleagues (Bioinformatics 2026;42(7):btag480). The copy in
`icpdag/_vendor/cscn/` is the unmodified official implementation; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Status

- Data preparation and Aim 1 (mean and dispersion contrasts of programme scores) are complete,
  together with a score reliability and power analysis.
- The held-out evaluation endpoints are pre-registered in
  [config/evaluation_endpoints.json](config/evaluation_endpoints.json) before any modelling code,
  with one dated amendment (2026-10-21) made before any model fitting.
- Aim 2 is complete: a 250-gene network panel, cell-specific CSCN graphs for all control cells, and
  knockdown predictions of six pre-specified SEM models for every held-out target, condition and
  cross-fit direction.
- Next: the pre-registered held-out evaluation. No held-out performance metric has been computed yet.
- Overview of all completed analyses: [results/summary/project_summary.md](results/summary/project_summary.md).

## Repository layout

| Path | Contents |
|---|---|
| `icpdag/` | Python package: vendored CSCN, intervention-assisted orientation, gene-space CKM, average effects, program scores, dispersion, expression contrasts, target split, power; network panel (`cscn_panel`), cell-specific CSCN graphs (`cscn_networks`), SEM inputs, comparator graphs, fold-specific orientation with local IDA, and SEM coefficients and knockdown simulation (`sem_inputs`, `sem_graphs`, `sem_orientation`, `sem_prediction`) |
| `scripts/` | Data download and conversion, data preparation, target split, program gene sets, Aim 1, reliability and power analysis, figures; network panel, CSCN timing and graph inference, SEM structure, lambda cross-validation, knockdown prediction and their summaries |
| `tests/` | `pytest` test suite |
| `config/` | Run configurations, frozen held-out target split, program gene sets, pre-registered held-out evaluation endpoints (with amendment), frozen network panel (`cscn_panel.json`), CSCN settings (`cscn_networks.json`), frozen SEM model definitions (`sem_prediction.json`) |
| `results/summary/` | Aggregate result tables and summaries (no cell-level data) |
| `reports/` | Figures and the earlier mechanism report |
| `data/` | Input data; see [data/README.md](data/README.md) |

## Environment

- Python 3.11 in a conda environment, then from the repository root:

  ```bash
  pip install -e ".[plots,io]" pytest
  ```

  The `plots` extra installs matplotlib; the `io` extra installs anndata and openpyxl.
- R >= 4.4 with the `Matrix` package, for converting the GEO RDS files.
- The full pipeline runs on a 32 GB RAM machine; count matrices are kept sparse.

## Reproduce

Run from the repository root (Windows paths shown).

1. Input data. The two small GEO files are in `data/raw/`. Get the RNA count matrix from the
   [data-v1 release](https://github.com/mghanw/BMI5109/releases/tag/data-v1) or download all three
   files from GEO (details and checksums in [data/README.md](data/README.md)):

   ```bash
   powershell -ExecutionPolicy Bypass -File scripts\download_data.ps1
   ```

2. Convert the RDS matrices to MatrixMarket:

   ```bash
   Rscript scripts\rds_to_mtx.R data\raw\GSE291147_Dual_omics_RNA_gene_count_matrix.RDS data\processed\rna_mtx
   ```

   ```bash
   Rscript scripts\rds_to_mtx.R data\raw\GSE291147_Dual_omics_sgRNA_count_matrix.RDS data\processed\sgrna_mtx
   ```

3. Build cells x features AnnData files:

   ```bash
   python scripts\mtx_to_h5ad.py data\processed\rna_mtx data\raw\GSE291147_Dual_omics_RNA_cell_metadata.csv.gz data\processed\rna.h5ad
   ```

   ```bash
   python scripts\mtx_to_h5ad.py data\processed\sgrna_mtx data\raw\GSE291147_Dual_omics_RNA_cell_metadata.csv.gz data\processed\sgrna.h5ad
   ```

4. Fetch the public reference files listed in [data/ref/SOURCES.md](data/ref/SOURCES.md) and check
   the GENCODE gene IDs against the matrix (this rewrites `SOURCES.md`; only access dates should
   change):

   ```bash
   python -I scripts\fetch_refs.py --download --gtf-check data\processed\rna_mtx\rows.txt
   ```

5. Data preparation and QC (perturbation assignment, knockdown QC, calibration folds;
   thresholds in `config/prepare_data.json`; smoke test with `--prefix <library prefix> --run-name NAME`):

   ```bash
   python scripts\prepare_data.py
   ```

6. Held-out target split. It is frozen; rerunning it must reproduce the committed file exactly:

   ```bash
   python scripts\make_target_split.py
   ```

   ```bash
   git diff --exit-code config/target_split.json
   ```

7. Program gene sets (Fufa 2019, Kanai 2024, Tsoi 2018 -> `config/programs.json`):

   ```bash
   python -I scripts\build_programs.py
   ```

8. Aim 1 (about 20 minutes; config `config/aim1.json`; smoke test with
   `--targets A,B --bootstrap 200 --run-name smoke_aim1`), then its figures:

   ```bash
   python scripts\aim1.py
   ```

   ```bash
   python scripts\aim1_figures.py
   ```

9. Reliability and power analysis (NTC cells only, about 5 minutes; config `config/power_analysis.json`;
   needs `outputs/aim1/scores.csv.gz` from step 8), then its figure:

   ```bash
   python scripts\power_analysis.py
   ```

   ```bash
   python scripts\power_figures.py
   ```

10. Network panel from NTC cells (about 1 minute; frozen in `config/cscn_panel.json`, so a rerun only
    verifies it), and the CSCN timing and sensitivity grids used to choose the settings (about 20 and
    15 minutes):

    ```bash
    python scripts\build_cscn_panel.py
    ```

    ```bash
    python scripts\cscn_panel_timing.py
    ```

    ```bash
    python scripts\cscn_panel_timing.py --sizes 250 --orders 2 --windows 0.35,0.5,0.65 --alphas 0.05,0.1
    ```

11. Cell-specific CSCN graphs for all NTC cells (about 1.8 h for DMSO and 1.4 h for vemurafenib with 14
    workers on a 16-thread machine; resumable, and shardable with `--shard i/n`), the merge per
    condition, then the summaries and figure:

    ```bash
    python scripts\infer_cscn_networks.py --condition DMSO
    ```

    ```bash
    python scripts\infer_cscn_networks.py --condition DMSO --merge
    ```

    Repeat both commands with `--condition Vem`, then:

    ```bash
    python scripts\summarize_cscn_networks.py
    ```

12. SEM structure (NTC only): shared and cluster graphs, PCs and neighbourhoods (about 11 minutes),
    then the fold-specific orientation diagnostics (about 1.7 h; cached per graph type):

    ```bash
    python scripts\build_sem_graphs.py
    ```

    ```bash
    python scripts\sem_structure_report.py
    ```

13. Shrinkage lambda by NTC-only cross-validation (about 8 minutes):

    ```bash
    python scripts\sem_lambda_cv.py
    ```

14. Knockdown predictions for each condition and outer fold (about 1.6 h for all 10 runs with 14
    workers; resumable and shardable), each followed by its merge, then the summaries:

    ```bash
    python scripts\predict_sem_knockdown.py --condition DMSO --fold 0
    ```

    ```bash
    python scripts\predict_sem_knockdown.py --condition DMSO --fold 0 --merge
    ```

    Repeat for folds 0-4 and both conditions, then:

    ```bash
    python scripts\summarize_sem_prediction.py
    ```

15. Tests:

    ```bash
    pytest
    ```

Long runs write to `outputs/<run_name>/` with a `config.json`; aggregate tables go to
`results/summary/` and figures to `reports/figures/`.

## Key results

From [results/summary](results/summary):

- 53,213 primary-rule singlets (authors' sgRNA assignment rule) and 56 eligible knockdown targets,
  split into a frozen 5-fold held-out design. Replicate labels cannot be recovered within the
  Dual-omics data, so whole-well calibration folds are used for calibration and evaluation.
- Mean shifts of the primary readout (Fufa 2019 pigmentation program) are rare: 3 targets in DMSO
  and 4 in vemurafenib (BH q < 0.05). Four targets have a significant gene x drug interaction, led
  by MED12 (-0.76 [-0.98, -0.53]).
- No target changes the dispersion of the Fufa score without a mean shift (0 primary dispersion
  calls; smallest q = 0.13, MED12 in vemurafenib).
- This null is uninformative about heterogeneity rather than evidence of homogeneity: the Fufa score's split-half reliability in NTC cells
  is 0.01 (DMSO) / 0.12 (vemurafenib), and at the median target size 80% power after Bonferroni
  correction requires a log variance ratio of about 0.48 / 0.57.
- The earlier mechanism report is reproduced: all 105 of its contrasts fall inside its 95% CIs.
- Aim 2: cell-specific CSCN graphs have a median of 64 (DMSO) and 61 (vemurafenib) edges per cell
  on the 250-gene panel. Target connectivity is weak and tracks the target's detection in control
  cells. Knockdown predictions of all six pre-specified models are complete. The control-only
  cross-validation found no gain from local coefficients, so the cell-specific versus shared contrast
  is mainly structural. See [project_summary.md](results/summary/project_summary.md); held-out
  performance has not been evaluated yet.

## Data source and citation

Xu Z, Lu Z, Ugurbil A, et al. Mapping convergent regulators of melanoma drug resistance by
PerturbFate. *Nature* 2026;654(8117):261-71.
[doi:10.1038/s41586-026-10367-0](https://www.nature.com/articles/s41586-026-10367-0).
Data: GEO [GSE291147](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291147).

## Previous analysis

[Mechanism report](reports/mechanism-report.zh.md) (in Chinese): research question, key results,
discussion and methods of the earlier analysis.

That analysis investigates differential differentiation responses within a shared YAP-associated
context in melanoma by integrating perturbation information, cell-specific networks, and average
treatment effects. It analyzes A375 RNA data from PerturbFate. Its results indicate that NF2,
MED12, and MED13 knockdowns share enhanced YAP-associated outputs under drug-free conditions,
whereas MED12 exhibits a weaker response regarding drug-induced differentiation programs.

- [Summary results](results/summary): average effects, CKM state annotations, factor annotations,
  and network directionality counts.
- [Core code](icpdag): official NMF/CSCN implementation, Meek-based pre-intervention directional
  assistance, gene-space CKM, and ATE estimation.

CSCN uses the fixed official implementation; NMF serves as a preprocessing step, with factor CKM
mapped back to gene space via the H matrix from the same fit. Per-cell graphs remain independent,
and intervention responses inform directionality prior to Meek propagation.
