# Reproducing the analysis

All commands run from the repository root, in the order below. Each step reads the outputs of the
previous ones. Intermediate files go to `outputs/` and `data/processed/` (not tracked); aggregate tables
go to `results/summary/` and figures to `reports/figures/`. Every script describes its options with
`--help`.

## Environment

- Python 3.11 in a conda environment, then:

  ```bash
  pip install -e ".[plots,io]" pytest
  ```

  The `plots` extra installs matplotlib; the `io` extra installs anndata and openpyxl.
- R >= 4.4 with the `Matrix` package, for converting the GEO RDS files.
- A machine with 32 GB RAM; count matrices are kept sparse.

## 1. Input data

The two small GEO files are in `data/raw/`. Get the RNA count matrix from the
[data-v1 release](https://github.com/mghanw/BMI5109/releases/tag/data-v1), or download all three files
from GEO with the PowerShell script (Windows); details and checksums are in
[data/README.md](../data/README.md).

```powershell
powershell -ExecutionPolicy Bypass -File scripts\download_data.ps1
```

## 2. Convert the RDS matrices to MatrixMarket

```bash
Rscript scripts/rds_to_mtx.R data/raw/GSE291147_Dual_omics_RNA_gene_count_matrix.RDS data/processed/rna_mtx
```

```bash
Rscript scripts/rds_to_mtx.R data/raw/GSE291147_Dual_omics_sgRNA_count_matrix.RDS data/processed/sgrna_mtx
```

## 3. Build cells x features AnnData files

```bash
python scripts/mtx_to_h5ad.py data/processed/rna_mtx data/raw/GSE291147_Dual_omics_RNA_cell_metadata.csv.gz data/processed/rna.h5ad
```

```bash
python scripts/mtx_to_h5ad.py data/processed/sgrna_mtx data/raw/GSE291147_Dual_omics_RNA_cell_metadata.csv.gz data/processed/sgrna.h5ad
```

## 4. Public reference files

Fetch the files listed in [data/ref/SOURCES.md](../data/ref/SOURCES.md) and check the GENCODE gene IDs
against the matrix (this rewrites `SOURCES.md`; only access dates should change):

```bash
python -I scripts/fetch_refs.py --download --gtf-check data/processed/rna_mtx/rows.txt
```

## 5. Data preparation and QC

Perturbation assignment, knockdown QC and calibration folds (thresholds in `config/prepare_data.json`):

```bash
python scripts/prepare_data.py
```

## 6. Held-out target split

The split is frozen; rerunning it must reproduce the committed file exactly:

```bash
python scripts/make_target_split.py
```

```bash
git diff --exit-code config/target_split.json
```

## 7. Program gene sets

Fufa 2019, Kanai 2024 and Tsoi 2018 programs, written to `config/programs.json`:

```bash
python -I scripts/build_programs.py
```

## 8. Aim 1: mean and distribution contrasts

Configuration `config/aim1.json`; then the figures:

```bash
python scripts/aim1.py
```

```bash
python scripts/aim1_figures.py
```

## 9. Score reliability and power

NTC cells only (configuration `config/power_analysis.json`; uses the Aim 1 scores from step 8); then the
figure:

```bash
python scripts/power_analysis.py
```

```bash
python scripts/power_figures.py
```

## 10. Network panel and settings selection grids

The network panel is built from NTC cells and frozen in `config/cscn_panel.json`, so a rerun only
verifies it:

```bash
python scripts/build_cscn_panel.py
```

The settings selection grids compare panel size x maximum conditioning set
(`results/summary/cscn_panel_timing.csv`) and CSCN window x CI alpha
(`results/summary/cscn_panel_sensitivity.csv`) on graph diagnostics and computational feasibility only:

```bash
python scripts/cscn_panel_timing.py
```

```bash
python scripts/cscn_panel_timing.py --sizes 250 --orders 2 --windows 0.35,0.5,0.65 --alphas 0.05,0.1
```

## 11. Cell-specific CSCN graphs

One graph per NTC cell and condition, then the merge per condition:

```bash
python scripts/infer_cscn_networks.py --condition DMSO
```

```bash
python scripts/infer_cscn_networks.py --condition DMSO --merge
```

Repeat both commands with `--condition Vem`, then write the summaries and figure:

```bash
python scripts/summarize_cscn_networks.py
```

## 12. SEM structure

Shared and cluster graphs, principal components and neighbourhoods (NTC cells only), then the
fold-specific orientation diagnostics:

```bash
python scripts/build_sem_graphs.py
```

```bash
python scripts/sem_structure_report.py
```

## 13. Shrinkage lambda

NTC-only cross-validation of the ridge-to-pooled shrinkage:

```bash
python scripts/sem_lambda_cv.py
```

## 14. Knockdown predictions

Model definitions are frozen in `config/sem_prediction.json`. For each condition and outer fold, run the
prediction and then its merge:

```bash
python scripts/predict_sem_knockdown.py --condition DMSO --fold 0
```

```bash
python scripts/predict_sem_knockdown.py --condition DMSO --fold 0 --merge
```

Repeat for folds 0-4 and both conditions, then write the summaries:

```bash
python scripts/summarize_sem_prediction.py
```

## 15. Held-out evaluation

The pre-registered evaluation (`config/evaluation_endpoints.json`, run rules in
`config/evaluation_run.json`). It runs once from committed code; a repeated run requires a dated
deviation note, which is then reported in `results/summary/evaluation_summary.md`.

```bash
python scripts/evaluate_heldout.py
```

```bash
python scripts/summarize_evaluation.py
```

## 16. Tests

```bash
pytest
```
