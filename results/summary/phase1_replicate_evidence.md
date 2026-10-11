# Phase 1: can the two screen replicates be recovered?

**Conclusion:** no replicate label can be recovered *within* the Dual-omics RNA data. Pending a
decision on Tri-omics, `calib_fold` (whole first-ligation wells, 2 folds per condition, seeded)
is the leakage-safe calibration split.

## Evidence

**Paper (Nature, open access via Europe PMC PMC13233327).** Methods state that "Two CRISPR screen
replicates were performed in parallel", with transduction, passaging and processing done
independently. Cells were taken on day 6 for PerturbFate and on days 6 and 14 for bulk screens.
Replicate reproducibility is shown only for the bulk screens (Extended Data Fig. 3f). Neither the
paper nor the figure legends say which single-cell experiment (Dual-omics: 99,354 cells;
Tri-omics: 218,678 cells) came from which replicate.

**GEO (GSE291147 SOFT).**
- The 64 Dual-omics GSMs (GSM8829735–GSM8829798, "PerturbFate pooled screen Dual-omics-01…64") are
  PCR wells. Every sample carries identical characteristics ("treatment: DMSO or Vemurafenib"), and
  no sample has a replicate, batch or date field.
- Tri-omics has 128 + 128 GSMs (pre-existing RNA + ATAC; nascent RNA), with the same characteristics.

**Barcode structure (data + authors' code at commit bafdfff).**
- Cell name = `<PCR well prefix>.<ligation-1 bc (10 nt)><ligation-2 bc (10 nt)>`; reverse
  transcription is done in bulk.
- Each of the 64 PCR prefixes mixes both conditions in a near-constant ratio, so prefixes carry no
  sample information.
- The authors assign condition from the ligation-1 barcode through a `condition_table`, which is
  not public.
- All 96 ligation-1 wells are pure in condition. Whitelist indices 1–48 are Vem and 49–96 are DMSO.
  So a replicate split, if it exists in Dual-omics, must be a split of the 48 wells of each
  condition.

**Data-driven checks within Dual-omics** (`phase1_batch_tests.csv`, `phase1_well_batch.csv`).
- *Guide composition per well:* separates DMSO from Vem (drug fitness effect) but shows no
  two-group structure within a condition.
- *Well pseudobulk expression* (NTC singlets and all singlets, top-2,000 variable genes, 5 PCs):
  - Four candidate 24/24 layouts were tested by permutation pseudo-F: contiguous halves, mod-12
    halves, parity and mod-8 halves.
  - No layout is significant consistently across conditions or populations.
  - The isolated p < 0.05 results (Vem NTC parity p = 0.001; DMSO NTC mod-8 p = 0.011) do not
    replicate in the all-singlet pseudobulk, and are expected by chance among 16 layout tests.
  - Unsupervised two-cluster splits isolate 1–5 outlier wells; the 5-well DMSO group is
    lower-depth (PC1–depth Spearman −0.36), not a half-plate.
- *Depth:* PC1 explains only 5–14% of variance; per-well median UMI and detected genes show no
  block structure.

## Interpretation
Three readings fit the evidence:
1. Dual-omics came from one replicate and Tri-omics from the other.
2. Both replicates were pooled before fixation.
3. Both are present in Dual-omics but the batch effect is below detection.

Only reading 1 would give a usable second replicate, and only through the Tri-omics files. Those
were **not downloaded**:

| File | Size |
|---|---|
| GSE291147_Tri_omics_RNA_cell_metadata.csv.gz | 5.8 MB |
| GSE291147_Tri_omics_sgRNA_count_matrix.RDS | 2.8 MB |
| GSE291147_UnlabeledRNA_gene_count_matrix.RDS (pre-existing RNA) | 663 MB |
| GSE291147_LabeledRNA_gene_count_matrix.RDS (nascent RNA) | 403 MB |
| GSE291147_Tri_omics_ATAC_cell_metadata.csv.gz | 6.5 MB |
| GSE291147_Tri_omics_ATAC_peak_count_matrix.h5ad | 3.8 GB |

Tri-omics RNA is split into nascent and pre-existing fractions, so its steady-state profile is not
directly comparable to Dual-omics RNA.
