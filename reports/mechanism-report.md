# MED12 knockdown attenuates the drug-induced differentiation response of melanoma cells within a shared YAP-associated context

PerturbFate preliminary research report · 2026-10-11

Translated from the Chinese original, [mechanism-report.zh.md](mechanism-report.zh.md).

## Abstract

Different gene perturbations can converge on similar drug-resistance-associated states, yet their effects on drug-induced differentiation may still differ. We integrated average treatment effect (ATE) analysis, cell-specific causal networks (CSCN) and a gene-space Causal Katz Matrix (CKM) to analyse the intervention responses of 7,125 A375 cells in PerturbFate. NF2, MED12 and MED13 knockdowns all enhance a YAP/TAZ-associated program without drug, whereas the drug-induced differentiation response is markedly attenuated after MED12 knockdown. The differentiation program increases by 1.128, 0.847 and 0.842 SD in the NTC, NF2 and MED13 backgrounds, respectively, and by 0.293 SD in the MED12 background. CKM further resolves a shared drug-associated redistribution of states, and CSCN establishes, in a subset of cells, an intervention-compatible orientation between NF2 and the YAP-associated program. These results support a unified interpretation: **a shared YAP-associated output and perturbation-specific differentiation responses can coexist, and MED12 is an important entry point linking these two layers of phenomena.**

## 1. Background and question

The response of melanoma to MAPK inhibition is closely tied to the cell differentiation state. Through multi-omic perturbation analysis, PerturbFate proposed that perturbations of Hippo/YAP repressors and of the Mediator kinase module can converge on a shared dedifferentiation and drug-resistance-associated program; within this shared context, MED12 knockdown additionally shows a stronger impairment of SOX10-associated differentiation function.[1]

This finding raises a specific question: **do perturbations that share a YAP-associated output retain different degrees of drug-induced differentiation capacity?** We chose NF2 as a control for the Hippo/YAP branch and MED13 as a control within the Mediator kinase module, and built a quantitative comparison around the differentiation response of MED12. ATE measures the magnitude of program responses, CKM organises cell states, and CSCN links local network structure with intervention information, forming a chain of evidence from population responses to cellular heterogeneity to candidate regulatory directions.

## 2. Results

### 2.1 Several classes of perturbation share a YAP-associated output without drug

We analysed NTC and seven eligible knockdown targets under DMSO and vemurafenib (Vem). Using an external melanoma YAP/TAZ-up program score,[3] the effects of NF2, MED12 and MED13 knockdown relative to NTC are +1.588, +1.028 and +0.527 SD, respectively. LATS2 and PTPN14 show the same direction, with effects of +0.828 and +0.699 SD. This multi-perturbation pattern is consistent with the shared Hippo–Mediator output proposed in the original paper.

The drug environment further separates these responses. Under Vem, the YAP program effects of NF2 and MED12 relative to NTC are +0.453 [0.308, 0.599] and +0.554 [0.438, 0.671] SD, respectively; MED13 is +0.044 [−0.064, 0.153] SD, close to the control level. The shared no-drug context thus shows different dependence on the drug condition.

### 2.2 The differentiation program increases less after MED12 knockdown than after NF2 and MED13 knockdown

We used an external SOX10-bound pigmentation differentiation program[2] and compared the mean Vem−DMSO response within each intervention background. The program rises by +1.128 SD in NTC, and by +0.847 and +0.842 SD in NF2 and MED13, respectively; in MED12 the increase is +0.293 SD.

| Intervention background | DMSO / Vem cells | No-drug YAP program: KD−NTC | Differentiation program: Vem−DMSO |
| --- | ---: | ---: | ---: |
| NTC | 1,989 / 1,864 | Reference | +1.128 [1.057, 1.199] |
| NF2 KD | 387 / 171 | +1.588 [1.442, 1.736] | +0.847 [0.639, 1.056] |
| MED12 KD | 162 / 203 | +1.028 [0.858, 1.205] | **+0.293 [0.082, 0.503]** |
| MED13 KD | 273 / 170 | +0.527 [0.380, 0.672] | +0.842 [0.630, 1.046] |

Values are in units of each readout's frozen DMSO NTC standard deviation; brackets are 95% cell bootstrap intervals (5,000 resamples).

The gene × drug interaction of MED12 is **−0.835 [−1.058, −0.613] SD**. Comparing response magnitudes directly, MED12−NF2 is **−0.554 [−0.850, −0.257] SD** and MED12−MED13 is **−0.549 [−0.839, −0.257] SD**. These differences show that the attenuated differentiation response of MED12 is distinct both from the Hippo branch and from MED13 in the same Mediator module.

This difference is especially clear at the program level. The gene × drug interaction of SOX10 RNA in MED12 is −0.012 [−0.207, 0.184] SD, whereas the differentiation program shows a large negative interaction. Interpreting the MED12 phenotype therefore requires considering both transcription factor expression and downstream program output; the drug response of the SOX10 transcript alone is not sufficient to summarise the whole differentiation response.

![YAP-associated output, differentiation program and gene × drug interaction](figures/average-effects.png)

**Figure 1 | Shared output context and differences in differentiation response.** a, KD−NTC effects on the YAP/TAZ program; b, KD−NTC effects on the differentiation program; c, gene × drug interaction of the differentiation program. Points and lines show effect estimates and 95% cell bootstrap intervals, respectively.

### 2.3 Gene-space CKM resolves a shared drug-associated redistribution of states

To place the mean responses back in the context of cellular heterogeneity, we built CKMs through the official NMF preprocessing and per-cell CSCN, then multiplied them back by the loading matrix H of the same NMF run to obtain a cell × gene representation.[5] The final matrix contains 2,005 gene coordinates, and 7,049 cells entered the state analysis. Six states were derived from this representation and biologically annotated with external expression programs outside the input features.[2–4]

S3 has high YAP/TAZ-up, undifferentiated and cell cycle readouts. In DMSO, the proportion of S3 is 20.8% in NTC and 46.4%, 36.4% and 34.2% in NF2, MED12 and MED13, respectively. Under Vem, the S3 proportion decreases in every group, and more cells are distributed in S1 and S4; S4 has a high SOX10-bound pigmentation differentiation readout.

| Intervention background | S3 in DMSO | S3 in Vem | S4 in Vem |
| --- | ---: | ---: | ---: |
| NTC | 20.8% | 0.8% | 52.3% |
| NF2 KD | 46.4% | 1.2% | 37.4% |
| MED12 KD | 36.4% | 3.0% | 43.3% |
| MED13 KD | 34.2% | 0.6% | 45.3% |

These results make the shared YAP-associated context and the drug response concrete as changes in the distribution of cell states. The ranking of S4 proportions is not exactly the same as that of the continuous differentiation program response, indicating that coarse-grained state composition and program response magnitude capture information at different levels. Combined with ATE, CKM provides a basis for further locating in which cell states the MED12 differentiation difference arises.

### 2.4 Intervention-assisted CSCN links NF2 with the YAP-associated program

Among the native NMF factors, the top 50 high-loading genes of `factor_002` include 8 external YAP/TAZ-up members, including CCN1, CCN2, TGM2 and ABL2. With the 1,998 input genes as background and BH correction over the annotations of 20 factors against two programs, the q value of this enrichment is 1.33×10⁻⁵.

After collider orientation and before Meek propagation, we introduced compatible orientations supported by intervention responses in the discovery set. Of the 1,920 valid NTC reference graphs in DMSO, 162 retain `NF2 → factor_002`: in 146 of them the direction is already given by the collider structure, and in 16 it is oriented by intervention responses on an existing adjacency. This direction is also observed in Vem, indicating that the connection has local support in both drug environments.

| Condition | Valid NTC reference graphs | NF2 → YAP-associated factor | Opposite direction | Undirected adjacency |
| --- | ---: | ---: | ---: | ---: |
| DMSO | 1,920 | 162 | 65 | 1 |
| Vem | 1,859 | 39 | 26 | 0 |

These counts retain the differences between cell-specific networks. NF2 knockdown is accompanied by an increase in this factor, which is compatible with the Hippo/YAP biological context; local directional information thus connects with the population program response.

## 3. Discussion

The central finding of this study is that a shared YAP-associated output does not correspond to an identical differentiation response. NF2, MED12 and MED13 share program changes without drug, whereas MED12 shows a weaker increase of the differentiation program in the drug environment. Compared with examining drug-resistance-associated markers alone, this result shifts attention to how cells respond to the drug and to their differentiation plasticity.

This interpretation continues and quantifies the mechanistic framework of the original PerturbFate paper. Combining transcriptome, chromatin and functional experiments, the original paper proposed a shared Mediator–YAP output and a MED12-associated impairment of SOX10 differentiation function.[1] Using externally defined programs and direct response comparisons, we further describe how MED12 differs from NF2 and MED13, and use CKM and CSCN to place this difference in the context of cell states and local regulation. The different dependence on MED12 and MED13 in the VEGFC functional experiments of the original paper also supports distinguishing perturbations within the same complex.

In this framework, ATE, CKM and CSCN take the roles of effect quantification, state organisation and directional interpretation, respectively. At this stage, the quantitative evidence for the MED12 differentiation difference is clearest, and the NF2–YAP branch has obtained local structural support; how MED12 affects SOX10-associated output through specific regulatory connections is the mechanistic question that naturally follows. A next step could first test whether the attenuated differentiation response mainly corresponds to changes in state proportions or also appears within similar CKM states, and then combine SOX10 restoration with YAP/TEAD inhibition to distinguish regulatory paths. Mechanisms such as MED12–TGF-β receptor are also worth considering as comparison directions.[6]

This report is based on cells captured on day 6 that passed quality control; bootstrap intervals describe variation within the current cell pool. The three analyses use the same RNA data, and CKM is derived from CSCN, so they provide complementary computational evidence. Direct mediation relationships and growth effects still need to be determined with independent replicates and functional experiments. Within this scope, the results support a working model in which a shared YAP-associated context and perturbation-specific differentiation responses coexist.

## 4. Methods overview

**Data and readouts.** The analysis includes 7,125 A375 cells, comprising NTC and seven targets that passed knockdown screening: DUSP6, LATS2, MED12, MED13, NF1, NF2 and PTPN14. CP10k was computed from the total of the full transcriptome and log1p-transformed. A program score is the mean of the target genes minus the mean of expression-matched background genes, then standardised to DMSO NTC. The differentiation program comes from Fufa's 501mel/AZD6244 SOX10-bound pigmentation differentiation gene set, with 44 of the 45 source genes mapped; the YAP and state programs come from Kanai and Tsoi, respectively.[2–4]

**ATE.** Let μ(g,d) be the mean readout under intervention g and drug d; the drug response is δg=μ(g,Vem)−μ(g,DMSO), and the gene × drug interaction is δg−δNTC. A 5,000-draw cell bootstrap was stratified by drug × target, with program membership, background genes and scale fixed, and all contrasts sharing the resampled group means. Estimation uses marginal between-group contrasts, without adjustment for post-intervention CKM states. The full tables report the estimates, their intervals and their identification scope.

**CSCN and CKM.** The official CSCN is pinned at version `d579c03d77a6435c570957188de356b72e369f00`. Native `use_nmf=True` was run on 1,998 selected genes, retaining 20 factors; the seven measured target genes were added as separate anchor nodes. The per-cell reference set is stratified by drug × target, with a CI significance level of 0.05 and a maximum conditioning set size of 4. Intervention directions come from discovery-set responses (global BH q<0.05) and knockdown quality screening, and are added before Meek propagation while keeping the structure compatible. The one-sided rank window fraction of 0.35 was chosen by a limited structural check on the same focal cells. Every graph is kept separately; no consensus voting graph was constructed.

The factor CKM uses the official Katz solver (α=0.05), followed by `C_gene = C_factor @ H`, keeping the seven anchor gene columns. The 76 cells for which no consistent DAG extension could be obtained are stored with a missing flag, and the remaining 7,049 cells are used for the state analysis. The final gene representation was standardised, reduced to 20 principal components and clustered by six-cluster KMeans into S1–S6, with random seed 20261011. S3 also contains cell cycle variation, and S1/S5 have lower library complexity; these factors are used in interpreting the states.

**Expression programs and network features.** CKM state annotation uses expression readouts of genes that did not enter the NMF or the anchor nodes. In the back-projected matrix the Fufa set is covered by 9 genes, and the direction of the drug change of their network feature mean differs from that of the full expression program; the related results are retained as well. Expression effects and CKM features are each interpreted on their own scale. Responses of MED12 and of the factor containing SOX10/GAS7 were detected, but the current NTC reference graphs do not retain the corresponding anchor adjacency; the direct structure of this differentiation branch is within the scope of later mechanistic testing.

## 5. Result materials

| Content | File |
| --- | --- |
| Full mean effects, drug responses and interactions | [Average effects table](../results/summary/average_effects.csv) |
| CKM state proportions and numbers of valid cells | [State summary](../results/summary/ckm_state_summary.csv) |
| Expression program and quality annotations of states | [Biological annotations](../results/summary/ckm_state_annotations.csv), [quality annotations](../results/summary/ckm_state_quality.csv) |
| Gene CKM program means | [Feature summary](../results/summary/ckm_gene_program_means.csv) |
| Native factor loading annotations | [All 40 annotation results](../results/summary/native_factor_annotations.csv) |
| Direction counts by condition and connection | [Direction summary](../results/summary/cscn_direction_summary.csv) |

## References

1. Xu et al. *Mapping convergent regulators of melanoma drug resistance by PerturbFate*. Nature (2026). [Main text, Fig. 5 and Extended Data Figs. 8–10](https://www.nature.com/articles/s41586-026-10367-0).
2. Fufa et al. (2019). [Original article and Table S2](https://pmc.ncbi.nlm.nih.gov/articles/PMC6688322/).
3. Kanai et al. (2024). [Original article and Table S3](https://pmc.ncbi.nlm.nih.gov/articles/PMC10930532/).
4. Tsoi et al. (2018). [Original article and Table S3](https://pmc.ncbi.nlm.nih.gov/articles/PMC5953834/).
5. CSCN, Bioinformatics (2026), DOI: 10.1093/bioinformatics/btag480. [Supplementary material S1.1.2, Fig. S2](https://academic.oup.com/bioinformatics/article/42/7/btag480/8722299#supplementary-data).
6. Huang et al. (2012). [Study of the MED12 and TGF-β receptor mechanism](https://pmc.ncbi.nlm.nih.gov/articles/PMC3672971/).
