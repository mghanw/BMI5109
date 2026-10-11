# Intervention-informed CSCN and CKM

This project investigates differential differentiation responses within a shared YAP-associated context in melanoma by integrating perturbation information, cell-specific networks, and average treatment effects.

The project analyzes A375 RNA data from PerturbFate. Current results indicate that NF2, MED12, and MED13 knockdowns share enhanced YAP-associated outputs under drug-free conditions, whereas MED12 exhibits a weaker response regarding drug-induced differentiation programs.

- [Mechanistic report](reports/mechanism-report.zh.md):Research question, key results, discussion, and methods.
- [Summary results](results/summary): Average effects, CKM state annotations, factor annotations, and network directionality counts.
- [Core code](icpdag):Official NMF/CSCN implementation, Meek-based pre-intervention directional assistance, gene-space CKM, and ATE estimation.

CSCN utilizes the fixed official implementation; NMF serves as a preprocessing step, with factor CKM mapped back to gene space via the H matrix from the same fit. Per-cell graphs remain independent, and intervention responses inform directionality prior to Meek propagation.

Data source: PerturbFate.(https://www.nature.com/articles/s41586-026-10367-0). See the third-party code documentation for CSCN sources and licensing information.(THIRD_PARTY_NOTICES.md)。
