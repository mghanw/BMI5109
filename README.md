# Intervention-informed CSCN and CKM

结合干预信息、细胞特异网络与平均处理效应，研究黑色素瘤中共同 YAP 相关背景下的分化响应差异。

本项目分析 PerturbFate 的 A375 RNA 数据。当前结果显示，NF2、MED12 和 MED13 敲低共享无药条件下增强的 YAP 相关输出，而 MED12 的药物诱导分化程序响应更弱。

- [机制报告](reports/mechanism-report.zh.md)：研究问题、主要结果、讨论及方法。
- [汇总结果](results/summary)：平均效应、CKM 状态注释、因子注释和网络方向计数。
- [核心代码](icpdag)：官方 NMF/CSCN、Meek 前干预辅助定向、基因空间 CKM 和 ATE 估计。

CSCN 使用固定的官方实现；NMF 作为预处理，因子 CKM 通过同一次拟合的 H 回到基因空间。每细胞图保持独立，干预响应在 Meek 传播前参与定向。

研究数据来源：[PerturbFate](https://www.nature.com/articles/s41586-026-10367-0)。CSCN 来源与许可证见 [第三方代码说明](THIRD_PARTY_NOTICES.md)。
