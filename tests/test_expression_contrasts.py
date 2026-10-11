"""Sparse pseudobulk contrasts equal dense calculations; filters and exclusions hold."""

import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ttest_ind

from icpdag.expression import detected_genes, expression_contrasts, group_moments, summarize_expression


def data():
    rng = np.random.default_rng(0)
    rows, dense = [], []
    for condition in ("DMSO", "Vem"):
        for target, n in (("NTC", 60), ("A", 25), ("B", 30)):
            block = rng.poisson(0.6, (n, 5)).astype(float)
            block[:, 4] = 0
            block[:3, 4] = 1  # rare gene: detected in < 5% of NTC
            if target == "A":
                block[:, 0] = 0  # perturbed gene
                block[:, 1] += 2
            dense.append(np.log1p(block))
            rows += [dict(condition=condition, target=target)] * n
    X = np.vstack(dense)
    genes = pd.DataFrame({"gene_id": [f"ENSG{i}.1" for i in range(5)], "symbol": list("VWXYZ")})
    return sp.csr_matrix(X), X, genes, pd.DataFrame(rows)


class Expression(unittest.TestCase):
    def test_moments_match_dense(self):
        Xs, X, _, meta = data()
        labels = list(zip(meta.condition, meta.target))
        levels = [("DMSO", "A"), ("Vem", "NTC")]
        n, mean, mean_sq, detection = group_moments(Xs, labels, levels)
        for k, level in enumerate(levels):
            rows = [i for i, lab in enumerate(labels) if lab == level]
            self.assertEqual(n[k], len(rows))
            np.testing.assert_allclose(mean[k], X[rows].mean(axis=0))
            np.testing.assert_allclose(mean_sq[k], (X[rows] ** 2).mean(axis=0))
            np.testing.assert_allclose(detection[k], (X[rows] > 0).mean(axis=0))

    def test_detection_filter(self):
        Xs, _, _, meta = data()
        flags = detected_genes(Xs, meta, 0.06)
        self.assertEqual(flags.any(axis=1).tolist(), [True, True, True, True, False])

    def test_contrasts_welch_and_exclusion(self):
        Xs, X, genes, meta = data()
        out = expression_contrasts(Xs, genes, meta, ["A", "B"], {"A": "ENSG0.7", "B": "ENSG9.1"})
        a = out[(out.condition == "DMSO") & (out.target == "A")]
        self.assertNotIn("ENSG0.1", set(a.gene_id))
        self.assertEqual(len(out[(out.condition == "DMSO") & (out.target == "B")]), 5)
        row = a.set_index("symbol").loc["W"]
        ta = X[(meta.condition == "DMSO").to_numpy() & (meta.target == "A").to_numpy(), 1]
        tn = X[(meta.condition == "DMSO").to_numpy() & (meta.target == "NTC").to_numpy(), 1]
        self.assertAlmostEqual(row["diff"], ta.mean() - tn.mean())
        ref = ttest_ind(ta, tn, equal_var=False)
        self.assertAlmostEqual(row.welch_t, ref.statistic)
        self.assertAlmostEqual(row.p_value, ref.pvalue)
        summary, tops = summarize_expression(out, 0.5)
        self.assertGreaterEqual(int(summary.set_index(["condition", "target"]).loc[("DMSO", "A"),
                                                                                  "n_absdiff_ge_cutoff_q_lt_005"]), 1)
        self.assertEqual(tops.iloc[0].symbol, "W")


if __name__ == "__main__":
    unittest.main()
