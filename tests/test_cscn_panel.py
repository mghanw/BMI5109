"""Network panel selection: NTC-only statistics, ranking, prefix structure and forced genes."""

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd
import scipy.sparse as sp

from icpdag.cscn_panel import (EXCLUDE_PATTERN, build_panel, excluded_symbol, log_cp10k, ntc_gene_stats,
                               panel_of_size)

ROOT = Path(__file__).resolve().parents[1]
PARAMS = dict(conditions=["DMSO", "Vem"], fill_min_detection=0.10, forced_low_detection=0.05, n_bins=4,
              exclude_pattern=EXCLUDE_PATTERN)


def load_aim1():
    spec = importlib.util.spec_from_file_location("aim1_script", ROOT / "scripts" / "aim1.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic(seed=7, n_genes=60):
    rng = np.random.default_rng(seed)
    symbols = ["MT-CO1", "RPL13A", "RPLP0", "RPSA", "TGT1", "TGT2", "FUFA1"] + [f"G{i}" for i in range(n_genes - 7)]
    var = pd.DataFrame({"symbol": symbols, "gene_type": "protein_coding", "is_par_y": False},
                       index=pd.Index([f"ENSG{i:011d}.1" for i in range(n_genes)], name="gene_id"))
    target = ["NTC"] * 160 + ["A"] * 80
    meta = pd.DataFrame({"condition": (["DMSO", "Vem"] * 120), "target": target})
    rates = rng.gamma(0.6, 1.5, n_genes)
    rates[5] = 0.0  # TGT2 never detected: forced, kept and flagged
    shape = rng.uniform(0.3, 3, n_genes)
    lam = rng.gamma(shape[None, :], rates[None, :] / shape[None, :], (len(meta), n_genes))
    counts = sp.csr_matrix(rng.poisson(lam))
    meta["total_counts"] = np.asarray(counts.sum(axis=1)).ravel() + 1000
    forced = pd.DataFrame({"gene_id": var.index[[4, 5, 6]], "symbol": ["TGT1", "TGT2", "FUFA1"],
                           "role": ["target", "target", "Fufa"], "target": ["TGT1", "TGT2", ""],
                           "in_fufa": [False, False, True]})
    return counts, meta, var, forced


class PanelSelection(unittest.TestCase):
    def test_log_cp10k_matches_aim1(self):
        rng = np.random.default_rng(1)
        counts = sp.csr_matrix(rng.poisson(0.7, (25, 12)))
        total = np.asarray(counts.sum(axis=1)).ravel() + 50.0
        ours, theirs = log_cp10k(counts, total), load_aim1().log_cp10k(counts, total)
        self.assertEqual(ours.dtype, theirs.dtype)
        np.testing.assert_array_equal(ours.toarray(), theirs.toarray())

    def test_deterministic_and_forced_first(self):
        counts, meta, var, forced = synthetic()
        a, b = build_panel(counts, meta, var, forced, PARAMS), build_panel(counts, meta, var, forced, PARAMS)
        pd.testing.assert_frame_equal(a["ranked"], b["ranked"])
        self.assertGreater(len(a["ranked"]), 10)
        self.assertEqual(list(a["forced"].gene_id), list(forced.gene_id))
        self.assertTrue(a["forced"].set_index("symbol").low_detection["TGT2"])
        self.assertFalse(set(a["ranked"].symbol) & {"MT-CO1", "RPL13A", "RPLP0", "RPSA", "TGT1", "FUFA1"})

    def test_non_ntc_expression_does_not_change_panel(self):
        counts, meta, var, forced = synthetic()
        reference = build_panel(counts, meta, var, forced, PARAMS)
        changed = counts.tolil()
        rows = np.flatnonzero(meta.target.ne("NTC").to_numpy())
        rng = np.random.default_rng(99)
        changed[rows] = rng.poisson(20, (len(rows), counts.shape[1]))
        other = build_panel(changed.tocsr(), meta, var, forced, PARAMS)
        pd.testing.assert_frame_equal(reference["ranked"], other["ranked"])
        pd.testing.assert_frame_equal(reference["forced"], other["forced"])

    def test_statistics_reject_non_ntc_cells(self):
        meta = pd.DataFrame({"condition": ["DMSO"] * 3, "target": ["NTC", "NTC", "A"]})
        with self.assertRaisesRegex(ValueError, "NTC"):
            ntc_gene_stats(sp.csr_matrix(np.ones((3, 2))), meta, ["DMSO"])

    def test_sizes_are_prefixes_and_keep_forced(self):
        counts, meta, var, forced = synthetic()
        result = build_panel(counts, meta, var, forced, PARAMS)
        ranked = list(result["ranked"].gene_id)
        panels = [panel_of_size(forced.gene_id, ranked, p) for p in (3, 6, 9, 12)]
        for small, large in zip(panels, panels[1:]):
            self.assertEqual(large[:len(small)], small)
        for panel in panels:
            self.assertEqual(panel[:3], list(forced.gene_id))
            self.assertEqual(panel[3:], ranked[:len(panel) - 3])
        with self.assertRaisesRegex(ValueError, "smaller"):
            panel_of_size(forced.gene_id, ranked, 2)

    def test_exclusion_pattern(self):
        excluded = ["MT-CO1", "MT-ND4", "RPL13A", "RPS4Y1", "RPS4X", "RPL36AL", "RPL26L1", "RPS27A",
                    "RPLP0", "RPLP1", "RPLP2", "RPSA", "RPL13AP5", "RPSAP58", "RPLP0P2"]
        kept = ["RPS6KA1", "RPS6KB1", "RPS6KA3", "RPRD1A", "MTOR", "TGT1"]
        np.testing.assert_array_equal(excluded_symbol(excluded), True)
        np.testing.assert_array_equal(excluded_symbol(kept), False)


if __name__ == "__main__":
    unittest.main()
