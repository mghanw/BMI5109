"""Committed Aim 1 tables carry the frozen fold of every eligible target and no cell identifiers."""

import re
import unittest
from pathlib import Path

import pandas as pd

from icpdag.target_split import load_split

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "results" / "summary"
PER_TARGET = ["phase2_mean_contrasts.csv", "phase2_group_summaries.csv", "phase2_dispersion_contrasts.csv",
              "phase2_dispersion_depth_sensitivity.csv", "phase2_cluster_bootstrap_fufa.csv",
              "phase2_highlight_targets.csv", "phase2_report_reproduction.csv",
              "phase2_target_in_program_sensitivity.csv", "phase2_expression_summary.csv",
              "phase2_expression_top_genes.csv"]


class Aim1Outputs(unittest.TestCase):
    def test_fold_columns_match_split(self):
        split = load_split(ROOT / "config" / "target_split.json").set_index("target")
        found = [name for name in PER_TARGET if (SUMMARY / name).exists()]
        if not found:
            self.skipTest("Aim 1 tables not generated")
        for name in found:
            frame = pd.read_csv(SUMMARY / name)
            self.assertIn("fold", frame.columns, name)
            eligible = frame[frame.target.isin(split.index)]
            self.assertTrue(len(eligible), name)
            expected = eligible.target.map(split.fold)
            self.assertTrue((eligible.fold.astype(int) == expected).all(), name)
            self.assertTrue(frame.loc[~frame.target.isin(split.index), "fold"].isna().all(), name)

    def test_no_barcodes_in_phase2_tables(self):
        pattern = re.compile(r"[ACGT]{8,}")
        for path in SUMMARY.glob("phase2_*.csv"):
            text = path.read_text(encoding="utf-8")
            self.assertIsNone(pattern.search(text), path.name)


if __name__ == "__main__":
    unittest.main()
