"""The held-out target split must be deterministic, disjoint, complete and response-free."""

import importlib.util
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from icpdag.target_split import add_fold, eligible_targets, load_split, make_split

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "results" / "summary" / "phase1_target_qc.csv"
SPLIT = ROOT / "config" / "target_split.json"
spec = importlib.util.spec_from_file_location("make_target_split", ROOT / "scripts" / "make_target_split.py")
mts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mts)


def qc():
    return pd.read_csv(QC, usecols=mts.QC_COLUMNS)


class TargetSplit(unittest.TestCase):
    def test_deterministic(self):
        pd.testing.assert_frame_equal(make_split(qc(), 5, 20261011), make_split(qc(), 5, 20261011))
        self.assertFalse(make_split(qc(), 5, 1).fold.equals(make_split(qc(), 5, 20261011).fold))

    def test_disjoint_and_complete(self):
        split = make_split(qc(), 5, 20261011)
        self.assertEqual(sorted(split.target), eligible_targets(qc()))
        self.assertFalse(split.target.duplicated().any())
        folds = [set(split.target[split.fold.eq(k)]) for k in range(5)]
        for i in range(5):
            for j in range(i + 1, 5):
                self.assertFalse(folds[i] & folds[j])
        self.assertEqual(set().union(*folds), set(split.target))

    def test_balanced_folds_and_tertiles(self):
        split = make_split(qc(), 5, 20261011)
        sizes = split.fold.value_counts()
        self.assertEqual(len(sizes), 5)
        self.assertLessEqual(sizes.max() - sizes.min(), 1)
        counts = pd.crosstab(split.kd_tertile, split.fold)
        self.assertLessEqual((counts.max(axis=1) - counts.min(axis=1)).max(), 1)
        means = split.groupby("kd_tertile").kd_strength.agg(["min", "max"])
        self.assertTrue((means["max"].to_numpy()[:-1] <= means["min"].to_numpy()[1:]).all())

    def test_committed_json_regenerates(self):
        with open(SPLIT, encoding="utf-8") as f:
            committed = json.load(f)
        regenerated = mts.build_spec(QC, committed["n_folds"], committed["seed"])
        self.assertEqual(committed["targets"], regenerated["targets"])
        self.assertEqual(committed["source_sha256"], regenerated["source_sha256"])
        self.assertEqual(sorted(committed["targets"]), eligible_targets(qc()))

    def test_rejects_response_columns(self):
        frame = qc().assign(Fufa2019_score=0.0)
        with self.assertRaisesRegex(ValueError, "Response-like"):
            make_split(frame, 5, 1)

    def test_add_fold_marks_ntc_missing(self):
        split = load_split(SPLIT)
        frame = pd.DataFrame({"target": ["NTC", split.target.iloc[0]], "x": [1, 2]})
        out = add_fold(frame, split)
        self.assertTrue(pd.isna(out.fold.iloc[0]))
        self.assertEqual(out.fold.iloc[1], split.fold.iloc[0])


if __name__ == "__main__":
    unittest.main()
