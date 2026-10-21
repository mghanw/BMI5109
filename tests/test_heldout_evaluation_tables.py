"""Held-out evaluation tables, summary and figures on synthetic inputs (no real responses)."""

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from icpdag import heldout_evaluation as ev
from icpdag.sem_prediction import MODELS

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import evaluate_heldout as eh  # noqa: E402
import summarize_evaluation as se  # noqa: E402

TARGETS = ["CAND1", "DUSP6", "NF2", "PTEN", "RB1", "TP53"]  # names present in the Aim 1 dispersion table


def synthetic(seed=0):
    rng = np.random.default_rng(seed)
    run = json.loads((ROOT / "config" / "evaluation_run.json").read_text(encoding="utf-8"))
    run["bootstrap_draws"] = 200
    endpoints = json.loads((ROOT / "config" / "evaluation_endpoints.json").read_text(encoding="utf-8"))
    split = pd.DataFrame(dict(target=TARGETS, fold=[0, 0, 1, 1, 2, 2], kd_tertile=[0, 1, 2, 0, 1, 2]))
    rows = []
    for m in MODELS:
        for i, t in enumerate(TARGETS):
            for condition in ("DMSO", "Vem"):
                constant = m == "zero_downstream_baseline" or (m == "shared_cscn_sem" and i < 2)
                for c in (0, 1):
                    z = np.nan if constant else rng.normal(0.2, 0.3)
                    rows.append(dict(model=m, target=t, condition=condition, direction=c, r=np.tanh(z), z=z,
                                     rmse=rng.uniform(0.1, 0.3), constant_prediction=constant, constant_observed=False))
    table = ev.combine_metrics(pd.DataFrame(rows))
    table["noise_ceiling_r"] = rng.uniform(-0.2, 0.9, len(table))
    table["fold"] = table.target.map(split.set_index("target").fold)
    stratum = {t: "connected" if t in {"NF2", "PTEN", "TP53"} else "weak" for t in TARGETS}
    table["stratum"] = table.target.map(stratum)

    meta_rows, preds, ntc_ids = [], {}, {}
    n_ref = 40
    for condition in ("DMSO", "Vem"):
        ids = [f"{condition}_ntc{k}" for k in range(n_ref)]
        meta_rows += [dict(cell=i, condition=condition, target="NTC", calib_fold=k % 2) for k, i in enumerate(ids)]
        for t in TARGETS:
            meta_rows += [dict(cell=f"{condition}_{t}{k}", condition=condition, target=t, calib_fold=k % 2)
                          for k in range(12)]
        fufa = rng.normal(size=(len(MODELS), n_ref, len(TARGETS), 2))
        for name in ("shared_cscn_sem", "zero_downstream_baseline"):
            fufa[MODELS.index(name)] = rng.normal(size=(1, len(TARGETS), 2))
        fufa[0, :3, 0, 0] = np.nan
        preds[condition] = dict(targets=TARGETS, fufa_scaled=fufa)
        ntc_ids[condition] = pd.Index(ids)
    meta = pd.DataFrame(meta_rows).set_index("cell")
    scores = rng.normal(size=len(meta))
    return run, endpoints, split, table, stratum, preds, meta, scores, ntc_ids


class Tables(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        run, endpoints, split, table, stratum, preds, meta, scores, ntc_ids = synthetic()
        cls.table = table
        cls.primary = eh.primary_tables(table, TARGETS, split, stratum, run, endpoints)
        cls.secondary = eh.secondary_table(table, TARGETS, run)
        cls.fufa = eh.fufa_table(preds, meta, scores, TARGETS, split, run)
        cls.dist = eh.distribution_table(preds, meta, scores, ntc_ids, TARGETS, split, run)
        cls.run_config = run

    def test_primary_sections(self):
        p = self.primary
        tests = p[p.section.eq("test")].set_index("analysis")
        self.assertEqual(list(tests.index), ["primary_decision", "complete_pairs", "per_condition_DMSO",
                                             "per_condition_Vem", "rmse"])
        self.assertEqual(int(tests.loc["primary_decision", "n_targets"]), 6)
        self.assertEqual(int(tests.loc["primary_decision", "n_nan_dropped"]), 0)
        self.assertEqual(int(tests.loc["complete_pairs", "n_nan_dropped"]), 2)
        self.assertIn(tests.loc["primary_decision", "method"], {"exact", "approx"})
        for section in ("median", "fold_median", "stratum_median", "coverage", "nan_flag"):
            self.assertTrue(p.section.eq(section).any(), section)
        cov = p[p.section.eq("coverage") & p.model.eq("shared_cscn_sem") & p.condition.eq("both")].iloc[0]
        self.assertEqual(int(cov.n_nonconstant), 4)
        self.assertTrue(bool(cov.nan_flag))
        med = p[p.section.eq("median") & p.model.eq("shared_cscn_sem") & p.condition.eq("pooled")]
        self.assertEqual(int(med[med.version.eq("nan_excluded")].n_targets.iloc[0]), 4)
        self.assertEqual(int(med[med.version.eq("r_zero")].n_targets.iloc[0]), 6)
        self.assertFalse(p[p.section.eq("median") & p.model.eq("zero_downstream_baseline")].metric.eq("r").any())

    def test_secondary_holm_families(self):
        s = self.secondary
        self.assertEqual(s.groupby("holm_family").size().to_dict(),
                         {"secondary_rmse": 4, "secondary_z_amended": 3, "secondary_z_complete_pairs": 3})
        for _, g in s.groupby("holm_family"):
            np.testing.assert_allclose(g.p_holm, ev.holm(g.p_value.to_numpy()))

    def test_fufa_and_distribution(self):
        f = self.fufa
        self.assertEqual(f[f.section.eq("test")].contrast.tolist(),
                         ["KD_minus_NTC_DMSO", "KD_minus_NTC_Vem", "interaction_Vem_minus_DMSO"])
        d = self.dist
        cal = d[d.section.eq("calibration") & d.condition.eq("both")]
        self.assertTrue(cal.n_rows.eq(12).all())
        shared = d[d.section.eq("target") & d.model.eq("shared_cscn_sem") & d.statistic.eq("log_var_ratio")]
        np.testing.assert_allclose(shared.predicted, 0.0, atol=1e-12)
        drop = d[d.section.eq("dropped_cells") & d.model.eq(MODELS[0]) & d.condition.eq("DMSO")]
        self.assertEqual(int(drop[drop.direction.eq(0)].n_dropped.iloc[0]), 1)   # 3 NaN cells, 1 in eval fold 1
        self.assertIn("aim1_all_cell_ci_high", d.columns)
        for frame in (self.primary, self.secondary, f, d):
            eh.assert_no_cell_ids(frame)

    def test_summary_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            names = dict(metrics=self.table, primary=self.primary, secondary=self.secondary, fufa=self.fufa,
                         distribution=self.dist)
            tables = {}
            for k, frame in names.items():
                frame.to_csv(tmp / f"{k}.csv", index=False)
                tables[k] = pd.read_csv(tmp / f"{k}.csv")
            tables["metrics"]["direction"] = tables["metrics"].direction.astype(str)
            se.write_summary(tables, dict(date="2026-10-21", head="0" * 40), tmp / "summary.md",
                             {"DMSO": 10.0, "Vem": 10.0})
            text = (tmp / "summary.md").read_text(encoding="utf-8")
            se.style()
            se.paired_z_figure(tables["metrics"], tmp / "a.png")
            se.ceiling_figure(tables["metrics"], tmp / "b.png")
            self.assertTrue((tmp / "a.png").stat().st_size > 0 and (tmp / "b.png").stat().st_size > 0)
        for phrase in ("Amendment (2026-10-21", "Sign convention", "mainly structural", "Spearman partial-correlation",
                       "only calib_fold c cells", "Coverage and accuracy", "NaN flag", "Per-condition", "RMSE",
                       "Per-fold medians", "connectivity stratum", "Holm", "Tertiary"):
            self.assertIn(phrase, text)


if __name__ == "__main__":
    unittest.main()
