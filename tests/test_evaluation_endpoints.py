"""The pre-registered held-out evaluation endpoints are complete and frozen as decided."""

import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "config" / "evaluation_endpoints.json"


class EvaluationEndpoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with open(PATH, encoding="utf-8") as f:
            cls.spec = json.load(f)

    def test_frozen_and_complete(self):
        self.assertTrue(self.spec["frozen"])
        for key in ("response", "evaluation_cells", "gene_set", "primary", "secondary_comparisons",
                    "secondary_endpoint", "tertiary_endpoints", "leakage"):
            self.assertIn(key, self.spec)
        self.assertEqual(self.spec["conditions"], ["DMSO", "Vem"])
        self.assertEqual(self.spec["evaluation_cells"]["scheme"], "cross_fit_calib_fold")
        self.assertIn("the perturbed gene", self.spec["gene_set"]["exclude"])

    def test_primary_comparison_and_decision_test(self):
        p = self.spec["primary"]
        self.assertEqual(p["metrics"], ["pearson_r_across_genes", "rmse"])
        self.assertEqual(p["comparison"], {"model": "cell_specific_cscn_sem", "reference": "shared_cscn_sem"})
        test = p["decision_test"]
        self.assertEqual(test["test"], "wilcoxon_signed_rank")
        self.assertEqual(test["alternative"], "two-sided")
        self.assertEqual(test["alpha"], 0.05)
        self.assertEqual(test["n_tests"], 1)

    def test_fisher_z_averaging(self):
        p = self.spec["primary"]
        self.assertEqual(p["averaging"]["method"], "fisher_z")
        self.assertEqual(p["averaging"]["transform"], "atanh")
        self.assertEqual(p["averaging"]["back_transform"], "tanh")
        self.assertEqual(p["decision_test"]["statistic"], "paired_z_difference")

    def test_hodges_lehmann_effect_size(self):
        e = self.spec["primary"]["effect_size"]
        self.assertEqual(e["method"], "hodges_lehmann")
        self.assertEqual(e["ci"]["method"], "percentile_bootstrap_over_targets")
        self.assertGreaterEqual(e["ci"]["draws"], 1000)
        self.assertIsInstance(e["ci"]["seed"], int)
        self.assertIn("every Wilcoxon test", e["applies_to"])
        self.assertIn("Hodges-Lehmann", self.spec["secondary_comparisons"]["test"])

    def test_noise_ceiling_is_context_only(self):
        n = self.spec["primary"]["noise_ceiling"]
        self.assertFalse(n["decision_test"])
        self.assertEqual(n["splits"], 100)
        self.assertIsInstance(n["seed"], int)
        self.assertEqual(n["correction"], "spearman_brown")
        self.assertIn("separately", n["halves"])

    def test_nan_rule(self):
        r = self.spec["primary"]["nan_rule"]
        self.assertEqual(r["paired_tests"], "complete_pairs")
        self.assertTrue(r["report_nan_counts"])
        self.assertEqual(r["flag_threshold"], 0.1)

    def test_amendment_before_model_fitting(self):
        amendments = self.spec["amendments"]
        self.assertEqual(len(amendments), 1)
        a = amendments[0]
        self.assertEqual(a["date"], "2026-10-21")
        self.assertTrue(a["before_model_fitting"])
        changes = {c["field"]: c["rule"] for c in a["changes"]}
        self.assertIn("r = 0", changes["primary.nan_rule"])
        self.assertIn("all 56", changes["primary.nan_rule"])
        self.assertIn("supporting analysis", changes["primary.nan_rule"])
        self.assertIn("zero_method 'wilcox'", changes["primary.decision_test"])
        self.assertIn("dropped pairs is reported", changes["primary.decision_test"])
        self.assertIn("40 of 56", a["reason"])

    def test_secondary_comparators(self):
        s = self.spec["secondary_comparisons"]
        models = [c["model"] for c in s["comparators"]]
        self.assertEqual(models, ["zero_downstream_baseline", "cluster_level_sem", "local_association_network_sem",
                                  "unweighted_propagation_ablation"])
        self.assertEqual(s["comparators"][0]["metrics"], ["rmse"])
        self.assertEqual(s["multiplicity"], "holm across the 4 comparators")


if __name__ == "__main__":
    unittest.main()
