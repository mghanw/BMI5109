"""Split-half reliability helpers and the dispersion power simulation on synthetic data."""

import unittest

import numpy as np

from icpdag import dispersion as ds
from icpdag import power as pw


class Reliability(unittest.TestCase):
    def test_split_halves_partition(self):
        genes = [f"G{i}" for i in range(45)]
        a, b = pw.split_halves(genes, np.random.default_rng(0))
        self.assertEqual((len(a), len(b)), (22, 23))
        self.assertFalse(set(a) & set(b))
        self.assertEqual(set(a) | set(b), set(genes))
        a2, _ = pw.split_halves(genes, np.random.default_rng(0))
        self.assertEqual(a, a2)

    def test_spearman_brown_known_values(self):
        np.testing.assert_allclose(pw.spearman_brown([0.0, 0.5, 1.0]), [0.0, 2 / 3, 1.0])

    def test_reliability_recovers_analytic_value(self):
        rng = np.random.default_rng(1)
        cells, genes, noise = 20_000, 40, 3.0
        latent = rng.normal(size=(cells, 1))
        expr = latent + rng.normal(scale=noise, size=(cells, genes))
        a, b = pw.split_halves(range(genes), rng)
        r = np.corrcoef(expr[:, a].mean(axis=1), expr[:, b].mean(axis=1))[0, 1]
        truth = 1 / (1 + noise ** 2 / genes)  # signal share of the full-score variance
        self.assertAlmostEqual(float(pw.spearman_brown(r)), truth, delta=0.02)


class Injection(unittest.TestCase):
    def test_mixture_preserves_mean_and_adds_variance(self):
        rng = np.random.default_rng(2)
        x = rng.normal(size=200_000)
        y = pw.inject_mixture(x, 0.2, 1.0, rng.permutation(len(x)))
        self.assertAlmostEqual(y.mean(), x.mean(), places=10)
        self.assertAlmostEqual(np.log(y.var() / x.var()), pw.mixture_true_lvr(0.2, 1.0), delta=0.01)
        with self.assertRaises(ValueError):
            pw.inject_mixture(x, 0.6, 1.0, rng.permutation(len(x)))

    def test_scale_gives_exact_lvr(self):
        x = np.random.default_rng(3).normal(size=500)
        y = pw.inject_scale(x, 1.3)
        self.assertAlmostEqual(y.mean(), x.mean(), places=12)
        self.assertAlmostEqual(ds.log_var_ratio(y, x), 2 * np.log(1.3), places=12)


class Bootstrap(unittest.TestCase):
    def test_counts_match_group_draws(self):
        values = np.random.default_rng(4).normal(size=137)
        g = ds.GroupDraws(values, 300, 11, "DMSO", "X")
        idx = ds.bootstrap_indices(len(values), 300, 11, "DMSO", "X")
        mean, var = pw.bootstrap_mean_var(values, pw.bootstrap_counts(idx, len(values)))
        np.testing.assert_allclose(mean, g.mean, atol=1e-12)
        np.testing.assert_allclose(var, g.var, rtol=1e-9)

    def test_resolution_check(self):
        with self.assertRaises(ValueError):
            pw.check_resolution(2000, [0.05 / 56])
        pw.check_resolution(5000, [0.05, 0.05 / 56])


class Simulation(unittest.TestCase):
    ALPHAS = {"a05": 0.05}

    def run_grid(self, lvrs, reps):
        ntc = np.random.default_rng(5).normal(size=1100)
        scenarios = pw.make_scenarios([], [], lvrs)
        res = pw.simulate(ntc, 100, scenarios, reps, 500, 7, "DMSO", self.ALPHAS.values())
        return pw.summarize(res, scenarios, self.ALPHAS)

    def test_null_power_near_alpha(self):
        out = self.run_grid([], 200)
        self.assertEqual(out.model.tolist(), ["no_injection"])
        self.assertGreaterEqual(out.power_lvr_a05.iloc[0], 0.01)
        self.assertLessEqual(out.power_lvr_a05.iloc[0], 0.10)

    def test_power_increases_with_effect(self):
        out = self.run_grid([0.3, 0.8], 100)
        power = out.power_lvr_a05.to_numpy()
        self.assertTrue(np.all(np.diff(power) >= 0), power)
        self.assertGreater(power[-1], 0.8)
        self.assertTrue((out.power_call_a05 <= out.power_lvr_a05).all())

    def test_min_detectable(self):
        self.assertAlmostEqual(pw.min_detectable([0, 0.2, 0.4], [0.05, 0.6, 1.0]), 0.3)
        self.assertTrue(np.isnan(pw.min_detectable([0, 0.2], [0.05, 0.5])))

    def test_signal_scale(self):
        self.assertAlmostEqual(float(pw.signal_scale_lvr(np.log(1.5), 0.5)), np.log(2.0))
        self.assertAlmostEqual(float(pw.signal_scale_lvr(0.3, 1.0)), 0.3)
        self.assertTrue(np.isnan(pw.signal_scale_lvr(0.3, -0.01)))
        self.assertTrue(np.isnan(pw.signal_scale_lvr(np.nan, 0.5)))


class Outputs(unittest.TestCase):
    def test_summary_tables(self):
        import re
        from pathlib import Path

        import pandas as pd

        summary = Path(__file__).resolve().parents[1] / "results" / "summary"
        paths = [summary / "score_reliability.csv", summary / "power_simulation.csv"]
        if not all(p.exists() for p in paths):
            self.skipTest("Power analysis tables not generated")
        for path in paths:
            self.assertIsNone(re.search(r"[ACGT]{8,}", path.read_text(encoding="utf-8")), path.name)
        rel, power = (pd.read_csv(p) for p in paths)
        self.assertEqual(rel.groupby(["readout", "condition", "variant"]).size().max(), 1)
        self.assertEqual(set(power.condition), {"DMSO", "Vem"})
        self.assertTrue(power.bootstrap.eq(5000).all())


if __name__ == "__main__":
    unittest.main()
