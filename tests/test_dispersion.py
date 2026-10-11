"""Dispersion and distribution statistics on synthetic data with known answers."""

import unittest

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, wasserstein_distance

from icpdag import dispersion as ds


class PointStatistics(unittest.TestCase):
    def test_variance_ratio_and_iqr(self):
        rng = np.random.default_rng(0)
        x, y = rng.normal(0, 2, 200_000), rng.normal(0, 1, 200_000)
        self.assertAlmostEqual(ds.log_var_ratio(x, y), np.log(4), delta=0.02)
        u = rng.uniform(0, 1, 200_000)
        self.assertAlmostEqual(ds.iqr(u), 0.5, delta=0.01)
        self.assertAlmostEqual(ds.iqr_diff(2 * u, u), 0.5, delta=0.01)

    def test_shift_gives_w1_and_zero_centered_w1(self):
        y = np.random.default_rng(1).normal(size=300)
        self.assertAlmostEqual(ds.wasserstein1(y + 0.7, y), 0.7, places=10)
        self.assertAlmostEqual(ds.centered_wasserstein1(y + 0.7, y), 0.0, places=10)

    def test_disjoint_supports_ks_one(self):
        self.assertEqual(ds.ks_stat(np.arange(10) + 100.0, np.arange(7.0)), 1.0)

    def test_exact_against_scipy_unequal_sizes_and_ties(self):
        rng = np.random.default_rng(2)
        for n, m in [(37, 211), (150, 150), (5, 1989)]:
            x = rng.normal(0.3, 1.4, n)
            y = np.round(rng.normal(size=m), 1)  # ties
            x[:3] = y[:3]
            self.assertAlmostEqual(ds.wasserstein1(x, y), wasserstein_distance(x, y), places=10)
            self.assertAlmostEqual(ds.ks_stat(x, y), ks_2samp(x, y).statistic, places=12)

    def test_vectorised_rows_match_scalar(self):
        rng = np.random.default_rng(3)
        xs, ys = np.sort(rng.normal(size=(4, 20)), axis=1), np.sort(rng.normal(size=(4, 33)), axis=1)
        for b in range(4):
            self.assertAlmostEqual(ds.wasserstein_sorted(xs, ys)[b], wasserstein_distance(xs[b], ys[b]))
            self.assertAlmostEqual(ds.ks_sorted(xs, ys)[b], ks_2samp(xs[b], ys[b]).statistic)


class Resampling(unittest.TestCase):
    def test_bootstrap_ci_covers_truth(self):
        rng = np.random.default_rng(4)
        x, y = rng.normal(0, 1.5, 3000), rng.normal(0, 1, 3000)
        a, b = ds.GroupDraws(x, 1000, 7, "c", "A"), ds.GroupDraws(y, 1000, 7, "c", "NTC")
        low, high = np.quantile(np.log(a.var / b.var), [0.025, 0.975])
        self.assertLess(low, np.log(2.25))
        self.assertGreater(high, np.log(2.25))
        self.assertLess(ds.bootstrap_p(np.log(a.var / b.var)), 0.01)

    def test_draws_are_keyed_by_stratum(self):
        np.testing.assert_array_equal(ds.bootstrap_indices(50, 3, 1, "DMSO", "NTC"),
                                      ds.bootstrap_indices(50, 3, 1, "DMSO", "NTC"))
        self.assertFalse(np.array_equal(ds.bootstrap_indices(50, 3, 1, "DMSO", "NTC"),
                                        ds.bootstrap_indices(50, 3, 1, "Vem", "NTC")))

    def test_permutation_p(self):
        rng = np.random.default_rng(5)
        y = rng.normal(size=400)
        self.assertLess(ds.permutation_p(rng.normal(0, 2.5, 120), y, 500, rng, center=True), 0.01)
        self.assertGreater(ds.permutation_p(rng.normal(size=120), y, 500, np.random.default_rng(6)), 0.05)
        # A pure location shift is not a shape difference.
        self.assertGreater(ds.permutation_p(y[:120] + 3, y[120:], 500, np.random.default_rng(7), center=True), 0.05)

    def test_bh(self):
        q = ds.bh([0.01, 0.04, np.nan, 0.03])
        np.testing.assert_allclose(q[[0, 1, 3]], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(q[2]))


class Contrasts(unittest.TestCase):
    def data(self):
        """Planted factorial on the SD: NTC 1 / 1.5 (DMSO / Vem), A 2 / 2 (interaction on log var)."""
        rng = np.random.default_rng(8)
        frames = []
        for condition, sd_ntc, sd_a in [("DMSO", 1.0, 2.0), ("Vem", 1.5, 2.0)]:
            for target, sd, n in [("NTC", sd_ntc, 1500), ("A", sd_a, 600)]:
                frames.append(pd.DataFrame({"condition": condition, "target": target,
                                            "score": rng.normal(0, sd, n)}))
        frame = pd.concat(frames, ignore_index=True)
        return frame[["score"]], frame[["condition", "target"]]

    def test_drug_and_interaction_arithmetic(self):
        scores, meta = self.data()
        out = ds.dispersion_contrasts(scores, meta, ["A"], 200, 99, 11).set_index(
            ["statistic", "kind", "target", "condition"])
        lvr = lambda c, t: ds.log_var_ratio(scores.score[(meta.condition == c) & (meta.target == t)],
                                            scores.score[(meta.condition == c) & (meta.target == "NTC")])
        gene_d = out.loc[("log_var_ratio", "gene", "A", "DMSO"), "estimate"]
        gene_v = out.loc[("log_var_ratio", "gene", "A", "Vem"), "estimate"]
        self.assertAlmostEqual(gene_d, lvr("DMSO", "A"))
        inter = out.loc[("log_var_ratio", "interaction", "A", "Vem_minus_DMSO"), "estimate"]
        self.assertAlmostEqual(inter, gene_v - gene_d)
        self.assertAlmostEqual(inter, np.log(4 / 2.25) - np.log(4), delta=0.15)
        drug_ntc = out.loc[("log_var_ratio", "drug", "NTC", "Vem_minus_DMSO"), "estimate"]
        self.assertAlmostEqual(drug_ntc, np.log(2.25), delta=0.1)
        self.assertTrue(out.q_value.notna().any())
        self.assertTrue((out.ci_low <= out.ci_high).all())


class DepthAndClusters(unittest.TestCase):
    def test_residualize_depth_recovers_slope(self):
        rng = np.random.default_rng(9)
        depth = rng.normal(8, 0.5, 4000)
        score = 0.8 * (depth - 8) + rng.normal(0, 0.3, 4000)
        fit = np.arange(4000) < 2000
        # A group with wider depth spread has inflated variance before adjustment only.
        depth[3000:] = rng.normal(8, 1.5, 1000)
        score[3000:] = 0.8 * (depth[3000:] - 8) + rng.normal(0, 0.3, 1000)
        adjusted, slope, r2 = ds.residualize_depth(score, depth, fit)
        self.assertAlmostEqual(slope, 0.8, delta=0.03)
        self.assertGreater(r2, 0.5)
        self.assertGreater(ds.log_var_ratio(score[3000:], score[:2000]), 1)
        self.assertLess(abs(ds.log_var_ratio(adjusted[3000:], adjusted[:2000])), 0.15)

    def test_cluster_bootstrap_sufficient_statistics_and_width(self):
        rng = np.random.default_rng(10)
        rows = []
        for condition in ("DMSO", "Vem"):
            for w in range(24):
                effect = rng.normal(0, 0.6)  # planted well effect on the target group only
                for target, n in (("NTC", 40), ("A", 10)):
                    for v in rng.normal((effect + 0.2) if target == "A" else 0, 1, n):
                        rows.append(dict(condition=condition, target=target, well=f"{condition}{w}", value=v))
        frame = pd.DataFrame(rows)
        draws = ds.cluster_bootstrap(frame.value, frame[["condition", "target", "well"]], ["A"], 400, 3)
        # Draw with every well once equals the direct statistics.
        d = draws["DMSO"]
        sub = frame[frame.condition.eq("DMSO")]
        np.testing.assert_equal(d["wells"], 24)
        counts = np.ones(24)
        wells = np.sort(sub.well.unique())
        stats = sub.assign(sq=sub.value ** 2).groupby(["well", "target"]).agg(n=("value", "size"),
                                                                                s1=("value", "sum"),
                                                                                s2=("sq", "sum"))
        n = stats.n.unstack().reindex(index=wells, columns=["NTC", "A"]).to_numpy()
        s1 = stats.s1.unstack().reindex(index=wells, columns=["NTC", "A"]).to_numpy()
        s2 = stats.s2.unstack().reindex(index=wells, columns=["NTC", "A"]).to_numpy()
        nn, t1, t2 = counts @ n, counts @ s1, counts @ s2
        direct = sub.groupby("target").value.agg(["mean", "var"]).loc[["NTC", "A"]]
        np.testing.assert_allclose(t1 / nn, direct["mean"])
        np.testing.assert_allclose((t2 - t1 ** 2 / nn) / (nn - 1), direct["var"])
        cluster = ds.cluster_contrasts(draws, ["A"], "score")
        cl = cluster.set_index(["statistic", "kind", "condition"]).loc[("mean", "gene", "DMSO")]
        cell = ds.GroupDraws(sub.value[sub.target.eq("A")], 400, 3, "DMSO", "A").mean - \
            ds.GroupDraws(sub.value[sub.target.eq("NTC")], 400, 3, "DMSO", "NTC").mean
        cell_width = np.subtract(*np.quantile(cell, [0.975, 0.025]))
        self.assertGreater(cl.cluster_ci_high - cl.cluster_ci_low, cell_width)
        self.assertEqual(int(cluster.cluster_draws_dropped.sum()), 0)


if __name__ == "__main__":
    unittest.main()
