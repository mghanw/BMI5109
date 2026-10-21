"""SEM coefficients and knockdown simulation: estimator limits, forward substitution, delta, Fufa weights."""

import unittest

import networkx as nx
import numpy as np

from icpdag.sem_prediction import (PooledCoefficients, coefficient_matrix, delta_values, fufa_weights,
                                   knockdown_response, pooled_ols, ridge_to_pooled, unweighted_coefficients)


def random_dag(n_nodes, p, seed):
    rng = np.random.default_rng(seed)
    nodes = [f"N{i}" for i in range(n_nodes)]
    order = rng.permutation(n_nodes)
    dag = nx.DiGraph()
    dag.add_nodes_from(nodes)
    for a in range(n_nodes):
        for b in range(a + 1, n_nodes):
            if rng.random() < p:
                dag.add_edge(nodes[order[a]], nodes[order[b]])
    return dag, nodes


class Estimator(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.X = rng.normal(size=(300, 3))
        self.y = self.X @ [1.0, -2.0, 0.5] + rng.normal(size=300) + 4
        self.pooled = np.array([0.3, 0.3, 0.3])

    def test_lambda_limits(self):
        ols = pooled_ols(self.X, self.y)
        np.testing.assert_allclose(ridge_to_pooled(self.X, self.y, self.pooled, 0.0), ols)
        np.testing.assert_allclose(ridge_to_pooled(self.X, self.y, self.pooled, 1e-12), ols, atol=1e-8)
        np.testing.assert_allclose(ridge_to_pooled(self.X, self.y, self.pooled, 1e12), self.pooled, atol=1e-8)
        np.testing.assert_array_equal(ridge_to_pooled(self.X, self.y, self.pooled, np.inf), self.pooled)

    def test_shrinkage_is_monotone_towards_pooled(self):
        gaps = [np.linalg.norm(ridge_to_pooled(self.X, self.y, self.pooled, lam) - self.pooled)
                for lam in (0.01, 0.1, 1, 10)]
        self.assertTrue(all(a > b for a, b in zip(gaps, gaps[1:])))

    def test_no_parents(self):
        self.assertEqual(ridge_to_pooled(np.zeros((5, 0)), np.ones(5), np.zeros(0), 1.0).shape, (0,))


class Simulation(unittest.TestCase):
    def test_forward_substitution_matches_dense_inverse(self):
        dag, nodes = random_dag(12, 0.35, 2)
        index = {n: i for i, n in enumerate(nodes)}
        values = np.random.default_rng(3).normal(size=(200, len(nodes)))
        coefficients = PooledCoefficients(values, index)
        B = coefficient_matrix(dag, coefficients, index)
        for source in nodes:
            dense = np.linalg.solve(np.eye(len(nodes)) - B, np.eye(len(nodes))[index[source]])
            np.testing.assert_allclose(knockdown_response(dag, source, coefficients, index), dense, atol=1e-10)

    def test_response_is_linear_in_delta(self):
        dag, nodes = random_dag(8, 0.4, 5)
        index = {n: i for i, n in enumerate(nodes)}
        B = coefficient_matrix(dag, unweighted_coefficients, index)
        unit = knockdown_response(dag, nodes[0], unweighted_coefficients, index)
        e = np.zeros(len(nodes))
        e[0] = -0.7
        np.testing.assert_allclose(unit * -0.7, np.linalg.solve(np.eye(len(nodes)) - B, e), atol=1e-12)

    def test_unweighted_rows_sum_to_one(self):
        dag, nodes = random_dag(10, 0.5, 6)
        index = {n: i for i, n in enumerate(nodes)}
        B = coefficient_matrix(dag, unweighted_coefficients, index)
        has_parents = np.array([dag.in_degree(n) > 0 for n in nodes])
        np.testing.assert_allclose(B[has_parents].sum(axis=1), 1.0)
        self.assertTrue((B[~has_parents] == 0).all())

    def test_delta_uses_calibration_fold_means(self):
        d = delta_values([1.0, 2.0, 5.0], [0, 0, 1], [3.0, 4.0, 4.0], [0, 1, 1])
        self.assertEqual(d, {0: 1.5 - 3.0, 1: 5.0 - 4.0})
        with self.assertRaisesRegex(ValueError, "perturbed gene"):
            delta_values(np.ones((3, 2)), [0, 0, 1], [1.0], [0])

    def test_fufa_weights_reproduce_the_linear_score_shift(self):
        rng = np.random.default_rng(7)
        genes = [f"G{i}" for i in range(12)]
        program, controls, scale, centre = ["G0", "G1", "G2"], ["G5", "G6", "G7", "G11"], 0.8, 0.1
        x = rng.normal(size=12)
        change = rng.normal(size=12)

        def score(v):
            col = {g: v[i] for i, g in enumerate(genes)}
            return (np.mean([col[g] for g in program]) - np.mean([col[g] for g in controls]) - centre) / scale

        panel = genes[:10]  # G11 is outside the panel and changes by 0
        w, counts = fufa_weights(program, controls, scale, panel)
        restricted = change.copy()
        restricted[10:] = 0
        self.assertAlmostEqual(w @ restricted[:10], score(x + restricted) - score(x))
        self.assertEqual(counts["controls_in_panel"], 3)


class Models(unittest.TestCase):
    """Model-level contracts: pooled-group models, zero baseline, CV errors, delta and seed leakage."""

    def setUp(self):
        import pandas as pd
        from icpdag.cscn_networks import reference_model
        from icpdag.sem_orientation import undirected_pdag

        rng = np.random.default_rng(11)
        self.nodes = ["A", "B", "C", "D"]
        self.index = {n: i for i, n in enumerate(self.nodes)}
        a = rng.gamma(2, 1, 300)
        b = a + rng.gamma(2, 1, 300)
        c = b + rng.gamma(2, 1, 300)
        self.values = np.c_[a, b, c, rng.gamma(2, 1, 300)]
        self.model = reference_model(self.values, self.nodes, dict(window=0.35, ci_alpha=0.05, max_conditioning=2))
        self.pdag = undirected_pdag(self.nodes, [("A", "B"), ("B", "C")])
        self.pdag.remove_edge("B", "A")  # A -> B fixed, B - C undirected
        self.pd = pd

    def test_single_cluster_equals_shared(self):
        from icpdag.sem_prediction import group_responses
        shared = group_responses({0: self.pdag}, np.zeros(300), self.values, [], ["A", "B"], self.model, self.index)
        cluster = group_responses({0: self.pdag}, np.zeros(300, dtype=int), self.values, [], ["A", "B"], self.model,
                                  self.index)
        for g in ("A", "B"):
            np.testing.assert_allclose(shared[0][g][0], cluster[0][g][0])
        pooled = PooledCoefficients(self.values, self.index)
        self.assertAlmostEqual(shared[0]["A"][0][self.index["B"]], pooled("B", ("A",))[0])

    def test_zero_baseline(self):
        from icpdag.sem_prediction import zero_response
        np.testing.assert_array_equal(zero_response("C", self.index), [0, 0, 1, 0])

    def test_cv_errors_at_infinite_lambda_use_pooled_slopes(self):
        from icpdag.sem_prediction import node_prediction_errors
        dag = nx.DiGraph([("A", "B"), ("B", "C")])
        dag.add_nodes_from(self.nodes)
        pooled = PooledCoefficients(self.values[:150], self.index)
        local = self.values[150:250]
        x = self.values[299]
        sse, n = node_prediction_errors(dag, x, local, pooled, [np.inf], self.index)
        self.assertEqual(n, 2)
        expected = 0.0
        for node, parent in (("B", "A"), ("C", "B")):
            b = pooled(node, (parent,))
            a = local[:, self.index[node]].mean() - local[:, self.index[parent]].mean() * b[0]
            expected += (x[self.index[node]] - a - b[0] * x[self.index[parent]]) ** 2
        self.assertAlmostEqual(sse[0], expected)

    def test_delta_reads_only_the_perturbed_gene(self):
        from icpdag.sem_prediction import target_delta
        rng = np.random.default_rng(2)
        meta = self.pd.DataFrame({"condition": "DMSO", "target": ["NTC"] * 40 + ["T"] * 20,
                                  "calib_fold": [0, 1] * 30})
        values = self.pd.DataFrame(rng.normal(size=(60, 3)), columns=["G", "X", "Y"])
        first = target_delta(values, meta, "T", "G", "DMSO")
        poisoned = values.copy()
        poisoned.loc[meta.target.eq("T"), ["X", "Y"]] = 1e9
        self.assertEqual(first, target_delta(poisoned, meta, "T", "G", "DMSO"))

    def test_held_out_cells_never_reach_the_seeds(self):
        from icpdag.sem_orientation import fold_seeds
        rng = np.random.default_rng(3)
        meta = self.pd.DataFrame({"condition": "DMSO", "target": ["NTC"] * 40 + ["T1"] * 20 + ["H"] * 20})
        values = self.pd.DataFrame(rng.normal(size=(80, 3)), columns=["T1", "H", "Z"])
        values.loc[meta.target.eq("T1"), "T1"] -= 3
        values.loc[meta.target.eq("T1"), "Z"] += 2
        seeds, evidence = fold_seeds(values, meta, ["T1"], {"T1": "T1"}, "DMSO")
        poisoned = values.copy()
        poisoned.loc[meta.target.eq("H")] = 1e9
        seeds2, evidence2 = fold_seeds(poisoned, meta, ["T1"], {"T1": "T1"}, "DMSO")
        self.pd.testing.assert_frame_equal(evidence, evidence2)
        self.assertEqual(seeds, seeds2)
        self.assertNotIn("H", set(evidence.target))


class SummaryTable(unittest.TestCase):
    def test_exclusion_flag_and_ida_statistics(self):
        import importlib.util
        from pathlib import Path

        import pandas as pd

        path = Path(__file__).resolve().parents[1] / "scripts" / "summarize_sem_prediction.py"
        spec = importlib.util.spec_from_file_location("summarize_sem_prediction", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        pred = pd.DataFrame(dict(model="m", condition="DMSO", fold=0, target=["A", "B"], included_cells=[90, 100],
                                 reference_cells=100, genes_nonzero=[0, 4], constant_prediction=[True, False]))
        stats = pd.DataFrame(dict(model="m", condition="DMSO", target=["A", "B"], units=[100, 100], excluded=[10, 2],
                                  n_undirected=[0, 50], n_candidates=[100, 150], n_failed=[10, 2], n_fallback=[0, 37]))
        row = module.model_summary(pred, stats).iloc[0]
        self.assertEqual((row.n_constant_prediction, row.constant_targets), (1, "A"))
        self.assertEqual((row.n_targets_exclusion_flag, row.exclusion_flag_targets), (1, "A"))
        self.assertAlmostEqual(row.excluded_share, 12 / 200)
        self.assertAlmostEqual(row.mean_undirected_per_unit, 50 / 200)
        self.assertAlmostEqual(row.fallback_share_of_valid_choices, 37 / 238)


if __name__ == "__main__":
    unittest.main()
