"""Comparator graphs: population PC on synthetic structures, coverage matching, neighbourhoods, association."""

import unittest

import numpy as np
from scipy.stats import norm

from icpdag.sem_graphs import (PopulationPC, association_pairs, clusters, coverage_matched_graph,
                               nearest_neighbours, spearman_matrix)

NODES = ["X", "Y", "Z"]


def population_graph(values, alpha=0.01, max_conditioning=2):
    corr, _ = spearman_matrix(values)
    return PopulationPC(NODES, corr, len(values), norm.isf(alpha / 2), max_conditioning).fit()


class PopulationPCTests(unittest.TestCase):
    def test_chain_drops_the_indirect_edge(self):
        rng = np.random.default_rng(1)
        x = rng.normal(size=3000)
        y = x + rng.normal(size=3000)
        z = y + rng.normal(size=3000)
        skeleton, _ = population_graph(np.c_[x, y, z])
        self.assertEqual({frozenset(e) for e in skeleton.edges}, {frozenset("XY"), frozenset("YZ")})

    def test_collider_is_oriented(self):
        rng = np.random.default_rng(2)
        x, y = rng.normal(size=3000), rng.normal(size=3000)
        z = x + y + rng.normal(scale=0.5, size=3000)
        _, pdag = population_graph(np.c_[x, y, z])
        self.assertTrue(pdag.has_edge("X", "Z") and not pdag.has_edge("Z", "X"))
        self.assertTrue(pdag.has_edge("Y", "Z") and not pdag.has_edge("Z", "Y"))

    def test_coverage_matching_is_the_sparsest_graph_with_enough_targets(self):
        rng = np.random.default_rng(4)
        base = rng.normal(size=(800, 3))
        weights = rng.normal(size=(3, 7)) * np.linspace(1.0, 0.05, 7)  # later genes connect only weakly
        values = np.c_[base, base @ weights + rng.normal(size=(800, 7))]
        nodes = [f"G{i}" for i in range(10)]
        targets = nodes[3:]
        for required in (2, 5):
            g = coverage_matched_graph(values, nodes, targets, required, max_conditioning=2)
            self.assertGreaterEqual(g["targets_with_edge"], required)
            sparser = PopulationPC(nodes, spearman_matrix(values)[0], 800, g["z_crit"] * 1.05, 2).fit()[0]
            self.assertLess(sum(sparser.degree(t) > 0 for t in targets), required)
        again = coverage_matched_graph(values, nodes, targets, 5, 2)
        self.assertEqual(set(again["pdag"].edges), set(g["pdag"].edges))
        self.assertEqual(coverage_matched_graph(values, nodes, targets, 0, 2)["edges"], 0)

    def test_coverage_guard(self):
        rng = np.random.default_rng(6)
        values = rng.normal(size=(300, 6))
        with self.assertRaisesRegex(RuntimeError, "runtime guard"):
            coverage_matched_graph(values, [f"G{i}" for i in range(6)], ["G4", "G5"], 2, 2, max_seconds=-1)

    def test_constant_gene_has_no_edges(self):
        rng = np.random.default_rng(5)
        values = np.c_[rng.normal(size=500), np.zeros(500), rng.normal(size=500)]
        values[:, 2] += values[:, 0]
        skeleton, _ = population_graph(values)
        self.assertEqual(skeleton.degree("Y"), 0)


class Neighbourhoods(unittest.TestCase):
    def test_nearest_neighbours_include_self_and_respect_pool(self):
        scores = np.random.default_rng(0).normal(size=(50, 3))
        nn = nearest_neighbours(scores, 10)
        self.assertEqual(nn.shape, (50, 10))
        self.assertTrue((nn[:, 0] == np.arange(50)).all())
        pool = np.arange(0, 50, 2)
        nn = nearest_neighbours(scores, 5, candidates=pool)
        self.assertTrue(np.isin(nn, pool).all())

    def test_single_cluster_is_everything(self):
        scores = np.random.default_rng(0).normal(size=(40, 3))
        self.assertTrue((clusters(scores, 1, 7) == 0).all())

    def test_association_pairs_density_and_ties(self):
        rng = np.random.default_rng(3)
        a = rng.normal(size=200)
        values = np.c_[a, a + 0.1 * rng.normal(size=200), rng.normal(size=200), np.ones(200)]
        pairs = association_pairs(values, 2)
        self.assertEqual(pairs[0], (0, 1))
        self.assertEqual(len(pairs), 2)
        self.assertTrue(all(3 not in p for p in association_pairs(values, 6)))


if __name__ == "__main__":
    unittest.main()
