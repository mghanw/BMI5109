"""Fold-specific orientation: training-only evidence and local IDA enumeration on toy graphs."""

import unittest

import networkx as nx
import numpy as np
import pandas as pd

from icpdag.cscn_networks import reference_model
from icpdag.sem_orientation import (LocalIDAError, local_ida_choices, orient, seed_records, training_evidence,
                                    undirected_pdag)

CFG = dict(window=0.35, ci_alpha=0.05, max_conditioning=2)


def toy_model(nodes):
    return reference_model(np.random.default_rng(0).gamma(2, 1, (20, len(nodes))), nodes, CFG)


def pdag(nodes, undirected=(), directed=()):
    graph = undirected_pdag(nodes, undirected)
    graph.add_edges_from(directed)
    return graph


def screen(seed=3):
    """Two conditions; targets A (knocks down node A, raises C) and H (held out)."""
    rng = np.random.default_rng(seed)
    groups = [("NTC", 40), ("A", 30), ("H", 30)]
    meta = pd.DataFrame([dict(condition=c, target=t) for c in ("DMSO", "Vem") for t, n in groups for _ in range(n)])
    values = pd.DataFrame(rng.normal(5, 1, (len(meta), 4)), columns=["A", "C", "H", "D"])
    hit = meta.target.eq("A")
    values.loc[hit, "A"] -= 3
    values.loc[hit, "C"] += 2
    return values, meta


class Evidence(unittest.TestCase):
    def test_training_rows_only(self):
        values, meta = screen()
        with self.assertRaisesRegex(ValueError, "non-training"):
            training_evidence(values, meta, ["A"], {"A": "A"}, "DMSO")

    def test_held_out_expression_does_not_change_evidence(self):
        values, meta = screen()
        rows = meta.target.isin(["NTC", "A"])
        first = training_evidence(values[rows], meta[rows], ["A"], {"A": "A"}, "DMSO")
        poisoned = values.copy()
        poisoned.loc[meta.target.eq("H")] = 1e6
        second = training_evidence(poisoned[rows], meta[rows], ["A"], {"A": "A"}, "DMSO")
        pd.testing.assert_frame_equal(first, second)
        seeds = {(s["source"], s["destination"]) for s in seed_records(first)}
        self.assertIn(("A", "C"), seeds)
        self.assertNotIn(("A", "A"), seeds)

    def test_seeds_orient_an_undirected_edge(self):
        nodes = ["A", "C", "D"]
        graph = pdag(nodes, undirected=[("A", "C"), ("C", "D")])
        result = orient(graph, [dict(source="A", destination="C", estimate=1.0, standardized_estimate=1.0,
                                     q_value=0.001)], toy_model(nodes))
        self.assertTrue(result["valid"])
        self.assertTrue(result["dag"].has_edge("A", "C"))
        self.assertTrue(result["pdag"].has_edge("C", "D") and not result["pdag"].has_edge("D", "C"))  # Meek R1
        self.assertEqual(result["records"][0]["status"], "accepted_seed")


class LocalIDA(unittest.TestCase):
    def choices(self, nodes, **edges):
        return local_ida_choices(pdag(nodes, **edges), "G", toy_model(nodes))

    def test_isolated_node(self):
        out = self.choices(["G", "A"])
        self.assertEqual((out["n_undirected"], len(out["dags"])), (0, 1))

    def test_chain(self):
        out = self.choices(["A", "G", "B"], undirected=[("A", "G"), ("G", "B")])
        self.assertEqual((out["n_undirected"], out["n_candidates"], len(out["dags"])), (2, 3, 3))
        parents = sorted(tuple(sorted(d.predecessors("G"))) for d in out["dags"])
        self.assertEqual(parents, [(), ("A",), ("B",)])

    def test_star(self):
        out = self.choices(["G", "A", "B", "C"], undirected=[("G", "A"), ("G", "B"), ("G", "C")])
        self.assertEqual((out["n_candidates"], len(out["dags"])), (4, 4))

    def test_undirected_triangle(self):
        out = self.choices(["G", "A", "B"], undirected=[("G", "A"), ("G", "B"), ("A", "B")])
        self.assertEqual((out["n_candidates"], len(out["dags"])), (4, 4))
        self.assertTrue(all(nx.is_directed_acyclic_graph(d) for d in out["dags"]))

    def test_directed_parent_forbids_new_collider(self):
        # P -> G with P not adjacent to A: A -> G would create the new collider A -> G <- P.
        out = self.choices(["P", "G", "A"], undirected=[("G", "A")], directed=[("P", "G")])
        self.assertEqual(out["n_candidates"], 1)
        self.assertTrue(all(d.has_edge("G", "A") for d in out["dags"]))
        # With P adjacent to A, both choices are allowed.
        out = self.choices(["P", "G", "A"], undirected=[("G", "A")], directed=[("P", "G"), ("P", "A")])
        self.assertEqual(out["n_candidates"], 2)

    def test_fallback_completion_when_colliders_cannot_be_preserved(self):
        # An undirected 4-cycle is not chordal: every completion adds a new collider.
        nodes = ["A", "B", "C", "G"]
        cycle = pdag(nodes, undirected=[("A", "B"), ("B", "C"), ("C", "G"), ("G", "A")])
        model = toy_model(nodes)
        result = orient(cycle, [], model)
        self.assertTrue(result["valid"] and result["fallback"])
        self.assertTrue(nx.is_directed_acyclic_graph(result["dag"]))
        self.assertEqual({frozenset(e) for e in result["dag"].edges}, {frozenset(e) for e in cycle.edges})
        out = local_ida_choices(result["pdag"], "G", model)
        self.assertGreater(len(out["dags"]), 0)
        # A chordal graph needs no fallback.
        chain = orient(pdag(["A", "G", "B"], undirected=[("A", "G"), ("G", "B")]), [], toy_model(["A", "G", "B"]))
        self.assertTrue(chain["valid"] and not chain["fallback"])

    def test_too_many_undirected_edges_raise(self):
        leaves = [f"L{i}" for i in range(11)]
        with self.assertRaises(LocalIDAError):
            self.choices(["G", *leaves], undirected=[("G", leaf) for leaf in leaves])


if __name__ == "__main__":
    unittest.main()
