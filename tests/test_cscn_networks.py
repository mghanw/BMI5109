"""Cell-specific CSCN wrapper: synthetic structure recovery, NTC guard, edge encoding and sharding."""

from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np
import pandas as pd

from icpdag._vendor.cscn.core import CSCN
from icpdag.cscn_networks import (decode_edges, encode_edges, infer_cell, infer_networks, merge_shards,
                                  reference_model, run_shard)

CFG = dict(window=0.35, ci_alpha=0.05, max_conditioning=2)
NODES = ["X", "Y", "Z"]
SCRATCH = Path(__file__).resolve().parents[1] / "outputs"


def adjacency_share(values, focal, pairs):
    model = reference_model(values, NODES, CFG)
    graphs = [decode_edges(*infer_cell(model, f, NODES)[1], NODES) for f in focal]
    return graphs, {p: np.mean([g.has_edge(*p) or g.has_edge(*p[::-1]) for g in graphs]) for p in pairs}


def ntc_meta(n, condition="DMSO"):
    return pd.DataFrame({"condition": condition, "target": "NTC", "calib_fold": np.arange(n) % 2,
                         "lig1_well": [f"L1_{i % 7:02d}" for i in range(n)]},
                        index=[f"cell{i}" for i in range(n)])


class SyntheticStructure(unittest.TestCase):
    def test_chain_removes_the_indirect_edge(self):
        rng = np.random.default_rng(1)
        n = 600
        x = rng.gamma(2, 1, n)
        y = x + rng.gamma(2, 2, n)
        z = y + rng.gamma(2, 2, n)
        _, share = adjacency_share(np.c_[x, y, z], rng.choice(n, 40, replace=False),
                                   [("X", "Z"), ("X", "Y"), ("Y", "Z")])
        self.assertLess(share[("X", "Z")], 0.5)
        self.assertGreater(share[("X", "Y")], 0.5)
        self.assertGreater(share[("Y", "Z")], 0.5)

    def test_collider_is_oriented(self):
        rng = np.random.default_rng(2)
        n = 600
        x, y = rng.gamma(2, 1, n), rng.gamma(2, 1, n)
        z = x + y + rng.gamma(1, 0.5, n)
        graphs, _ = adjacency_share(np.c_[x, y, z], rng.choice(n, 40, replace=False), [])
        oriented = [g.has_edge("X", "Z") and not g.has_edge("Z", "X") and g.has_edge("Y", "Z")
                    and not g.has_edge("Z", "Y") and not g.has_edge("X", "Y") for g in graphs]
        self.assertGreater(np.mean(oriented), 0.5)


class Contracts(unittest.TestCase):
    def test_non_ntc_input_raises(self):
        values = np.random.default_rng(0).gamma(2, 1, (30, 3))
        meta = ntc_meta(30)
        meta.loc["cell3", "target"] = "MED12"
        with self.assertRaisesRegex(ValueError, "NTC"):
            infer_networks(values, meta, NODES, CFG)
        with self.assertRaisesRegex(ValueError, "NTC"):
            run_shard(values, meta, NODES, CFG, SCRATCH / "never_written", 0, 1)

    def test_ci_error_raises_when_strict_and_is_recorded_otherwise(self):
        model = reference_model(np.random.default_rng(0).gamma(2, 1, (30, 3)), NODES, CFG)

        def failing(self, local):
            print("ICT ERROR: synthetic failure")
            return nx.complete_graph(NODES), {}

        with patch.object(CSCN, "_estimate_pc_skeleton", failing):
            with self.assertRaisesRegex(RuntimeError, "ICT ERROR"):
                infer_cell(model, 0, NODES, strict=True)
            row, edges = infer_cell(model, 0, NODES, strict=False)
        self.assertIn("synthetic failure", row["error"])
        self.assertEqual((row["n_edges"], len(edges[0])), (0, 0))

    def test_edge_round_trip(self):
        nodes = list("ABCDE")
        pdag = nx.DiGraph()
        pdag.add_nodes_from(nodes)
        pdag.add_edges_from([("A", "C"), ("B", "C"), ("C", "D"), ("D", "C"), ("E", "A")])
        src, dst, directed = encode_edges(pdag, nodes)
        self.assertEqual(len(src), 4)
        self.assertEqual(int((~directed).sum()), 1)
        back = decode_edges(src, dst, directed, nodes)
        self.assertEqual(set(back.edges), set(pdag.edges))
        self.assertEqual(list(back.nodes), nodes)


class Sharding(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(5)
        n = 90
        a = rng.gamma(2, 1, n)
        b = a + rng.gamma(2, 1, n)
        c = rng.gamma(2, 1, n)
        d = b + c + rng.gamma(1, 1, n)
        self.values, self.meta, self.nodes = np.c_[a, b, c, d], ntc_meta(n), list("ABCD")
        SCRATCH.mkdir(exist_ok=True)
        self.tmp = Path(tempfile.mkdtemp(dir=SCRATCH))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_all(self, name, n_shards, workers=1):
        out = self.tmp / name
        for shard in range(n_shards):
            run_shard(self.values, self.meta, self.nodes, CFG, out, shard, n_shards, workers=workers,
                      chunk_size=7, log=lambda *_: None)
        return merge_shards(out, self.meta)

    @staticmethod
    def stable(cells):
        return cells.drop(columns=["seconds", "seconds_pc"])

    def test_results_do_not_depend_on_sharding_or_workers(self):
        one_cells, one_edges = self.run_all("one", 1)
        three_cells, three_edges = self.run_all("three", 3, workers=2)
        direct_cells, direct_edges = infer_networks(self.values, self.meta, self.nodes, CFG, chunk_size=11)
        pd.testing.assert_frame_equal(self.stable(one_cells), self.stable(three_cells))
        pd.testing.assert_frame_equal(self.stable(one_cells)[direct_cells.columns.drop(["seconds", "seconds_pc"])],
                                      self.stable(direct_cells), check_dtype=False)
        for key in one_edges:
            np.testing.assert_array_equal(one_edges[key], three_edges[key])
            np.testing.assert_array_equal(one_edges[key], direct_edges[key])
        self.assertGreater(len(one_edges["source"]), 0)
        self.assertEqual(list(one_cells.cell_id), list(self.meta.index))
        self.assertEqual(list(one_cells.calib_fold), list(self.meta.calib_fold))

    def test_resume_recomputes_missing_chunks_identically(self):
        cells, edges = self.run_all("resume", 2)
        shard = self.tmp / "resume" / "shard_001_of_002"
        (shard / "chunk_00001.npz").unlink()
        (shard / "done").unlink()
        again_cells, again_edges = self.run_all("resume", 2)
        pd.testing.assert_frame_equal(self.stable(cells), self.stable(again_cells))
        for key in edges:
            np.testing.assert_array_equal(edges[key], again_edges[key])

    def test_merge_refuses_incomplete_or_mismatched_shards(self):
        out = self.tmp / "partial"
        run_shard(self.values, self.meta, self.nodes, CFG, out, 0, 2, chunk_size=7, log=lambda *_: None)
        with self.assertRaisesRegex(ValueError, "shards"):
            merge_shards(out, self.meta)
        with self.assertRaisesRegex(ValueError, "different configuration"):
            run_shard(self.values, self.meta, self.nodes, dict(CFG, ci_alpha=0.01), out, 0, 2, chunk_size=7,
                      log=lambda *_: None)


class SummaryTables(unittest.TestCase):
    def setUp(self):
        import importlib.util

        path = Path(__file__).resolve().parents[1] / "scripts" / "summarize_cscn_networks.py"
        spec = importlib.util.spec_from_file_location("summarize_cscn_networks", path)
        self.s = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.s)
        # Cell 0: T0 -> G2, G2 - G3 (undirected); cell 1: G3 -> G2; cell 2: CI error (ignored)
        self.edges = dict(cell_index=np.array([0, 0, 1, 2], np.int32), source=np.array([0, 2, 3, 0], np.int16),
                          destination=np.array([2, 3, 2, 1], np.int16), directed=np.array([True, False, True, True]))
        self.cells = pd.DataFrame(dict(cell_index=[0, 1, 2], error=["", "", "ICT ERROR: x"]))
        self.nodes = ["T0", "T1", "G2", "G3"]

    def test_edge_frequency_counts_each_orientation(self):
        f = self.s.edge_frequency("DMSO", self.cells, self.edges, self.nodes).set_index(["gene_a", "gene_b"])
        self.assertEqual(set(f.index), {("T0", "G2"), ("G2", "G3")})
        row = f.loc[("G2", "G3")]
        self.assertEqual((row.n_cells, row.share_cells_adjacent), (2, 1.0))
        self.assertEqual((row.share_cells_a_to_b, row.share_cells_b_to_a, row.share_cells_undirected), (0, 0.5, 0.5))
        self.assertEqual(f.loc[("T0", "G2")].share_cells_a_to_b, 0.5)

    def test_target_connectivity_and_two_step_reach(self):
        t = self.s.target_connectivity("DMSO", self.cells, self.edges, self.nodes, [("T0", 0), ("T1", 1)])
        t = t.set_index("target")
        self.assertEqual(t.loc["T0", "share_cells_with_edge"], 0.5)
        self.assertEqual(t.loc["T0", "median_genes_within_2_steps"], 1.0)  # cells: {G2, G3} and {}
        self.assertTrue(t.loc["T1", "isolated_majority"])
        self.assertFalse(t.loc["T0", "isolated_majority"])


if __name__ == "__main__":
    unittest.main()
