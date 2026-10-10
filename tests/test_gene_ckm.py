"""Contracts at the official CSCN, orientation, and gene-projection boundaries."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np
import pandas as pd

from icpdag._vendor.cscn.core import CSCN
from icpdag.gene_ckm import fit_gene_ckm, project_to_genes, save_result, _infer_group
from icpdag.responses import response_table


class GeneCKMContracts(unittest.TestCase):
    def test_upstream_bytes_match_recorded_commit(self):
        root = Path(__file__).resolve().parents[1] / "icpdag/_vendor/cscn"
        record = json.loads((root / "upstream.json").read_text())
        for filename, metadata in record["files"].items():
            content = (root / filename).read_bytes()
            blob = b"blob " + str(len(content)).encode() + b"\0" + content
            self.assertEqual(hashlib.sha1(blob).hexdigest(), metadata["git_blob_sha"])

    def test_projection_aligns_factor_names_and_preserves_anchors(self):
        nodes = pd.DataFrame({"F2": [3.], "A": [7.], "F1": [2.]}, index=["cell"])
        h = pd.DataFrame([[1., 4.], [2., 5.]], index=["F1", "F2"], columns=["X", "Y"])
        projected = project_to_genes(nodes, h, ["A"])
        self.assertEqual(list(projected), ["X", "Y", "A"])
        np.testing.assert_array_equal(projected.to_numpy(), [[8., 23., 7.]])
        with self.assertRaisesRegex(ValueError, "excluded"):
            project_to_genes(nodes, h.rename(columns={"X": "A"}), ["A"])

    def test_native_nmf_zero_alpha_reconstructs_gene_coordinates(self):
        rng = np.random.default_rng(41)
        x = pd.DataFrame(rng.gamma(2, 1, (64, 9)), columns=[f"G{i}" for i in range(9)])
        meta = pd.DataFrame({"condition": ["DMSO"]*32+["Vem"]*32, "target": ["NTC"]*64})
        result = fit_gene_ckm(x, meta.iloc[::-1], anchors=["G0"], components=3,
                              katz_alpha=0, n_states=3, workers=2)
        valid = result["cells"].valid_completion.to_numpy()
        self.assertGreater(valid.sum(), 0)
        expected = result["W"].to_numpy(dtype=float) @ result["H"].to_numpy(dtype=float)
        np.testing.assert_allclose(result["ckm_genes"].loc[valid, x.columns[1:]], expected[valid])
        np.testing.assert_allclose(result["ckm_genes"].loc[valid, "G0"], x.loc[valid, "G0"], rtol=1e-6)
        self.assertEqual(result["ckm_genes"].shape, x.shape)
        self.assertEqual(len(result["states"]), valid.sum())
        scratch = Path(__file__).resolve().parents[1] / "outputs"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as directory:
            self.assertTrue(Path(directory).resolve().is_relative_to(scratch.resolve()))
            save_result(result, directory)
            saved = pd.read_csv(Path(directory) / "ckm_genes.csv.gz", index_col=0)
            np.testing.assert_allclose(saved, result["ckm_genes"], equal_nan=True)
            with self.assertRaisesRegex(ValueError, "empty"):
                save_result(result, directory)

    def job(self, nodes):
        return dict(nodes=nodes, data=np.random.default_rng(3).gamma(2, 1, (30, len(nodes))),
                    config=dict(window=.35, ci_alpha=.05, max_conditioning=4, katz_alpha=.05),
                    condition="DMSO", target="NTC", indices=np.arange(30), focals=[0], evidence=[])

    def test_seed_precedes_meek_propagation(self):
        nodes = ["A", "F1", "F2"]
        skeleton = nx.Graph([("A", "F1"), ("F1", "F2")])
        job = self.job(nodes)
        job["evidence"] = [dict(source="A", destination="F1", q_value=.001)]
        with patch.object(CSCN, "_estimate_pc_skeleton", return_value=(
                skeleton, {frozenset(("A", "F2")): {"F1"}})):
            result = _infer_group(job)
        p = result["pdag"][0]
        self.assertEqual((p[0, 1], p[1, 0], p[1, 2], p[2, 1]), (1, 0, 1, 0))
        self.assertEqual(result["attempts"][0]["status"], "accepted_seed")

    def test_invalid_extension_keeps_ckm_missing(self):
        nodes = list("ABCD")
        cycle = nx.cycle_graph(nodes)
        pdag = nx.DiGraph(cycle)
        with patch.object(CSCN, "_estimate_pc_skeleton", return_value=(cycle, {})), \
             patch.object(CSCN, "_orient_unshielded_colliders", return_value=pdag), \
             patch.object(CSCN, "_apply_meek_rules", return_value=pdag):
            result = _infer_group(self.job(nodes))
        self.assertFalse(result["rows"][0]["valid_completion"])
        self.assertTrue(np.isnan(result["ckm"]).all())

    def test_validation_responses_do_not_select_discovery_seeds(self):
        rng = np.random.default_rng(3)
        meta = pd.DataFrame({"condition": ["DMSO"]*120,
                             "target": ["NTC"]*30+["A"]*30+["NTC"]*30+["A"]*30,
                             "dictionary_fold": ["discovery"]*60+["validation"]*60})
        values = pd.DataFrame(rng.normal(size=(120, 2)), columns=["A", "F1"])
        values.loc[30:59, "A"] -= 5
        values.loc[30:59, "F1"] += 3
        qc = pd.DataFrame([dict(condition="DMSO", target="A", knockdown_remaining_ratio=.1,
                                eligible_both_conditions=True)])
        first = response_table(values, meta, ["A"], qc, "discovery")
        values.loc[60:] = rng.normal(0, 1e6, (60, 2))
        second = response_table(values, meta, ["A"], qc, "discovery")
        pd.testing.assert_frame_equal(first, second)
        self.assertTrue(first.iloc[0].eligible_seed)


if __name__ == "__main__":
    unittest.main()
