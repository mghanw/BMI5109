"""Phase 1 preparation: construct parsing, assignment, knockdown QC and output safety."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prepare_data", ROOT / "scripts" / "prepare_data.py")
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)


def config():
    return prep.load_config(ROOT / "config" / "prepare_data.json")


class ConstructsAndGenes(unittest.TestCase):
    def genes(self):
        return pd.DataFrame({"gene_id": ["ENSG01.3", "ENSG02.1", "ENSG03.7", "ENSG04.2", "ENSG04.2_PAR_Y"],
                             "symbol": ["AAA", "BBB", "ATP5MC3", "DDD", "DDD"],
                             "gene_type": ["protein_coding"] * 5, "chrom": ["chr1"] * 4 + ["chrY"]})

    def test_construct_names_use_global_index(self):
        self.assertEqual(prep.parse_construct("ACTN4_1"), "ACTN4")
        self.assertEqual(prep.parse_construct("NO-TARGET_144"), "NTC")
        self.assertEqual(prep.parse_construct("C7orf26_11"), "C7orf26")
        with self.assertRaises(ValueError):
            prep.parse_construct("ACTN4_P1")

    def test_targets_resolve_by_id_symbol_and_alias(self):
        library = pd.DataFrame({"id": ["AAA_1", "OLD_2", "BBB_3", "ATP5G3_4", "NO-TARGET_5", "NO-TARGET_6"],
                                "sgRNA_pair_id": ["1_AAA_P1_ENSG01", "2_OLD_P1_ENSG04", "BBB_+_1.23-P1",
                                                  "ATP5G3_-_1.23-P1", "x", "y"]})
        genes = self.genes()
        genes.loc[2, "gene_id"] = "ENSG00000154518.10"
        out = prep.resolve_targets(library, genes).set_index("guide")
        self.assertEqual(out.loc["AAA_1", "source"], "library_ensembl_id")
        self.assertEqual(out.loc["OLD_2", "target_gene_id"], "ENSG04.2")  # PAR_Y copy is never chosen
        self.assertEqual(out.loc["BBB_3", "source"], "symbol")
        self.assertEqual(out.loc["ATP5G3_4", "target_symbol"], "ATP5MC3")
        self.assertEqual(out.loc["NO-TARGET_6", "target"], "NTC")

    def test_unresolvable_target_fails(self):
        library = pd.DataFrame({"id": ["ZZZ_1"], "sgRNA_pair_id": ["ZZZ_+_1.23-P1"]})
        with self.assertRaises(ValueError):
            prep.resolve_targets(library, self.genes())

    def test_alias_table_is_unique(self):
        ids = [v[1] for v in prep.TARGET_ALIASES.values()]
        self.assertEqual(len(ids), len(set(ids)))

    def test_missing_matrix_gene_fails(self):
        with self.assertRaises(ValueError):
            prep.annotate_genes(["ENSG01.3", "ENSG01.4"], self.genes())
        out = prep.annotate_genes(["ENSG04.2_PAR_Y", "ENSG01.3"], self.genes())
        self.assertEqual(out.is_par_y.tolist(), [True, False])


class Assignment(unittest.TestCase):
    def test_barcode_split(self):
        whitelist = ["AAAAAAAAAA", "CCCCCCCCCC"]
        out = prep.parse_barcode(["240613_2lig_RNA_a05.CCCCCCCCCCGGGGGGGGGG"], whitelist)
        self.assertEqual(out.iloc[0].tolist(), ["240613_2lig_RNA_a05", "a05", 1])
        with self.assertRaises(ValueError):
            prep.parse_barcode(["p.TTTTTTTTTTGGGGGGGGGG"], whitelist)

    def test_threshold_at_valley_of_bimodal_counts(self):
        rng = np.random.default_rng(0)
        ambient = rng.integers(1, 4, 5000)
        real = np.round(np.exp(rng.normal(np.log(40), 0.6, 5000))).astype(int)
        self.assertIn(prep.choose_threshold(np.r_[ambient, real], (3, 10)), range(4, 9))

    def test_classes(self):
        X = sp.csr_matrix(np.array([[50, 1, 0], [30, 20, 0], [2, 1, 0], [6, 1, 1], [40, 9, 0], [8, 8, 0],
                                    [5, 0, 0], [4, 1, 0]]))
        cls, n = prep.classify_threshold(X, 5)
        self.assertEqual(cls.tolist(), ["singlet", "multiplet", "unassigned_low", "singlet", "multiplet", "multiplet",
                                        "singlet", "unassigned_low"])
        total, top, arg, second = prep.top_two(X)
        self.assertEqual(second.tolist(), [1, 20, 1, 1, 9, 8, 0, 1])
        paper = prep.classify_paper(total, top, second, 5, 0.5, 0.3)
        self.assertEqual(paper.tolist(), ["singlet", "multiplet", "unassigned_low", "singlet", "singlet", "multiplet",
                                          "unassigned_low", "unassigned_low"])
        inclusive = prep.classify_paper(total, top, second, 5, 0.5, 0.3, inclusive=True)
        self.assertEqual(inclusive.tolist()[6:], ["singlet", "singlet"])
        self.assertEqual(inclusive.tolist()[:6], paper.tolist()[:6])

    def test_default_rule_is_strict_paper_rule(self):
        a = config()["assignment"]
        self.assertEqual(a["rule"], "paper")
        self.assertFalse(a["paper_min_total_inclusive"])
        rules = prep.rule_columns(config(), 5)
        self.assertEqual(list(rules.items())[0], ("class_paper", "paper_gt5_max0.5_sec0.3"))
        self.assertEqual(rules["class_paper_alt"], "paper_ge5_max0.5_sec0.3")

    def test_cells_without_sgrna_are_no_sgrna(self):
        cfg = config()
        cfg["assignment"]["threshold"] = 5
        X = sp.csr_matrix(np.array([[50, 1], [0, 30], [4, 1]]))
        out, T, _ = prep.assign_cells(["a", "b", "c", "d"], ["a", "c", "d"], X, ["G_1", "NO-TARGET_2"], cfg)
        self.assertEqual(out.assignment.tolist(), ["singlet", "no_sgrna", "singlet", "unassigned_low"])
        self.assertTrue(out.assignment.equals(out.class_paper))
        self.assertEqual(out.class_paper_alt.tolist(), ["singlet", "no_sgrna", "singlet", "singlet"])
        self.assertEqual(out.guide.tolist(), ["G_1", None, "NO-TARGET_2", None])

    def assigned(self):
        obs = pd.DataFrame({"condition": ["DMSO", "DMSO", "DMSO", "Vem", "Vem"],
                            "top_guide": ["G_1", "G_1", "NO-TARGET_2", "G_1", "G_1"],
                            "class_paper": ["singlet", "multiplet", "singlet", "singlet", "unassigned_low"],
                            "class_paper_alt": ["singlet", "multiplet", "singlet", "singlet", "singlet"],
                            "class_threshold": ["singlet", "singlet", "singlet", "singlet", "singlet"]})
        return obs, prep.rule_columns(config(), 5)

    def test_assignment_summary_marks_one_primary_rule(self):
        obs, rules = self.assigned()
        summ = prep.assignment_summary(obs[list(rules)], obs.condition, rules, {3: np.array(["singlet"] * 5)},
                                       pd.Index(range(5)))
        roles = summ.drop_duplicates("rule").set_index("rule").role
        self.assertEqual(roles[roles.eq("primary")].index.tolist(), ["paper_gt5_max0.5_sec0.3"])
        self.assertEqual(roles[["paper_ge5_max0.5_sec0.3", "threshold_T5_auto", "threshold_T3"]].tolist(),
                         ["sensitivity"] * 3)
        self.assertEqual(summ.rule.iloc[0], "paper_gt5_max0.5_sec0.3")

    def test_report_comparison(self):
        obs, rules = self.assigned()
        comp = prep.report_comparison(obs, pd.Series({"G_1": "G", "NO-TARGET_2": "NTC"}),
                                      {"G": [1, 2], "NTC": [1, 0]}, rules).set_index(["target", "condition"])
        self.assertEqual(comp.loc[("G", "DMSO"), "threshold_T5_auto_n"], 2)
        self.assertEqual(comp.loc[("G", "Vem"), "paper_gt5_max0.5_sec0.3_n"], 1)
        self.assertEqual(comp.loc[("G", "Vem"), "paper_ge5_max0.5_sec0.3_n"], 2)
        self.assertEqual(comp.matches_report.tolist(), [True, False, True, True])
        self.assertEqual(comp.loc[("G", "Vem"), "primary_minus_report"], -1)

    def test_calibration_folds_keep_wells_whole_and_balanced(self):
        wells = np.repeat(np.arange(8), 5)
        condition = np.where(wells < 4, "DMSO", "Vem")
        fold = prep.calibration_folds(wells, condition, 2, 1)
        frame = pd.DataFrame({"w": wells, "c": condition, "f": fold})
        self.assertTrue(frame.groupby("w").f.nunique().eq(1).all())
        self.assertTrue(frame.drop_duplicates("w").groupby("c").f.sum().eq(2).all())


class Knockdown(unittest.TestCase):
    def table(self, n_kd=60):
        rng = np.random.default_rng(1)
        obs, cols = [], []
        for condition in ["DMSO", "Vem"]:
            for target, scale in [("NTC", None), ("A", 0.2), ("B", 0.9), ("C", 0.2)]:
                n = 200 if target == "NTC" else n_kd
                obs += [dict(condition=condition, target=target)] * n
                base = rng.poisson(10, (n, 3)).astype(float)
                if scale is not None:
                    j = "ABC".index(target)
                    base[:, j] = rng.poisson(10 * scale, n)
                cols.append(base)
        cp = np.vstack(cols)
        cp[:, 2] = 0  # target C: never detected in NTC
        obs = pd.DataFrame(obs)
        return prep.knockdown_table(cp, obs, {"A": "gA", "B": "gB", "C": "gC"}, config())

    def test_ratio_status_and_eligibility(self):
        qc = self.table().set_index(["condition", "target"])
        self.assertAlmostEqual(qc.loc[("DMSO", "A"), "knockdown_remaining_ratio"], 0.2, delta=0.05)
        self.assertEqual(qc.loc[("DMSO", "A"), "kd_status"], "adequate")
        self.assertEqual(qc.loc[("DMSO", "B"), "kd_status"], "inadequate")
        self.assertEqual(qc.loc[("Vem", "C"), "kd_status"], "not_evaluable")
        self.assertTrue(qc.loc[("DMSO", "A"), "eligible"])
        self.assertFalse(qc.loc[("DMSO", "B"), "eligible"])
        self.assertFalse(qc.loc[("DMSO", "C"), "eligible"])
        self.assertEqual(qc.loc[("DMSO", "NTC"), "n_cells"], 200)
        self.assertIn("eligible_paper_kd_rule", qc.columns)

    def test_min_cells_applies_to_each_condition(self):
        qc = self.table(n_kd=50).set_index(["condition", "target"])
        self.assertFalse(qc.loc[("DMSO", "A"), "eligible"])  # 50 is not > 50

    def test_bh(self):
        q = prep.bh([0.01, 0.04, np.nan, 0.03])
        np.testing.assert_allclose(q[[0, 1, 3]], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(q[2]))


class ConfigAndSafety(unittest.TestCase):
    def test_config_override(self):
        cfg = prep.load_config(ROOT / "config" / "prepare_data.json", {"knockdown.max_ratio": 0.4})
        self.assertEqual(cfg["knockdown"]["max_ratio"], 0.4)
        self.assertEqual(config()["knockdown"]["max_ratio"], 0.5)
        with self.assertRaises(KeyError):
            prep.load_config(ROOT / "config" / "prepare_data.json", {"knockdown.typo": 1})

    def test_aggregate_outputs_reject_barcodes(self):
        prep.assert_no_cell_ids(pd.DataFrame({"well": ["L1_01"], "target": ["NF2"]}))
        with self.assertRaises(ValueError):
            prep.assert_no_cell_ids(pd.DataFrame({"cell": ["240613_2lig_RNA_a05.ACGTACGTACGTACGTACGT"]}))

    def test_batch_check_detects_planted_split(self):
        rng = np.random.default_rng(3)
        n, g = 48, 300
        counts = rng.poisson(50, (n, g)).astype(float)
        counts[:24, :30] *= 3  # contiguous-halves batch
        info = pd.DataFrame({"lig1_index": np.arange(n), "condition": "DMSO", "median_umi": rng.normal(3000, 300, n)})
        cfg = config()
        cfg["batch_check"].update(n_hvg=100, n_permutations=200)
        _, tests = prep.well_batch_check(counts, info, cfg, rng)
        tests = tests.set_index("split")
        self.assertLess(tests.loc["contiguous_halves", "perm_p"], 0.01)
        self.assertGreater(tests.loc["parity", "perm_p"], 0.01)


if __name__ == "__main__":
    unittest.main()
