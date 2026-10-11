"""Program gene-set parsing, symbol mapping and target overlap on synthetic tables."""

import unittest

import numpy as np
import pandas as pd

from icpdag import gene_sets as gs


class Parsers(unittest.TestCase):
    def test_fufa_requires_up_and_sox10_bound(self):
        table = pd.DataFrame({"Gene Symbol": ["A", "B", "C", "D", None],
                              "Log2 Fold Change": [1.0, 2.0, -1.0, 0.5, 1.0],
                              "loss": [None] * 5, "gain": [None] * 5,
                              "SOX10": ["X", None, "X*", " x* ", "X"]})
        self.assertEqual(gs.parse_fufa(table), ["A", "D"])

    def test_kanai_up_and_down_totals(self):
        rows = [["Gene"] + [None] * 14, ["Symbol"] + [None] * 14]
        for name, up, down in [("U1", 4, 0), ("U2", 3, 0), ("D1", 0, 3)]:
            rows.append([name] + [None] * 12 + [up, down])
        out = gs.parse_kanai(pd.DataFrame(rows))
        self.assertEqual(out[gs.KANAI_UP], ["U1", "U2"])
        self.assertEqual(out[gs.KANAI_DOWN], ["D1"])
        rows.append(["BAD"] + [None] * 12 + [3, 3])
        with self.assertRaisesRegex(ValueError, "both up and down"):
            gs.parse_kanai(pd.DataFrame(rows))

    def test_tsoi_uses_own_labels_only(self):
        table = pd.DataFrame({"Gene": ["u", "unc", "nc", "nct", "t", "tm", "m"],
                              "Signature": list(gs.TSOI_LABELS)})
        out = gs.parse_tsoi(table)
        self.assertEqual(out["Tsoi_undifferentiated"], ["u"])
        self.assertEqual(out["Tsoi_neural_crest_like"], ["nc"])
        self.assertEqual(out["Tsoi_transitory"], ["t"])
        self.assertEqual(out["Tsoi_melanocytic"], ["m"])
        with self.assertRaisesRegex(ValueError, "Unknown"):
            gs.parse_tsoi(pd.DataFrame({"Gene": ["x"], "Signature": ["Other"]}))


class Mapping(unittest.TestCase):
    def genes(self):
        return pd.DataFrame({"gene_id": ["ENSG1.1", "ENSG2.4", "ENSG3.2", "ENSG4.1"],
                             "symbol": ["CCN1", "GRK3", "AAA", "BBB"]})

    def hgnc(self):
        return pd.DataFrame({"symbol": ["CCN1", "GRK3", "AAA", "BBB", "ZZZ"],
                             "prev_symbol": ["CYR61", "ADRBK2", None, None, None],
                             "alias_symbol": [None, None, "SHARED|ONLYA", "SHARED", "GONE"],
                             "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG9"]})

    def test_exact_previous_alias_ambiguous_missing(self):
        out = gs.map_symbols(["AAA", "CYR61", "ONLYA", "SHARED", "GONE", "NOPE"], self.genes(), self.hgnc())
        self.assertEqual(out.match_type.tolist(),
                         ["exact", "hgnc_previous", "hgnc_alias", "ambiguous", "missing", "missing"])
        self.assertEqual(out.symbol.tolist()[:3], ["AAA", "CCN1", "AAA"])
        self.assertEqual(out.gene_id.iloc[1], "ENSG1.1")

    def test_duplicates_are_flagged_and_dropped_for_scoring(self):
        sources = {name: [] for name in gs.ORDER}
        sources[gs.FUFA] = ["AAA", "ONLYA", "NOPE"]
        mapping = gs.mapping_table(sources, self.genes(), self.hgnc())
        fufa = mapping[mapping.program.eq(gs.FUFA)]
        self.assertEqual(fufa.duplicate_in_program.tolist(), [False, True, False])
        programs = gs.programs_for_scoring(mapping)
        self.assertEqual(programs[gs.FUFA]["genes"], ["AAA", "NOPE"])
        self.assertNotIn("role", programs[gs.FUFA])

    def test_scoring_genes_drop_par_y_and_prefer_protein_coding(self):
        var = pd.DataFrame({"symbol": ["A", "A", "B", "B"],
                            "gene_type": ["lncRNA", "protein_coding", "protein_coding", "protein_coding"],
                            "is_par_y": [False, False, False, True]},
                           index=pd.Index(["ENSG2.1", "ENSG1.1", "ENSG3.1", "ENSG3.1_PAR_Y"], name="gene_id"))
        np.testing.assert_array_equal(gs.scoring_genes(var), [False, True, True, False])

    def test_target_overlap(self):
        programs = {gs.FUFA: {"genes": ["PTEN", "MITF"]}, gs.KANAI_UP: {"genes": ["CCN1"]}}
        out = gs.target_overlap(programs, {"PTEN": "PTEN", "NF2": "NF2"})
        self.assertEqual(out.to_dict("records"), [dict(target="PTEN", program=gs.FUFA, target_symbol="PTEN")])


if __name__ == "__main__":
    unittest.main()
