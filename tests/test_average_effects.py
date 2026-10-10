"""Factorial contrasts must preserve drug effects and cancel shared controls."""

import unittest

import numpy as np
import pandas as pd

from icpdag.average_effects import define_contrasts, estimate_contrasts, reference_sensitivity


class AverageEffectContracts(unittest.TestCase):
    def data(self):
        metadata, values = [], []
        for condition, drug in [("DMSO", 0), ("Vem", 10)]:
            for target, gene, interaction in [("NTC", 0, 0), ("A", 2, -3), ("B", -1, 0)]:
                noise = np.linspace(-1, 1, 20) * (100 if target == "NTC" else 1)
                for i, error in enumerate(noise):
                    metadata.append(dict(condition=condition, target=target, guide=target))
                    values.append(gene + drug + (interaction if drug else 0) + error)
        return pd.DataFrame({"readout": values}), pd.DataFrame(metadata)

    def run_estimation(self):
        scores, meta = self.data()
        contrasts = define_contrasts(["A", "B"], [["A", "B"]])
        effects, _, groups = estimate_contrasts(scores, meta, contrasts, 400, 7)
        return effects, groups

    def test_factorial_effects_and_pair_interaction(self):
        effects, _ = self.run_estimation()
        actual = effects.set_index(["kind", "target", "condition"]).estimate
        self.assertAlmostEqual(actual["gene", "A", "DMSO"], 2)
        self.assertAlmostEqual(actual["gene", "A", "Vem"], -1)
        self.assertAlmostEqual(actual["drug", "A", "Vem_minus_DMSO"], 7)
        self.assertAlmostEqual(actual["drug", "B", "Vem_minus_DMSO"], 10)
        self.assertAlmostEqual(actual["interaction", "A", "Vem_minus_DMSO"], -3)
        self.assertAlmostEqual(actual["interaction", "B", "Vem_minus_DMSO"], 0)
        self.assertAlmostEqual(actual["pair_interaction", "A", "Vem_minus_DMSO"], -3)

    def test_pair_comparison_cancels_noisy_shared_control(self):
        effects, groups = self.run_estimation()
        pair = effects.loc[effects.kind.eq("pair") & effects.condition.eq("DMSO")].iloc[0]
        expected = np.sqrt(sum(groups["DMSO", t]["variance"][0] / groups["DMSO", t]["n"] for t in ("A", "B")))
        self.assertAlmostEqual(pair.iid_cell_se, expected)
        self.assertLess(pair.ci_high - pair.ci_low, 1)
        self.assertNotIn("NTC", pair.group_cell_counts)

    def test_metadata_alignment_and_invalid_cells(self):
        scores, meta = self.data()
        contrasts = define_contrasts(["A", "B"], [])
        a, _, _ = estimate_contrasts(scores, meta, contrasts, 100, 4)
        b, _, _ = estimate_contrasts(scores, meta.iloc[::-1], contrasts, 100, 4)
        pd.testing.assert_frame_equal(a, b)
        with self.assertRaisesRegex(ValueError, "identical cells"):
            estimate_contrasts(scores.iloc[:-1], meta, contrasts, 100, 4)

    def test_reference_weights_do_not_confuse_guide_and_cell_average(self):
        values, metadata = [], []
        for condition in ("DMSO", "Vem"):
            for target, guide, n, mean in [("NTC", "c1", 2, 1), ("NTC", "c2", 6, 3), ("A", "a1", 4, 4)]:
                for _ in range(n):
                    metadata.append(dict(condition=condition, target=target, guide=guide))
                    values.append(mean)
        scores, meta = pd.DataFrame({"readout": values}), pd.DataFrame(metadata)
        contrasts = define_contrasts(["A"], [])
        effects, _, groups = estimate_contrasts(scores, meta, contrasts, 100, 1)
        contrast = effects.loc[effects.kind.eq("gene") & effects.condition.eq("DMSO")].iloc[0]
        self.assertAlmostEqual(contrast.estimate, 1.5)
        sensitivity, _, _ = reference_sensitivity(scores, meta, contrasts, groups)
        equal = sensitivity.loc[sensitivity.contrast_id.eq(contrast.contrast_id)
                                & sensitivity.reference_scheme.eq("equal_construct_weight")].iloc[0]
        self.assertAlmostEqual(equal.estimate, 2)


if __name__ == "__main__":
    unittest.main()
