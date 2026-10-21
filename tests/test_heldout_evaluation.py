"""Held-out evaluation: statistics with known answers, NaN rules, verification and leakage guards (synthetic)."""

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from icpdag import heldout_evaluation as ev
from icpdag.dispersion import centered_wasserstein1, iqr_diff, log_var_ratio
from icpdag.sem_prediction import MODELS

ROOT = Path(__file__).resolve().parents[1]
CLIP = 0.999999999999


class PearsonFisher(unittest.TestCase):
    def test_clip_and_nan_rule(self):
        x = np.arange(10.0)
        self.assertAlmostEqual(ev.pearson(x, 2 * x + 1), 1.0)
        self.assertTrue(np.isfinite(ev.fisher_z(1.0, CLIP)))
        self.assertAlmostEqual(float(ev.fisher_z(1.0, CLIP)), np.arctanh(CLIP))
        self.assertAlmostEqual(float(ev.fisher_z(-1.0, CLIP)), -np.arctanh(CLIP))
        self.assertTrue(np.isnan(ev.pearson(np.zeros(10), x)))
        self.assertTrue(np.isnan(ev.pearson(x, np.full(10, 3.0))))
        self.assertTrue(np.isnan(ev.fisher_z(np.nan, CLIP)))
        m = ev.direction_metrics(np.zeros(10), x, CLIP)
        self.assertTrue(m["constant_prediction"] and np.isnan(m["z"]))
        self.assertAlmostEqual(m["rmse"], np.sqrt(np.mean(x ** 2)))

    def test_constant_tolerance(self):
        self.assertTrue(ev.is_constant(np.full(5, 1.0) + 1e-13 * np.arange(5)))
        self.assertFalse(ev.is_constant(1e-10 * np.arange(5)))
        self.assertFalse(ev.is_constant(np.array([0.0, np.nan])))


def metric_rows(spec):
    rows = []
    for (model, target, condition, direction), (z, cp, co) in spec.items():
        rows.append(dict(model=model, target=target, condition=condition, direction=direction, r=np.tanh(z), z=z,
                         rmse=float(direction + 1), constant_prediction=cp, constant_observed=co))
    return pd.DataFrame(rows)


class Averaging(unittest.TestCase):
    def setUp(self):
        self.rows = metric_rows({
            ("a", "T1", "DMSO", 0): (0.2, False, False), ("a", "T1", "DMSO", 1): (0.4, False, False),
            ("a", "T1", "Vem", 0): (1.0, False, False), ("a", "T1", "Vem", 1): (0.0, False, False),
            ("a", "T2", "DMSO", 0): (np.nan, True, False), ("a", "T2", "DMSO", 1): (np.nan, True, False),
            ("a", "T2", "Vem", 0): (0.6, False, False), ("a", "T2", "Vem", 1): (0.2, False, False),
            ("a", "T3", "DMSO", 0): (np.nan, True, True), ("a", "T3", "DMSO", 1): (0.3, False, False),
            ("a", "T3", "Vem", 0): (0.1, False, False), ("a", "T3", "Vem", 1): (0.1, False, False),
        })
        self.out = ev.combine_metrics(self.rows)

    def get(self, target, condition, direction="combined"):
        o = self.out
        return o[o.target.eq(target) & o.condition.eq(condition) & o.direction.eq(direction)].iloc[0]

    def test_direction_then_condition_on_z_scale(self):
        self.assertAlmostEqual(self.get("T1", "DMSO").z, 0.3)
        self.assertAlmostEqual(self.get("T1", "Vem").z, 0.5)
        self.assertAlmostEqual(self.get("T1", "pooled").z, 0.4)
        self.assertAlmostEqual(self.get("T1", "pooled").r, np.tanh(0.4))
        self.assertAlmostEqual(self.get("T1", "pooled").rmse, 1.5)
        self.assertAlmostEqual(self.get("T1", "DMSO", 0).r, np.tanh(0.2))

    def test_r_zero_substitution_vs_complete_pairs(self):
        t2 = self.get("T2", "pooled")
        self.assertTrue(np.isnan(t2.z))                        # pre-registered: NaN propagates
        self.assertAlmostEqual(t2.z_amended, (0.0 + 0.4) / 2)   # amended: constant prediction -> z = 0
        self.assertTrue(bool(t2.constant_prediction))

    def test_constant_observed_is_nan_in_every_version(self):
        t3 = self.get("T3", "pooled")
        self.assertTrue(np.isnan(t3.z) and np.isnan(t3.z_amended))

    def test_paired_test_drops_incomplete_pairs(self):
        a = np.array([0.5, np.nan, 0.2, 0.4, 0.3])
        b = np.array([0.1, 0.2, np.nan, 0.1, 0.1])
        out = ev.paired_test(a, b, 200, 1)
        self.assertEqual(out["n_targets"], 5)
        self.assertEqual(out["n_nan_dropped"], 2)
        self.assertEqual(out["n_effective"], 3)


class Wilcoxon(unittest.TestCase):
    def test_zero_differences_are_dropped_and_counted(self):
        d = np.array([0.0, 0.0, 0.3, -0.1, 0.5, 0.7, -0.2, 0.9, 1.1, 0.4])
        out = ev.wilcoxon_test(d)
        self.assertEqual(out["n_zero_dropped"], 2)
        self.assertEqual(out["n_effective"], 8)
        self.assertEqual(out["method"], "approx")
        ref = stats.wilcoxon(d[d != 0], zero_method="wilcox", method="asymptotic")
        self.assertAlmostEqual(out["p_value"], ref.pvalue)

    def test_method_rule(self):
        self.assertEqual(ev.wilcoxon_method(np.arange(1, 11) * 0.1), "exact")
        self.assertEqual(ev.wilcoxon_method(np.array([0.1, -0.1, 0.3])), "approx")      # tie in |d|
        self.assertEqual(ev.wilcoxon_method(np.array([0.0, 0.2, 0.3])), "approx")       # zero difference
        self.assertEqual(ev.wilcoxon_method(np.arange(1, 52) * 0.1), "approx")          # n_eff > 50
        d = np.arange(1, 11) * 0.1 - 0.35
        out = ev.wilcoxon_test(d)
        self.assertEqual(out["method"], "exact")
        self.assertAlmostEqual(out["p_value"], stats.wilcoxon(d, method="exact").pvalue)

    def test_all_zero(self):
        out = ev.wilcoxon_test(np.zeros(4))
        self.assertEqual(out["n_effective"], 0)
        self.assertTrue(np.isnan(out["p_value"]))


class HodgesLehmann(unittest.TestCase):
    def test_known_samples(self):
        self.assertAlmostEqual(ev.hodges_lehmann([1.0, 2.0, 3.0]), 2.0)
        # Walsh averages of (0, 0, 4): 0, 0, 2, 0, 2, 4 -> median 1; zeros are included
        self.assertAlmostEqual(ev.hodges_lehmann([0.0, 0.0, 4.0]), 1.0)
        self.assertAlmostEqual(ev.hodges_lehmann([5.0]), 5.0)
        self.assertTrue(np.isnan(ev.hodges_lehmann([])))

    def test_bootstrap_ci_is_deterministic(self):
        d = np.random.default_rng(3).normal(0.2, 1, 40)
        a = ev.hodges_lehmann_ci(d, 1000, 20261011)
        self.assertEqual(a, ev.hodges_lehmann_ci(d, 1000, 20261011))
        self.assertNotEqual(a, ev.hodges_lehmann_ci(d, 1000, 7))
        self.assertLess(a[0], ev.hodges_lehmann(d))
        self.assertGreater(a[1], ev.hodges_lehmann(d))
        # chunking does not change the draws
        self.assertEqual(a, ev.hodges_lehmann_ci(d, 1000, 20261011, chunk=37))


class Holm(unittest.TestCase):
    def test_known_values(self):
        np.testing.assert_allclose(ev.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
        np.testing.assert_allclose(ev.holm([0.5, 0.9]), [1.0, 1.0])
        out = ev.holm([0.01, np.nan, 0.04])
        np.testing.assert_allclose(out[[0, 2]], [0.02, 0.04])
        self.assertTrue(np.isnan(out[1]))


class NoiseCeiling(unittest.TestCase):
    def test_spearman_brown(self):
        self.assertAlmostEqual(float(ev.spearman_brown(0.5)), 2 / 3)
        self.assertAlmostEqual(float(ev.spearman_brown(0.0)), 0.0)
        self.assertAlmostEqual(float(ev.spearman_brown(1.0)), 1.0)

    def test_split_half(self):
        rng = np.random.default_rng(0)
        signal = rng.normal(size=30) * 3
        ntc = rng.normal(size=(200, 30))
        target = rng.normal(size=(120, 30)) + signal
        a = ev.split_half_ceiling(target, ntc, 20, ev.noise_ceiling_rng(1, "DMSO", "T", 0), CLIP)
        b = ev.split_half_ceiling(target, ntc, 20, ev.noise_ceiling_rng(1, "DMSO", "T", 0), CLIP)
        np.testing.assert_array_equal(a, b)
        self.assertGreater(np.median(a), 0.9)
        null = ev.split_half_ceiling(rng.normal(size=(120, 30)), ntc, 50, np.random.default_rng(2), CLIP)
        self.assertLess(abs(np.median(null)), 0.5)
        self.assertTrue(np.all(np.abs(null) <= CLIP))
        self.assertAlmostEqual(ev.combine_directions_z([0.5, 0.7], CLIP),
                               np.tanh((np.arctanh(0.5) + np.arctanh(0.7)) / 2))


def toy_meta():
    rows = []
    for target, n in (("NTC", 8), ("T1", 6)):
        for fold in (0, 1):
            rows += [dict(condition="DMSO", target=target, calib_fold=fold)] * n
    return pd.DataFrame(rows)


class ObservedResponse(unittest.TestCase):
    def test_uses_only_calib_fold_one_minus_c(self):
        meta = toy_meta()
        values = np.zeros((len(meta), 3))
        t, n = meta.target.eq("T1").to_numpy(), meta.target.eq("NTC").to_numpy()
        f0, f1 = meta.calib_fold.eq(0).to_numpy(), meta.calib_fold.eq(1).to_numpy()
        values[t & f1] = [1.0, 2.0, 3.0]
        values[t & f0] = 1000.0                                 # calibration cells of direction 0
        values[n & f1] = [0.5, 0.5, 0.5]
        np.testing.assert_allclose(ev.observed_change(values, meta, "T1", "DMSO", 0), [0.5, 1.5, 2.5])
        np.testing.assert_allclose(ev.observed_change(values, meta, "T1", "DMSO", 1), [1000.0, 1000.0, 1000.0])
        score = values[:, 0]
        self.assertAlmostEqual(ev.observed_mean_contrast(score, meta, "T1", "DMSO", 0), 0.5)   # target - NTC
        with self.assertRaises(ValueError):
            ev.observed_change(values, meta, "T9", "DMSO", 0)


class Normalisation(unittest.TestCase):
    def test_log_cp10k_paths_are_identical(self):
        sys.path.insert(0, str(ROOT / "scripts"))
        import aim1
        import anndata as ad
        from icpdag.cscn_panel import load_ntc_expression, log_cp10k
        from icpdag.sem_inputs import load_panel_cells

        rng = np.random.default_rng(5)
        counts = sp.csr_matrix(rng.poisson(0.8, (40, 12)).astype(np.float32))
        total = np.asarray(counts.sum(axis=1)).ravel() + 50.0
        np.testing.assert_array_equal(aim1.log_cp10k(counts, total).toarray(), log_cp10k(counts, total).toarray())
        obs = pd.DataFrame(dict(condition=["DMSO", "Vem"] * 20, target=(["NTC"] * 3 + ["T1"]) * 10,
                                calib_fold=[0, 1, 1, 0] * 10, lig1_well=["L1_1"] * 40, total_counts=total),
                           index=[f"cell{i}" for i in range(40)])
        var = pd.DataFrame(index=[f"G{i}" for i in range(12)])
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "toy.h5ad"
            ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(path)
            genes = ["G3", "G0", "G7"]
            panel, meta = load_panel_cells(path, genes, ["DMSO"], ["T1"])
            ntc, ntc_meta = load_ntc_expression(path, genes, ["DMSO"])
        rows = meta.target.eq("NTC").to_numpy()
        np.testing.assert_array_equal(panel[rows], ntc)
        self.assertEqual(list(meta.index[rows]), list(ntc_meta.index))
        expected = log_cp10k(counts[np.flatnonzero(obs.condition.eq("DMSO") & obs.target.isin(["NTC", "T1"]))][:, [3, 0, 7]],
                             total[(obs.condition.eq("DMSO") & obs.target.isin(["NTC", "T1"])).to_numpy()])
        np.testing.assert_array_equal(panel, expected.toarray())


class FufaContrasts(unittest.TestCase):
    def test_sign_and_interaction(self):
        rows = ev.fufa_contrast_rows({"DMSO": [1.0, 1.0], "Vem": [3.0, 3.0]}, {"DMSO": [0.0, 2.0], "Vem": [1.0, 1.0]})
        get = {(r["contrast"], r["direction"]): r for r in rows}
        self.assertAlmostEqual(get["KD_minus_NTC_DMSO", "combined"]["abs_error"], 1.0)
        self.assertAlmostEqual(get["KD_minus_NTC_Vem", 0]["abs_error"], 2.0)
        inter = get["interaction_Vem_minus_DMSO", 1]
        self.assertAlmostEqual(inter["predicted"], 2.0)
        self.assertAlmostEqual(inter["observed"], -1.0)
        self.assertAlmostEqual(get["interaction_Vem_minus_DMSO", "combined"]["abs_error"], 2.0)   # mean of 1 and 3


class Tertiary(unittest.TestCase):
    def setUp(self):
        self.x = np.random.default_rng(4).normal(size=500)

    def test_known_distributions(self):
        s, n_drop = ev.predicted_distribution_stats(self.x, np.full(500, 0.7))
        self.assertEqual(n_drop, 0)
        self.assertAlmostEqual(s["log_var_ratio"], 0.0)
        self.assertAlmostEqual(s["iqr_diff"], 0.0, places=12)
        self.assertAlmostEqual(s["centered_wasserstein"], 0.0, places=12)
        shift = (self.x - self.x.mean()) * 1.0                   # doubles every deviation
        s, _ = ev.predicted_distribution_stats(self.x, shift)
        self.assertAlmostEqual(s["log_var_ratio"], np.log(4.0))
        self.assertAlmostEqual(s["iqr_diff"], iqr_diff(2 * self.x, self.x))
        self.assertAlmostEqual(s["centered_wasserstein"], centered_wasserstein1(2 * self.x, self.x))
        o = ev.distribution_stats(2 * self.x, self.x)
        self.assertAlmostEqual(o["log_var_ratio"], log_var_ratio(2 * self.x, self.x))

    def test_dropped_cells_use_the_same_subset(self):
        shift = np.full(500, 0.3)
        shift[self.x > 1] = np.nan                                # drop the upper tail
        s, n_drop = ev.predicted_distribution_stats(self.x, shift)
        self.assertEqual(n_drop, int((self.x > 1).sum()))
        self.assertAlmostEqual(s["log_var_ratio"], 0.0)            # no artificial dispersion difference
        self.assertAlmostEqual(s["centered_wasserstein"], 0.0, places=12)

    def test_interval_and_flags(self):
        y = np.random.default_rng(8).normal(size=300)
        x = np.random.default_rng(9).normal(size=100) * 1.5
        a = ev.observed_lvr_interval(x, y, 500, 1, "DMSO", "T1", 1)
        self.assertEqual(a, ev.observed_lvr_interval(x, y, 500, 1, "DMSO", "T1", 1))
        self.assertLess(a[0], log_var_ratio(x, y))
        self.assertGreater(a[1], log_var_ratio(x, y))
        self.assertEqual(ev.calibration_flags(2.0, -1.0, 1.0), dict(over_dispersion=True, covered=False))
        self.assertEqual(ev.calibration_flags(0.0, -1.0, 1.0), dict(over_dispersion=False, covered=True))


# ---------------------------------------------------------------- prediction verification
NODES = ["A", "B", "C", "D"]


def toy_predictions(n_ref=6, seed=0):
    rng = np.random.default_rng(seed)
    M = len(MODELS)
    preds = {}
    for condition in ("DMSO", "Vem"):
        T = 2
        per_unit = rng.normal(size=(M, T, len(NODES)))
        per_unit[MODELS.index("zero_downstream_baseline")] = 0.0
        delta = rng.normal(size=(T, 2))
        scaled = per_unit[:, :, None, :] * delta[None, :, :, None]
        fufa = rng.normal(size=(M, n_ref, T, 2))
        fufa[MODELS.index("shared_cscn_sem")] = 0.25
        preds[condition] = dict(targets=["T1", "T2"], per_unit=per_unit, scaled=scaled, delta=delta,
                                included=np.full((M, T), n_ref), fufa_scaled=fufa)
    return preds


class PredictionChecks(unittest.TestCase):
    def test_constant_shift_assertion(self):
        preds = toy_predictions()
        preds["DMSO"]["fufa_scaled"][MODELS.index("shared_cscn_sem"), 0, 0, 0] = np.nan   # excluded cell is fine
        ev.check_constant_shift(preds, MODELS, ["shared_cscn_sem"])
        with self.assertRaises(ValueError):
            ev.check_constant_shift(preds, MODELS, ["cell_specific_cscn_sem"])

    def test_constant_flags_against_structure_report(self):
        preds = toy_predictions()
        preds["Vem"]["per_unit"][MODELS.index("shared_cscn_sem"), 1, :] = [0.0, 0.0, 0.0, 1.0]   # only D changes
        tn = {"T1": "A", "T2": "D"}
        flags = ev.constant_prediction_flags(preds, MODELS, NODES, tn)
        rows = []
        for gt, model in ev.STRUCTURE_MODELS.items():
            for condition in ("DMSO", "Vem"):
                for target in ("T1", "T2"):
                    const = flags[(flags.model == model) & (flags.condition == condition)
                                  & (flags.target == target)].constant.iloc[0]
                    rows.append(dict(graph_type=gt, condition=condition, target=target,
                                     constant_prediction_expected=bool(const)))
        structure = pd.DataFrame(rows)
        counts = ev.check_constant_flags(flags, structure)
        self.assertEqual(counts.loc["shared_cscn_sem", "any_condition"], 1)
        self.assertEqual(counts.loc["zero_downstream_baseline", "any_condition"], 2)
        structure.loc[0, "constant_prediction_expected"] = True
        with self.assertRaises(ValueError):
            ev.check_constant_flags(flags, structure)

    def test_included_zero_and_fufa_pseudobulk(self):
        preds = toy_predictions()
        self.assertEqual(len(ev.included_zero_rows(preds, MODELS)), 0)
        preds["DMSO"]["included"][2, 1] = 0
        self.assertEqual(ev.included_zero_rows(preds, MODELS).to_dict("records"),
                         [dict(model=MODELS[2], condition="DMSO", target="T2")])
        w = np.array([0.5, -0.25, 0.0, 1.0])
        for p in preds.values():
            p["fufa_scaled"] = np.repeat((p["scaled"] @ w)[:, None], 6, axis=1)
        ev.check_fufa_pseudobulk(preds, w)
        preds["Vem"]["fufa_scaled"][0, 0, 0, 0] += 1.0
        with self.assertRaises(ValueError):
            ev.check_fufa_pseudobulk(preds, w)


def write_fold(base, condition, fold, targets, config_hash, n_ref=5, directions=2, lam=10.0):
    d = base / condition / f"fold{fold}"
    shard = d / "shard_000_of_001"
    shard.mkdir(parents=True)
    refs = np.arange(n_ref)
    config = dict(condition=condition, fold=fold, shard=0, n_shards=1, targets=targets, lam=lam,
                  reference_sha256=ev.index_digest(refs), chunk_size=8, config_sha256=config_hash)
    (shard / "config.json").write_text(json.dumps(config), encoding="utf-8")
    M, T, P = len(MODELS), len(targets), len(NODES)
    np.savez_compressed(d / "predictions.npz", models=np.asarray(MODELS), targets=np.asarray(targets),
                        nodes=np.asarray(NODES), per_unit=np.zeros((M, T, P)), scaled=np.zeros((M, T, directions, P)),
                        delta=np.zeros((T, 2)), included=np.ones((M, T)), fufa_per_unit=np.zeros((M, n_ref, T)),
                        fufa_scaled=np.zeros((M, n_ref, T, 2)), reference_cells=refs, lam=lam)


class VerifyPredictions(unittest.TestCase):
    split = pd.DataFrame(dict(target=["T1", "T2", "T3", "T4"], fold=[0, 0, 1, 1], kd_tertile=[0, 1, 0, 1]))
    good = "a" * 64

    def build(self, base, hashes=None, directions=2, skip=None):
        for condition in ("DMSO", "Vem"):
            for fold, targets in ((0, ["T1", "T2"]), (1, ["T3", "T4"])):
                if skip == (condition, fold):
                    continue
                h = (hashes or {}).get((condition, fold), self.good)
                write_fold(base, condition, fold, targets, h, directions=directions if fold == 1 else 2)

    def test_passes_and_stacks(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(Path(d))
            out = ev.verify_predictions(d, self.good, self.split, NODES, MODELS)
        self.assertEqual(out["Vem"]["targets"], ["T1", "T2", "T3", "T4"])
        self.assertEqual(out["DMSO"]["scaled"].shape, (len(MODELS), 4, 2, len(NODES)))
        np.testing.assert_array_equal(out["DMSO"]["fold"], [0, 0, 1, 1])

    def test_hash_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(Path(d), hashes={("Vem", 1): "b" * 64})
            with self.assertRaisesRegex(ValueError, "config_sha256"):
                ev.verify_predictions(d, self.good, self.split, NODES, MODELS)

    def test_missing_direction_fails(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(Path(d), directions=1)
            with self.assertRaisesRegex(ValueError, "scaled has shape"):
                ev.verify_predictions(d, self.good, self.split, NODES, MODELS)

    def test_missing_fold_fails(self):
        with tempfile.TemporaryDirectory() as d:
            self.build(Path(d), skip=("DMSO", 1))
            with self.assertRaises(ValueError):
                ev.verify_predictions(d, self.good, self.split, NODES, MODELS)

    def test_text_hash_ignores_line_endings(self):
        self.assertEqual(ev.text_sha256(b"a\r\nb\n"), ev.text_sha256(b"a\nb\n"))


class FrozenRunConfig(unittest.TestCase):
    def setUp(self):
        self.run = json.loads((ROOT / "config" / "evaluation_run.json").read_text(encoding="utf-8"))
        self.endpoints = json.loads((ROOT / "config" / "evaluation_endpoints.json").read_text(encoding="utf-8"))

    def test_constants_match_the_endpoints(self):
        e = self.endpoints["primary"]
        self.assertTrue(self.run["frozen"])
        self.assertEqual(self.run["seed"], e["effect_size"]["ci"]["seed"])
        self.assertEqual(self.run["bootstrap_draws"], e["effect_size"]["ci"]["draws"])
        self.assertEqual(self.run["noise_ceiling_splits"], e["noise_ceiling"]["splits"])
        self.assertEqual(self.run["clip_abs_r"], e["averaging"]["clip_abs_r"])
        self.assertEqual(self.run["models"], MODELS)
        self.assertEqual(ev.file_sha256(ROOT / "config" / "evaluation_endpoints.json"),
                         self.run["inputs"]["endpoints_sha256"])
        self.assertEqual(ev.file_sha256(ROOT / "config" / "sem_prediction.json"),
                         self.run["inputs"]["prediction_config_sha256"])

    def test_rules(self):
        r = self.run["rules"]
        t = r["tertiary"]
        self.assertEqual(sorted(t["cell_level_models"] + t["mean_shift_models"]), sorted(MODELS))
        fam = r["holm"]["families"]
        self.assertEqual(len(fam["secondary_z_amended"]), 3)
        self.assertEqual(len(fam["secondary_rmse"]), 4)
        self.assertIn("passed explicitly", r["wilcoxon"]["method_rule"])
        self.assertIn("stop", r["verification"]["included_zero"])


if __name__ == "__main__":
    unittest.main()
