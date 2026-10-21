"""Held-out evaluation of the frozen SEM knockdown predictions (config/evaluation_run.json).

Prediction-side verification (hashes, coverage, constant predictions, constant shifts) runs before
any observed response is read. Observed changes use only evaluation cells (calib_fold 1 - c for the
prediction calibrated on calib_fold c). Statistics: Pearson r / Fisher z with a clip, RMSE, the
pre-registered and amended NaN rules, Wilcoxon signed-rank with an explicit method, Hodges-Lehmann
with a percentile bootstrap over targets, Holm, split-half noise ceilings, Fufa score contrasts and
distribution-level statistics.
"""

import hashlib
import json
import subprocess
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import dispersion as ds

CONSTANT_TOL = 1e-12
SHIFT_TOL = 1e-10
STRUCTURE_MODELS = {"cell_specific": "cell_specific_cscn_sem", "shared": "shared_cscn_sem",
                    "cluster": "cluster_level_sem", "local_association": "local_association_network_sem"}


# ---------------------------------------------------------------- hashes
def text_sha256(data):
    """SHA-256 of file bytes with CRLF replaced by LF (the rule used for the prediction config)."""
    return hashlib.sha256(data.replace(b"\r\n", b"\n")).hexdigest()


def file_sha256(path):
    return text_sha256(Path(path).read_bytes())


def git_file_sha256(root, commit, path):
    data = subprocess.run(["git", "show", f"{commit}:{path}"], cwd=root, capture_output=True, check=True).stdout
    return text_sha256(data)


def index_digest(array):
    """Same digest as scripts/predict_sem_knockdown.py records for the reference cells."""
    return hashlib.sha256(np.asarray(array, dtype=np.int64).tobytes()).hexdigest()


# ---------------------------------------------------------------- prediction verification
def verify_predictions(base, expected_hash, split, nodes, models, conditions=("DMSO", "Vem")):
    """Check every predictions.npz and its shards; return stacked arrays per condition.

    Returns condition -> dict(targets, fold, per_unit (M, T, P), scaled (M, T, 2, P), delta (T, 2),
    included (M, T), fufa_scaled (M, n_ref, T, 2), reference_cells, lam), targets sorted by name.
    Raises on any hash, linkage, shape or coverage mismatch.
    """
    base = Path(base)
    folds = sorted(split.fold.unique())
    out = {}
    for condition in conditions:
        parts, refs = [], None
        for fold in folds:
            d = base / condition / f"fold{fold}"
            shards = sorted(d.glob("shard_*_of_*"))
            if not shards:
                raise ValueError(f"{d}: no prediction shards")
            configs = [json.loads((s / "config.json").read_text(encoding="utf-8")) for s in shards]
            for s, c in zip(shards, configs):
                if c.get("config_sha256") != expected_hash:
                    raise ValueError(f"{s}: config_sha256 {c.get('config_sha256')} != {expected_hash}")
            common = [{k: v for k, v in c.items() if k != "shard"} for c in configs]
            if any(c != common[0] for c in common) or \
                    sorted(c["shard"] for c in configs) != list(range(common[0]["n_shards"])):
                raise ValueError(f"{d}: shards with different configurations or missing shards")
            config = common[0]
            if config["condition"] != condition or config["fold"] != fold:
                raise ValueError(f"{d}: shard config is for another condition or fold")
            held = sorted(split.target[split.fold.eq(fold)])
            with np.load(d / "predictions.npz") as z:
                arrays = {k: z[k] for k in z.files}
            targets = [str(t) for t in arrays["targets"]]
            if targets != held or config["targets"] != held:
                raise ValueError(f"{d}: held-out targets differ from config/target_split.json")
            if [str(m) for m in arrays["models"]] != list(models):
                raise ValueError(f"{d}: unexpected model list")
            if [str(n) for n in arrays["nodes"]] != list(nodes):
                raise ValueError(f"{d}: nodes differ from the frozen panel")
            if float(arrays["lam"]) != float(config["lam"]):
                raise ValueError(f"{d}: lambda differs from the shard config")
            if index_digest(arrays["reference_cells"]) != config["reference_sha256"]:
                raise ValueError(f"{d}: reference cells differ from the shard config")
            M, T, P, n_ref = len(models), len(held), len(nodes), len(arrays["reference_cells"])
            expected = dict(per_unit=(M, T, P), scaled=(M, T, 2, P), delta=(T, 2), included=(M, T),
                            fufa_scaled=(M, n_ref, T, 2))
            for key, shape in expected.items():
                if arrays[key].shape != shape:
                    raise ValueError(f"{d}: {key} has shape {arrays[key].shape}, expected {shape}")
            if refs is None:
                refs = arrays["reference_cells"]
            elif not np.array_equal(refs, arrays["reference_cells"]):
                raise ValueError(f"{condition}: reference cells differ between folds")
            parts.append((fold, targets, arrays))
        targets = [t for _, ts, _ in parts for t in ts]
        if sorted(targets) != sorted(split.target) or len(set(targets)) != len(targets):
            raise ValueError(f"{condition}: held-out targets do not cover the split exactly once")
        order = np.argsort(targets, kind="stable")
        cat = lambda key, axis: np.concatenate([a[key] for _, _, a in parts], axis=axis)  # noqa: E731
        out[condition] = dict(
            targets=[targets[i] for i in order],
            fold=np.concatenate([[f] * len(ts) for f, ts, _ in parts])[order],
            per_unit=cat("per_unit", 1)[:, order], scaled=cat("scaled", 1)[:, order],
            delta=cat("delta", 0)[order], included=cat("included", 1)[:, order],
            fufa_scaled=cat("fufa_scaled", 2)[:, :, order], reference_cells=refs,
            lam=float(parts[0][2]["lam"]))
    return out


def evaluation_mask(nodes, node):
    """Panel genes evaluated for a target: every node except the perturbed one."""
    mask = np.asarray([n != node for n in nodes])
    if mask.sum() != len(nodes) - 1:
        raise ValueError(f"Target node {node!r} must appear exactly once in the panel")
    return mask


def is_constant(vector, tol=CONSTANT_TOL):
    vector = np.asarray(vector, float)
    return bool(np.all(np.isfinite(vector)) and np.ptp(vector) <= tol)


def included_zero_rows(predictions, models):
    """(model, condition, target) rows without any included reference cell."""
    rows = []
    for condition, p in predictions.items():
        for i, m in enumerate(models):
            for t, target in enumerate(p["targets"]):
                if p["included"][i, t] == 0:
                    rows.append(dict(model=m, condition=condition, target=target))
    return pd.DataFrame(rows, columns=["model", "condition", "target"])


def constant_prediction_flags(predictions, models, nodes, target_node):
    """Per (model, condition, target): per-unit prediction constant over the evaluation genes."""
    rows = []
    for condition, p in predictions.items():
        for i, m in enumerate(models):
            for t, target in enumerate(p["targets"]):
                mask = evaluation_mask(nodes, target_node[target])
                rows.append(dict(model=m, condition=condition, target=target,
                                 constant=is_constant(p["per_unit"][i, t, mask])))
    return pd.DataFrame(rows)


def check_constant_flags(flags, structure_targets):
    """Constant-prediction flags must reproduce the structure report; returns counts per model."""
    expected = structure_targets.assign(model=structure_targets.graph_type.map(STRUCTURE_MODELS))
    if expected.model.isna().any():
        raise ValueError("Unknown graph_type in the structure report")
    merged = expected.merge(flags, on=["model", "condition", "target"], how="left", validate="one_to_one")
    if merged.constant.isna().any():
        raise ValueError("Structure report rows without a prediction")
    bad = merged[merged.constant.astype(bool) != merged.constant_prediction_expected.astype(bool)]
    if len(bad):
        raise ValueError("Constant predictions differ from the structure report:\n"
                         + bad[["model", "condition", "target", "constant", "constant_prediction_expected"]]
                         .to_string(index=False))
    return constant_counts(flags)


def constant_counts(flags):
    """Targets constant in DMSO, in Vem and in at least one condition, per model."""
    wide = flags.pivot_table(index=["model", "target"], columns="condition", values="constant", aggfunc="first")
    out = wide.groupby(level="model").sum().astype(int)
    out["any_condition"] = wide.any(axis=1).groupby(level="model").sum().astype(int)
    return out


def check_constant_shift(predictions, models, model_names, tol=SHIFT_TOL):
    """Models without cell-level output: the stored per-cell shift must be constant over included cells."""
    for condition, p in predictions.items():
        for name in model_names:
            i = list(models).index(name)
            for t, target in enumerate(p["targets"]):
                for c in range(2):
                    values = p["fufa_scaled"][i, :, t, c]
                    values = values[np.isfinite(values)]
                    if not len(values) or np.ptp(values) > tol:
                        raise ValueError(f"{name} {condition} {target} direction {c}: per-cell shift is not constant")


def check_fufa_pseudobulk(predictions, w, tol=1e-8):
    """Mean stored Fufa shift over included cells equals w . pseudobulk prediction."""
    for condition, p in predictions.items():
        with np.errstate(invalid="ignore"):
            mean_shift = np.nanmean(p["fufa_scaled"], axis=1)       # (M, T, 2)
        linear = p["scaled"] @ np.asarray(w, float)                  # (M, T, 2)
        if not np.allclose(mean_shift, linear, rtol=0, atol=tol):
            raise ValueError(f"{condition}: mean per-cell Fufa shift differs from w . pseudobulk")


# ---------------------------------------------------------------- observed responses
def evaluation_rows(meta, group, condition, direction):
    """Boolean rows of a group's evaluation cells for a prediction calibrated on calib_fold direction."""
    rows = (meta.condition.eq(condition) & meta.target.eq(group) & meta.calib_fold.eq(1 - direction)).to_numpy()
    if not rows.any():
        raise ValueError(f"No {group} cells of {condition} in calib_fold {1 - direction}")
    return rows


def observed_change(values, meta, target, condition, direction):
    """Target - NTC mean over the evaluation cells (calib_fold 1 - direction), all columns."""
    values = np.asarray(values)
    t = evaluation_rows(meta, target, condition, direction)
    n = evaluation_rows(meta, "NTC", condition, direction)
    return values[t].astype(np.float64).mean(axis=0) - values[n].astype(np.float64).mean(axis=0)


def observed_mean_contrast(scores, meta, target, condition, direction):
    """Target - NTC mean of a per-cell score over the evaluation cells (calib_fold 1 - direction)."""
    scores = np.asarray(scores, float)
    return float(scores[evaluation_rows(meta, target, condition, direction)].mean()
                 - scores[evaluation_rows(meta, "NTC", condition, direction)].mean())


# ---------------------------------------------------------------- per-target metrics
def pearson(x, y):
    """Pearson r; NaN if either vector is constant."""
    if is_constant(x) or is_constant(y):
        return np.nan
    return float(np.corrcoef(np.asarray(x, float), np.asarray(y, float))[0, 1])


def fisher_z(r, clip):
    r = np.asarray(r, float)
    return np.arctanh(np.clip(r, -clip, clip))


def rmse(predicted, observed):
    return float(np.sqrt(np.mean((np.asarray(predicted, float) - np.asarray(observed, float)) ** 2)))


def direction_metrics(predicted, observed, clip):
    """r, z, RMSE and constant flags for one (model, target, condition, direction)."""
    r = pearson(predicted, observed)
    return dict(r=r, z=float(fisher_z(r, clip)), rmse=rmse(predicted, observed),
                constant_prediction=is_constant(predicted), constant_observed=is_constant(observed))


def amended_z(z, constant_prediction, constant_observed):
    """Amended rule: constant prediction -> z = 0; a constant observed vector stays NaN."""
    z = np.asarray(z, float).copy()
    zero = np.asarray(constant_prediction, bool) & ~np.asarray(constant_observed, bool)
    z[zero] = 0.0
    return z


def _nanprop_mean(s):
    return float(np.mean(s.to_numpy(dtype=float)))


def combine_metrics(rows):
    """Add direction-combined and condition-pooled rows to per-direction metric rows.

    rows: model, target, condition, direction, r, z, rmse, constant_prediction, constant_observed.
    Adds z_amended / r_amended. z is averaged on the Fisher z scale with NaN propagation; RMSE by mean;
    constant flags by any. Combined rows have direction 'combined'; pooled rows condition 'pooled'.
    """
    rows = rows.copy()
    rows["direction"] = rows.direction.astype(object)
    rows["z_amended"] = amended_z(rows.z, rows.constant_prediction, rows.constant_observed)
    agg = dict(z=_nanprop_mean, z_amended=_nanprop_mean, rmse="mean", constant_prediction="any",
               constant_observed="any")
    combined = rows.groupby(["model", "target", "condition"], sort=False).agg(agg).reset_index()
    combined["direction"] = "combined"
    pooled = combined.groupby(["model", "target"], sort=False).agg(agg).reset_index()
    pooled["condition"], pooled["direction"] = "pooled", "combined"
    out = pd.concat([rows, combined, pooled], ignore_index=True)
    joint = out.direction.eq("combined").to_numpy()
    out.loc[joint, "r"] = np.tanh(out.loc[joint, "z"].to_numpy(float))
    out["r_amended"] = np.tanh(out.z_amended.to_numpy(float))
    return out


# ---------------------------------------------------------------- paired tests
def wilcoxon_method(d):
    """Frozen rule: exact iff n_eff <= 50 with no zero differences and no ties among |d|."""
    d = np.asarray(d, float)
    nonzero = d[d != 0]
    exact = len(nonzero) <= 50 and len(nonzero) == len(d) and len(np.unique(np.abs(d))) == len(d)
    return "exact" if exact else "approx"


def wilcoxon_test(d):
    """Two-sided Wilcoxon signed-rank on complete differences (zero_method 'wilcox')."""
    d = np.asarray(d, float)
    n_zero = int((d == 0).sum())
    n_eff = len(d) - n_zero
    method = wilcoxon_method(d)
    if n_eff == 0:
        return dict(n_zero_dropped=n_zero, n_effective=0, method=method, statistic=np.nan, p_value=np.nan)
    res = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided",
                         method="exact" if method == "exact" else "asymptotic")
    return dict(n_zero_dropped=n_zero, n_effective=n_eff, method=method, statistic=float(res.statistic),
                p_value=float(res.pvalue))


def hodges_lehmann(d):
    """Median of the Walsh averages (d_i + d_j) / 2, i <= j."""
    d = np.asarray(d, float)
    if not len(d):
        return np.nan
    i, j = np.triu_indices(len(d))
    return float(np.median((d[i] + d[j]) / 2))


def hodges_lehmann_ci(d, draws, seed, level=0.95, chunk=500):
    """Percentile bootstrap CI of Hodges-Lehmann over units (fresh generator per call)."""
    d = np.asarray(d, float)
    n = len(d)
    if not n:
        return np.nan, np.nan
    idx = np.random.default_rng(seed).integers(0, n, (draws, n))
    i, j = np.triu_indices(n)
    values = np.empty(draws)
    for s in range(0, draws, chunk):
        x = d[idx[s:s + chunk]]
        values[s:s + chunk] = np.median((x[:, i] + x[:, j]) / 2, axis=1)
    tail = (1 - level) / 2 * 100
    low, high = np.percentile(values, [tail, 100 - tail])
    return float(low), float(high)


def paired_test(model_values, reference_values, draws, seed, back_transform=False):
    """Paired Wilcoxon + Hodges-Lehmann on model - reference over complete pairs."""
    a, b = np.asarray(model_values, float), np.asarray(reference_values, float)
    complete = np.isfinite(a) & np.isfinite(b)
    d = a[complete] - b[complete]
    out = dict(n_targets=len(a), n_nan_dropped=int((~complete).sum()), **wilcoxon_test(d))
    hl = hodges_lehmann(d)
    low, high = hodges_lehmann_ci(d, draws, seed)
    out.update(hl=hl, hl_ci_low=low, hl_ci_high=high,
               median_model=float(np.median(a[complete])) if complete.any() else np.nan,
               median_reference=float(np.median(b[complete])) if complete.any() else np.nan)
    if back_transform:
        out.update(hl_back=float(np.tanh(hl)), hl_ci_low_back=float(np.tanh(low)),
                   hl_ci_high_back=float(np.tanh(high)))
    return out


def holm(p):
    """Holm step-down adjusted p-values over the finite entries; NaN stays NaN."""
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    m = len(finite)
    if not m:
        return out
    order = finite[np.argsort(p[finite], kind="stable")]
    adjusted = np.minimum(1.0, np.maximum.accumulate((m - np.arange(m)) * p[order]))
    out[order] = adjusted
    return out


# ---------------------------------------------------------------- noise ceiling
def spearman_brown(r):
    r = np.asarray(r, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * r / (1 + r)


def noise_ceiling_rng(seed, condition, target, direction):
    key = zlib.crc32(f"noise_ceiling|{condition}|{target}|{direction}".encode())
    return np.random.default_rng([seed, key])


def split_half_ceiling(target_values, ntc_values, splits, rng, clip):
    """Per split: Spearman-Brown corrected r between half-sample target - NTC change vectors (clipped)."""
    T, N = np.asarray(target_values, np.float64), np.asarray(ntc_values, np.float64)
    n_t, n_n = len(T), len(N)
    if n_t < 2 or n_n < 2:
        raise ValueError("A split half needs at least one cell of each group")
    out = np.empty(splits)
    for s in range(splits):
        pt, pn = rng.permutation(n_t), rng.permutation(n_n)
        at, bt, an, bn = pt[:n_t // 2], pt[n_t // 2:], pn[:n_n // 2], pn[n_n // 2:]
        va = T[at].mean(axis=0) - N[an].mean(axis=0)
        vb = T[bt].mean(axis=0) - N[bn].mean(axis=0)
        r = pearson(va, vb)
        out[s] = np.nan if np.isnan(r) else float(np.clip(spearman_brown(r), -clip, clip))
    return out


def combine_directions_z(r_values, clip):
    """tanh(mean of atanh(clipped r)) over directions (NaN propagates)."""
    return float(np.tanh(np.mean(fisher_z(r_values, clip))))


# ---------------------------------------------------------------- Fufa secondary endpoint
def fufa_contrast_rows(predicted, observed):
    """KD - NTC per condition and the Vem - DMSO interaction, per direction and combined.

    predicted, observed: condition -> length-2 array (direction c). Returns rows with
    contrast, direction, predicted, observed, abs_error; the combined row's abs_error is the mean
    of the two direction errors (predicted / observed: mean over directions, descriptive).
    """
    contrasts = {"KD_minus_NTC_DMSO": (np.asarray(predicted["DMSO"], float), np.asarray(observed["DMSO"], float)),
                 "KD_minus_NTC_Vem": (np.asarray(predicted["Vem"], float), np.asarray(observed["Vem"], float))}
    contrasts["interaction_Vem_minus_DMSO"] = (contrasts["KD_minus_NTC_Vem"][0] - contrasts["KD_minus_NTC_DMSO"][0],
                                               contrasts["KD_minus_NTC_Vem"][1] - contrasts["KD_minus_NTC_DMSO"][1])
    rows = []
    for name, (p, o) in contrasts.items():
        err = np.abs(p - o)
        for c in range(2):
            rows.append(dict(contrast=name, direction=c, predicted=p[c], observed=o[c], abs_error=err[c]))
        rows.append(dict(contrast=name, direction="combined", predicted=p.mean(), observed=o.mean(),
                         abs_error=err.mean()))
    return rows


# ---------------------------------------------------------------- tertiary endpoints
def distribution_stats(x, y):
    """Aim 1 estimators of x against reference y."""
    return dict(log_var_ratio=ds.log_var_ratio(x, y), iqr_diff=ds.iqr_diff(x, y),
                centered_wasserstein=ds.centered_wasserstein1(x, y))


def predicted_distribution_stats(ntc_scores, shifts):
    """Statistics of (included NTC + own shift) against the same included NTC cells unshifted."""
    ntc_scores, shifts = np.asarray(ntc_scores, float), np.asarray(shifts, float)
    keep = np.isfinite(shifts)
    if keep.sum() < 2:
        raise ValueError("Fewer than two evaluation NTC cells with a defined shift")
    base = ntc_scores[keep]
    return distribution_stats(base + shifts[keep], base), int((~keep).sum())


def observed_lvr_interval(target_scores, ntc_scores, draws, seed, condition, target, eval_fold, ntc_draws=None,
                          level=0.95):
    """Percentile bootstrap interval of the observed log variance ratio on evaluation cells."""
    t = ds.GroupDraws(target_scores, draws, seed, condition, f"{target}|eval{eval_fold}")
    r = ntc_draws if ntc_draws is not None else ds.GroupDraws(ntc_scores, draws, seed, condition,
                                                              f"NTC|eval{eval_fold}")
    lvr = np.log(t.var / r.var)
    tail = (1 - level) / 2 * 100
    low, high = np.percentile(lvr, [tail, 100 - tail])
    return float(low), float(high)


def calibration_flags(predicted, low, high):
    return dict(over_dispersion=bool(predicted > high), covered=bool(low <= predicted <= high))
